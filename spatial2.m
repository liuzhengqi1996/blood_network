clc; clear; close all;


filename = "CoW.cs31";
%filename = "S58_bb7";
[faceMx, ptCoordMx, dia, BC, np, nf, nt] = caseReaderMJ(filename);
[alpha, l, d4] = Resistance(ptCoordMx, faceMx, dia, nf);
C1 = ConnectivityMx(nf, np, faceMx);
boundary_nodes = BC(:,1);
all_nodes = (1:np)';
interior_nodes = setdiff(all_nodes, boundary_nodes);
C1_I = C1(:, interior_nodes);
C1_B = C1(:, boundary_nodes);
C1_I_T = C1_I';
C1_B_T = C1_B';
p_B = sparse(np, 1);
p_B(boundary_nodes) = BC(:,3);
p_B = p_B(boundary_nodes);
alpha_inv = spdiags(1 ./ diag(alpha), 0, nf, nf);
RHS = -C1_I_T * alpha_inv * C1_B * p_B;
LHS = C1_I_T * alpha_inv * C1_I;
tol = 1e-3; max_iter = 1000;
p_I = GaussSeidel(LHS, RHS, tol, max_iter);
f = alpha_inv * (C1_I * p_I + C1_B * p_B);
p = sparse(np,1); p(interior_nodes) = p_I; p(boundary_nodes) = p_B;

%% === Supernode  ===
maxSize = 5;
supernodeLabels = buildSupernodesFromLeaves(faceMx, np, maxSize);
[faceMx_coarse, AdjC, C1_coarse] = buildCoarseGraph(faceMx, supernodeLabels);
nc = size(faceMx_coarse,1);
alpha_coarse_diag = zeros(nc,1);
edgeMap = containers.Map();
for i = 1:size(faceMx,1)
    v1 = faceMx(i,2); v2 = faceMx(i,3);
    s1 = supernodeLabels(v1); s2 = supernodeLabels(v2);
    if s1 ~= s2
        key = sprintf('%d-%d', min(s1,s2), max(s1,s2));
        if ~isKey(edgeMap, key)
            edgeMap(key) = alpha(i,i);
        else
            edgeMap(key) = edgeMap(key) + alpha(i,i);
        end
    end
end
for i = 1:nc
    s1 = faceMx_coarse(i,2);
    s2 = faceMx_coarse(i,3);
    key = sprintf('%d-%d', min(s1,s2), max(s1,s2));
    alpha_coarse_diag(i) = edgeMap(key);
end
alpha_coarse = sparse(diag(alpha_coarse_diag));
nsuper = max(supernodeLabels);
p_b = zeros(nsuper,1);
D_coarse = sparse(nsuper, nsuper);
for i = 1:size(BC,1)
    node = BC(i,1);
    super = supernodeLabels(node);
    D_coarse(super,super) = 1;
    p_b(super) = p_b(super) + BC(i,3);
end
for i = 1:nsuper
    if D_coarse(i,i) > 0
        p_b(i) = p_b(i) / D_coarse(i,i);
    end
end
C1T = C1_coarse'; D2 = speye(nsuper) - D_coarse;
rhs = [zeros(nc,1); D_coarse * p_b];
M = [alpha_coarse, -C1_coarse; D2 * C1T, D_coarse];
xx = M \ rhs;
f_coarse = xx(1:nc);
p_coarse = xx(nc+1:end);

%% === Visualizations ===
save('results.mat', 'p', 'f');
figure;
RenderTubesAsGraph2024(ptCoordMx, faceMx, dia, p, [], jet(255), 'Pressure Distribution');
figure;
RenderTubesAsGraph2024(ptCoordMx, faceMx, dia, abs(f), [], hot(255), 'Flow Magnitude');
superColors = lines(nsuper);
nodeColors = superColors(supernodeLabels,:);
figure; scatter3(ptCoordMx(:,1), ptCoordMx(:,2), ptCoordMx(:,3), 60, nodeColors, 'filled');
title('Node Grouping by Supernodes'); axis equal;
superCoords = zeros(nsuper, 3);
for i = 1:nsuper
    superCoords(i,:) = mean(ptCoordMx(supernodeLabels == i,:), 1);
end
figure;
G = graph(AdjC);
plot(G, 'XData', superCoords(:,1), 'YData', superCoords(:,2), 'ZData', superCoords(:,3), ...
    'NodeLabel', {}, 'EdgeAlpha', 0.8, 'LineWidth', 2);
title('Coarse Graph (Supernodes)'); axis equal;

%% === Residual ===
residual = norm(C1 * p - alpha * f);
fprintf('Residual Norm: %.10f\n', residual);

%% === Functions ===

function [faceMx, ptCoordMx, dia, BC, np, nf, nt] = caseReaderMJ(filename)
    faceMx = load(strcat(filename, '.fMx'));
    ptCoordMx = load(strcat(filename, '.pMx'));
    dia = load(strcat(filename, '.dia'));
    BC = load(strcat(filename, '.BC'));
    np = length(ptCoordMx);
    nf = size(faceMx,1);
    nt = np + nf;
end

function [alpha, l, d4] = Resistance(ptCoordMx, faceMx, dia, nf)
    mu = 5.3317E-7;
    f = 128 * mu / pi;
    l = vecnorm(ptCoordMx(faceMx(:,2),:) - ptCoordMx(faceMx(:,3),:), 2, 2);
    d4 = dia.^4;
    r = f * (l ./ d4) * 1000;
    alpha = sparse(diag(r));
end

function C1 = ConnectivityMx(nf, np, faceMX)
    C1 = sparse(nf, np);
    for i = 1:nf
        C1(i, faceMX(i,2)) = 1;
        C1(i, faceMX(i,3)) = -1;
    end
end

function x = GaussSeidel(A, b, tol, max_iter)
    n = length(b);
    x = zeros(n,1);
    for iter = 1:max_iter
        x_old = x;
        for i = 1:n
            sum1 = A(i,1:i-1) * x(1:i-1);
            sum2 = A(i,i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) / A(i,i);
        end
        if norm(x - x_old, inf) < tol
            fprintf('Gauss-Seidel converged in %d iterations.\n', iter);
            return;
        end
    end
    warning('Gauss-Seidel did not converge within the maximum iterations.');
end

function [p,G] = RenderTubesAsGraph2024(ptCoordMx, faceMx, dia, faceProp, pointProp, cmap, figureTitle)
    colormap = cmap;
    nf = length(dia);
    G = digraph(faceMx(:,2), faceMx(:,3), 1:nf);
    p = plot(G, 'XData', ptCoordMx(:,1), 'YData', ptCoordMx(:,2), 'ZData', ptCoordMx(:,3));
    p.MarkerSize = 4;
    p.LineWidth = 6;
    p.ArrowPosition = 1;
    p.EdgeAlpha = 1;
    if length(cmap(:,1)) == 1
        p.EdgeColor = cmap(1,:);
    else
        if size(faceProp) ~= 0
            p.EdgeCData = faceProp(G.Edges.Weight);
            colorbar
        end
    end
    title(figureTitle, 'interpreter', 'none');
    colorbar; axis equal; view(2);
end

function supernodeLabels = buildSupernodesFromLeaves(faceMx, np, maxSize)
    adj = cell(np, 1);
    for i = 1:size(faceMx,1)
        u = faceMx(i,2); v = faceMx(i,3);
        adj{u}(end+1) = v;
        adj{v}(end+1) = u;
    end
    deg = cellfun(@length, adj);
    leaves = find(deg == 1);
    supernodeLabels = zeros(np,1);
    visited = false(np,1);
    label = 1;
    for i = 1:length(leaves)
        if visited(leaves(i)), continue; end
        queue = leaves(i);
        group = [];
        while ~isempty(queue) && length(group) < maxSize
            n = queue(1); queue(1) = [];
            if visited(n), continue; end
            visited(n) = true;
            group(end+1) = n;
            for neighbor = adj{n}
                if ~visited(neighbor)
                    queue(end+1) = neighbor;
                end
            end
        end
        supernodeLabels(group) = label;
        label = label + 1;
    end
    for i = 1:np
        if ~visited(i)
            supernodeLabels(i) = label;
            label = label + 1;
        end
    end
end

function [faceMx_coarse, AdjC, C1_coarse] = buildCoarseGraph(faceMx, supernodeLabels)
    nf = size(faceMx,1);
    superEdges = [];
    superEdgeMap = containers.Map();
    index = 1;
    for i = 1:nf
        v1 = faceMx(i,2);
        v2 = faceMx(i,3);
        s1 = supernodeLabels(v1);
        s2 = supernodeLabels(v2);
        if s1 ~= s2
            key = sprintf('%d-%d', min(s1,s2), max(s1,s2));
            if ~isKey(superEdgeMap, key)
                superEdgeMap(key) = index;
                superEdges(index,:) = [s1, s2]; %#ok<AGROW>
                index = index + 1;
            end
        end
    end
    nc = size(superEdges,1);
    faceMx_coarse = [(1:nc)', superEdges];
    nsuper = max(supernodeLabels);
    AdjC = sparse(superEdges(:,1), superEdges(:,2), 1, nsuper, nsuper);
    AdjC = AdjC + AdjC';
    C1_coarse = sparse(nc, nsuper);
    for i = 1:nc
        s1 = superEdges(i,1);
        s2 = superEdges(i,2);
        if s1 < s2
            C1_coarse(i,s1) = 1;
            C1_coarse(i,s2) = -1;
        else
            C1_coarse(i,s1) = -1;
            C1_coarse(i,s2) = 1;
        end
    end
end
