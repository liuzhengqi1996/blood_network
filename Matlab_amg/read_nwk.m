
function [LHS,RHS,interior_nodes]=read_nwk(filename)

% filename = "C:\Users\weams\OneDrive\reports and talks\misc memos and technical reports\Zikatanov collab\nwk_Cartesian_100x100x10";
tic;
[faceMx, ptCoordMx, dia, BC, np, nf, nt] = caseReaderMJ(filename);
tload=toc;

fprintf('Loaded Network Data: %.4fs\n',tload);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Compute Resistance Matrix A
% [alpha, l, d4] = Resistance(ptCoordMx, faceMx, dia, nf);
tic;
% 3. Generate Connectivity Matrix C1
C1 = ConnectivityMx(nf, np, faceMx);

% 4. Split Nodes into Interior and Boundary Nodes
boundary_nodes = BC(:,1); % Boundary node indices
all_nodes = (1:np)'; % All node indices
interior_nodes = setdiff(all_nodes, boundary_nodes); % Interior node indices

% 5. Extract Submatrices
C1_I = C1(:, interior_nodes); % Connectivity matrix for interior nodes
C1_B = C1(:, boundary_nodes); % Connectivity matrix for boundary nodes
C1_I_T = C1_I'; % Transpose of C1_I
C1_B_T = C1_B';
% 6. Extract Boundary Pressures
p_B = sparse(np, 1); % Initialize all pressures
p_B(boundary_nodes) = BC(:,3); % Assign boundary node pressures
p_B = p_B(boundary_nodes); % Keep only boundary pressures

% 7. Compute Right-Hand Side (RHS)
% alpha_inv = spdiags(1 ./ diag(alpha), 0, nf, nf); % Compute alpha inverse efficiently
alpha_inv = makeDiffusiveConductance(faceMx,ptCoordMx);
RHS = -C1_I_T * alpha_inv * C1_B * p_B;

LHS = C1_I_T * alpha_inv * C1_I;


%% uses pcg with ichol
tmatrix=toc;
fprintf('Formed the matrix: %.2fs\n',tmatrix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% leave all commented below for record purposes.
% [x, iter, resvec] = solver_spd_ilu0(LHS, RHS, 1e-10, 1000,1e-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function W = makeDiffusiveConductance(faceMx,ptCoordMx)
%     nf = size(faceMx,1);
%     faces_vec = ptCoordMx(faceMx(:,3),:) - ptCoordMx(faceMx(:,2),:);
%     faces_x = find(faces_vec(:,1)>0); hx = vecnorm(faces_vec(faces_x(1),:),2,2);
%     faces_y = find(faces_vec(:,2)>0); hy = vecnorm(faces_vec(faces_y(1),:),2,2);
%     faces_z = find(faces_vec(:,3)>0); hz = vecnorm(faces_vec(faces_z(1),:),2,2);
%     conductance = nan(nf,1);
%     conductance(faces_x) = hy*hz/hx;
%     conductance(faces_y) = hx*hz/hy;
%     conductance(faces_z) = hx*hy/hz;
%     W = spdiags(conductance,0,nf,nf);
% end

% function [faceMx,ptCoordMx,dia,BC,np,nf,nt]=caseReaderMJ(filename)
% myfileName = strcat(filename,'.fMx');
% faceMx = load(myfileName);
% myfileName = strcat(filename,'.pMx');
% ptCoordMx = load(myfileName);
% myfileName = strcat(filename,'.dia');
% dia = load(myfileName);
% myfileName = strcat(filename,'.BC');
% BC = load(myfileName);
% np= length(ptCoordMx); nf= length(faceMx(:,2)); nt= np+nf;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Function to Compute Resistance Matrix
% function [alpha, l, d4] = Resistance(ptCoordMx, faceMx, dia, nf)
%     mu = 5.3317E-7; 
%     f = 128 * mu / pi;
%     l = vecnorm((ptCoordMx(faceMx(:,2),:) - ptCoordMx(faceMx(:,3),:)), 2, 2);
%     d4 = dia.^4;
%     r = f * (l ./ d4) * 1000; % Convert to ml
%     alpha = sparse(diag(r));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Function to Compute Connectivity Matrix
% % function C1 = ConnectivityMx(nf, np, faceMX)
% %     C1 = sparse(nf, np);
% %     for i = 1:nf
% %         C1(i, faceMX(i,2)) = 1; 
% %         C1(i, faceMX(i,3)) = -1;
% %     end
% % end
% function C1 = ConnectivityMx(nf,np,faceMx)
%     fi = (1:nf);
%     p1 = faceMx(:,2)'; p2 = faceMx(:,3)'; % P1IDx   % P2 Idx
%     C1=sparse([fi, fi], [p1, p2],[ones(1,nf), (-1)*ones(1,nf)],nf,np);
% end
