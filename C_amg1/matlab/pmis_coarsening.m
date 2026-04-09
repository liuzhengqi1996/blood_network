function cf = pmis_coarsening(S)
% PMIS coarsening on strength matrix S.
% Returns cf: +1 = C-point, -1 = F-point.
n = size(S,1);

% Symmetric strength graph G = (S | S') with no diag
Gp = spones(S) | spones(S');
Gp = Gp - spdiags(diag(Gp),0,n,n);

% lambda_i = number of points i influences = nnz of column i of S = nnz row i of S'
lambda = full(sum(spones(S), 1)).';

rng(12345);
w = lambda + rand(n,1);

cf = zeros(n,1);          % 0 = undecided
[I,J] = find(Gp);

while true
    u = (cf == 0);
    if ~any(u), break; end

    % Find undecided edges where neighbor has strictly larger weight
    valid = u(I) & u(J) & (w(J) > w(I));
    not_max = false(n,1);
    not_max(I(valid)) = true;

    new_C = u & ~not_max;
    if ~any(new_C)
        % safety: mark all remaining undecided as C
        cf(u) = 1;
        break;
    end
    cf(new_C) = 1;

    % Mark F: undecided nodes with at least one C neighbor in G
    c_mask = double(cf == 1);
    has_c = (Gp * c_mask) > 0;
    new_F = (cf == 0) & has_c;
    cf(new_F) = -1;
end

% Second pass: every F-point must have a strong C-neighbor (in S, not G)
Sp = spones(S);
c_mask = double(cf == 1);
nC_per_row = Sp * c_mask;
promote = (cf == -1) & (nC_per_row == 0);
cf(promote) = 1;
end
