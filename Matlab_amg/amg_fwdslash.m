function x = amg_fwdslash(amg, lev, b, x)
% AMG_FWDSLASH  Fwdslash half-cycle (recursive).
%
%   x = amg_fwdslash(amg, lev, b, x)
%
%   Coarsest level: direct solve (A \ b)
%   Otherwise: residual -> restrict -> recurse -> prolongate + correct
%              -> post-smooth (backward GS)

nlev = amg.num_levels;
L = amg.levels{lev};
A = L.A;

if lev == nlev
    % Coarsest level: direct solve
    x = A \ b;
    return;
end

nu = amg.nu;

% Residual: r = b - A*x
r = b - A * x;

% Restrict: rc = P' * r
P = L.P;
rc = P' * r;

% Recurse on coarse grid
nc = size(rc, 1);
ec = zeros(nc, 1);
ec = amg_fwdslash(amg, lev + 1, rc, ec);

% Prolongate and correct: x = x + P * ec
x = x + P * ec;

% Post-smoothing: backward CF-ordered GS
x = cf_gauss_seidel(A, b, x, L.cf_order, nu, 'backward');
end
