function x = amg_backslash(amg, lev, b, x)
% AMG_BACKSLASH  Backslash half-cycle (recursive).
%
%   x = amg_backslash(amg, lev, b, x)
%
%   Coarsest level: direct solve (A \ b)
%   Otherwise: pre-smooth (forward GS) -> residual -> restrict -> recurse
%              -> prolongate + correct

nlev = amg.num_levels;
L = amg.levels{lev};
A = L.A;

if lev == nlev
    % Coarsest level: direct solve
    x = A \ b;
    return;
end

nu = amg.nu;

% Pre-smoothing: forward CF-ordered GS
x = cf_gauss_seidel(A, b, x, L.cf_order, nu, 'forward');

% Residual: r = b - A*x
r = b - A * x;

% Restrict: rc = P' * r
P = L.P;
rc = P' * r;

% Recurse on coarse grid
nc = size(rc, 1);
ec = zeros(nc, 1);
ec = amg_backslash(amg, lev + 1, rc, ec);

% Prolongate and correct: x = x + P * ec
x = x + P * ec;
end
