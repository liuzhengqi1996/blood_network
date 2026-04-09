function x = amg_vcycle(amg, b)
% Symmetric V-cycle: backslash + fwdslash, l1-Jacobi smoother.
n = size(amg.A{1},1);
x = zeros(n,1);
x = backslash(amg, 1, b, x);
r = b - amg.A{1}*x;
corr = zeros(n,1);
corr = fwdslash(amg, 1, r, corr);
x = x + corr;
end

function x = backslash(amg, lev, b, x)
nu = amg.nu;
if lev == amg.nlev
    x = coarse_solve(amg, b);
    return;
end
x = l1_smooth(amg.A{lev}, amg.l1{lev}, b, x, nu);
r = b - amg.A{lev}*x;
rc = amg.P{lev}.' * r;
ec = zeros(size(rc));
ec = backslash(amg, lev+1, rc, ec);
x = x + amg.P{lev}*ec;
end

function x = fwdslash(amg, lev, b, x)
nu = amg.nu;
if lev == amg.nlev
    x = coarse_solve(amg, b);
    return;
end
r = b - amg.A{lev}*x;
rc = amg.P{lev}.' * r;
ec = zeros(size(rc));
ec = fwdslash(amg, lev+1, rc, ec);
x = x + amg.P{lev}*ec;
x = l1_smooth(amg.A{lev}, amg.l1{lev}, b, x, nu);
end

function x = l1_smooth(A, l1d, b, x, nu)
for s = 1:nu
    x = x + (b - A*x) ./ l1d;
end
end

function x = coarse_solve(amg, b)
A = amg.A{amg.nlev}; L = amg.Lc;
[x,~] = pcg(A, b, 1e-12, size(A,1), L, L.');
end
