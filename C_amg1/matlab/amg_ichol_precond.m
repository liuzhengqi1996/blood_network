function x = amg_ichol_precond(amg, A, L, g)
% (1) u1 = backslash AMG cycle on g
% (2) u2 = u1 + L\(L'\(g - A*u1))
% (3) x  = u2 + fwdslash AMG cycle on (g - A*u2)
n = size(A,1);
u1 = zeros(n,1);
u1 = vcycle_backslash(amg, 1, g, u1);

r1 = g - A*u1;
ic = L.' \ (L \ r1);
u2 = u1 + ic;

r2 = g - A*u2;
corr = zeros(n,1);
corr = vcycle_fwdslash(amg, 1, r2, corr);
x = u2 + corr;
end

function x = vcycle_backslash(amg, lev, b, x)
nu = amg.nu;
if lev == amg.nlev
    [x,~] = pcg(amg.A{lev}, b, 1e-12, size(amg.A{lev},1), amg.Lc, amg.Lc.');
    return;
end
x = x + (b - amg.A{lev}*x) ./ amg.l1{lev};
for s = 2:nu
    x = x + (b - amg.A{lev}*x) ./ amg.l1{lev};
end
r = b - amg.A{lev}*x;
rc = amg.P{lev}.' * r;
ec = zeros(size(rc));
ec = vcycle_backslash(amg, lev+1, rc, ec);
x = x + amg.P{lev}*ec;
end

function x = vcycle_fwdslash(amg, lev, b, x)
nu = amg.nu;
if lev == amg.nlev
    [x,~] = pcg(amg.A{lev}, b, 1e-12, size(amg.A{lev},1), amg.Lc, amg.Lc.');
    return;
end
r = b - amg.A{lev}*x;
rc = amg.P{lev}.' * r;
ec = zeros(size(rc));
ec = vcycle_fwdslash(amg, lev+1, rc, ec);
x = x + amg.P{lev}*ec;
for s = 1:nu
    x = x + (b - amg.A{lev}*x) ./ amg.l1{lev};
end
end
