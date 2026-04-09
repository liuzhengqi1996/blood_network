function z = amg_precondition(amg, g)
% AMG_PRECONDITION  Combined AMG + ichol preconditioner (one application).
%
%   z = amg_precondition(amg, g)
%
%   Three-stage preconditioner:
%     (1) u1 = backslash(g)                         AMG backslash cycle
%     (2) u2 = u1 + L \ (L' \ (g - A*u1))          ichol correction
%     (3) z  = u2 + fwdslash(g - A*u2)              AMG fwdslash cycle

A = amg.A;
n = size(A, 1);

% Step 1: u1 = backslash(g)
u1 = amg_backslash(amg, 1, g, zeros(n, 1));

% Step 2: ichol correction
r1 = g - A * u1;
u2 = u1 + amg.L \ (amg.L' \ r1);

% Step 3: fwdslash correction
r2 = g - A * u2;
z = u2 + amg_fwdslash(amg, 1, r2, zeros(n, 1));
end
