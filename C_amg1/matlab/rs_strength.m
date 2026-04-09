function S = rs_strength(A, theta)
% RS strength of connection: -a_ij >= theta * max_{k!=i}(-a_ik)
n = size(A,1);
Aoff = A - spdiags(spdiags(A,0), 0, n, n);
[I,J,V] = find(Aoff);
negV = -V;
rmax = accumarray(I, negV, [n 1], @max, 0);
thr  = theta .* rmax;
mask = negV >= thr(I) & negV > 0;
S = sparse(I(mask), J(mask), V(mask), n, n);
end
