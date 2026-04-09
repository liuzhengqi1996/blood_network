function x = cf_gauss_seidel(A, b, x, cf_order, nu, direction)
% CF_GAUSS_SEIDEL  CF-ordered Gauss-Seidel smoother.
%
%   x = cf_gauss_seidel(A, b, x, cf_order, nu, direction)
%
%   direction: 'forward' or 'backward'
%   cf_order: permutation vector (C-points first, then F-points)
%   nu: number of sweeps
%
%   Requires compiled MEX: mex cf_gauss_seidel_mex.c

if strcmp(direction, 'forward')
    x = cf_gauss_seidel_mex(A, b, x, cf_order, nu, 0);
else
    x = cf_gauss_seidel_mex(A, b, x, cf_order, nu, 1);
end
end
