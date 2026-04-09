function [cf, nc] = rs_coarsening(S, lam_max)
% RS_COARSENING  Classical Ruge-Stuben C/F splitting.
%
%   [cf, nc] = rs_coarsening(S, lam_max)
%
%   First pass: greedy independent-set selection using bucket sort.
%   Second pass: promote any F-point with no strong C-neighbor to C.
%
%   Inputs:
%     S       : sparse n x n strength-of-connection matrix
%     lam_max : scalar, max column count of A (from rs_strength)
%
%   Output:
%     cf(i) =  1 for C-points, -1 for F-points (length n)
%     nc    = number of C-points
%
%   Requires compiled MEX: mex rs_coarsening_mex.c

[cf, nc] = rs_coarsening_mex(S, lam_max);
end
