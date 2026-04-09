function [S, lam_max] = rs_strength(A, theta)
% RS_STRENGTH  Strength of connection (Ruge-Stuben criterion).
%
%   [S, lam_max] = rs_strength(A, theta)
%
%   For each row i, off-diagonal j is "strong" if
%       -a_{ij} >= theta * max_{k ~= i}( -a_{ik} )
%
%   Returns sparse S with the same dimensions as A, containing only
%   strong off-diagonal entries with their original values from A.
%   lam_max is the max column count of A (for bucket allocation).
%
%   Requires compiled MEX: mex rs_strength_mex.c

[S, lam_max] = rs_strength_mex(A, theta);
end
