function P = rs_interpolation(A, S, cf)
% RS_INTERPOLATION  Standard (direct) interpolation for classical RS AMG.
%
%   P = rs_interpolation(A, S, cf)
%
%   C-points: injection (P(i, cmap(i)) = 1)
%   F-points: standard interpolation weights
%
%   Returns sparse P of size n x nc.
%
%   Requires compiled MEX: mex rs_interpolation_mex.c

P = rs_interpolation_mex(A, S, cf);
end
