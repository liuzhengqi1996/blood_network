% compile_mex.m - Compile MEX files for AMG solver
%
% Usage: >> compile_mex

fprintf('Compiling rs_strength_mex ...\n');
mex rs_strength_mex.c

fprintf('Compiling rs_coarsening_mex ...\n');
mex rs_coarsening_mex.c

fprintf('Compiling rs_interpolation_mex ...\n');
mex rs_interpolation_mex.c

fprintf('Compiling cf_gauss_seidel_mex ...\n');
mex cf_gauss_seidel_mex.c

fprintf('All MEX files compiled.\n');
