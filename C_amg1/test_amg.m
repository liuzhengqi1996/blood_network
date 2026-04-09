% test_amg.m - AMG+ichol PCG on problem_case.mat (pure MATLAB port of C_amg1)

addpath(fullfile(fileparts(mfilename('fullpath')), 'matlab'));

matfile = '/Users/david/Desktop/agg_pcg_ba1x/blood_networks_new/problem_case.mat';

fprintf('Loading %s ...\n', matfile);
t0 = tic;
S = load(matfile);
A = S.A; b = S.b; clear S;
if ~issparse(A), A = sparse(A); end
A = (A + A.') * 0.5;
n = size(A,1);
fprintf('  n=%d nnz=%d ||b||=%.3e (%.2fs)\n', n, nnz(A), norm(b), toc(t0));

opts.threshold  = 0.25;
opts.min_size   = 512;
opts.max_levels = 25;
opts.nu         = 2;
tol     = 1e-6;
maxit   = 500;

fprintf('\nBuilding AMG hierarchy ...\n');
t0 = tic;
amg = amg_setup(A, opts);
t_setup_amg = toc(t0);
fprintf('  AMG setup: %.2fs\n', t_setup_amg);

fprintf('Computing ichol(0) of A ...\n');
t0 = tic;
try
    L = ichol(A, struct('type','nofill','diagcomp',1e-6));
catch
    L = ichol(A, struct('type','nofill','diagcomp',1e-3));
end
t_ichol = toc(t0);
fprintf('  ichol: nnz(L)=%d (%.2fs)\n', nnz(L), t_ichol);

fprintf('\nPCG with AMG+ichol (tol=%.0e maxit=%d) ...\n', tol, maxit);
M = @(g) amg_ichol_precond(amg, A, L, g);
t0 = tic;
[x, flag, relres, iter] = pcg(A, b, tol, maxit, M);
t_solve = toc(t0);
fprintf('  flag=%d iter=%d relres=%.3e time=%.2fs\n', flag, iter, relres, t_solve);

fprintf('\n  %-20s %6s %12s %12s %12s\n','Solver','Iter','Setup(s)','Solve(s)','Relres');
fprintf('  %-20s %6d %12.2f %12.2f %12.2e\n','AMG+ichol',iter, ...
        t_setup_amg+t_ichol, t_solve, relres);
