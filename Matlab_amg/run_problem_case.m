% run_problem_case.m  — solve problem_case.mat with AMG+ichol PCG
cd('/Users/david/Desktop/agg_pcg_ba1x/blood_networks_new/C_amg_release/src');

fprintf('Loading problem_case.mat ...\n');
t0 = tic;
S = load('/Users/david/Desktop/agg_pcg_ba1x/blood_networks_new/problem_case.mat');
A = S.A; b = S.b; clear S;
fprintf('  loaded in %.1fs:  n=%d  nnz(A)=%d  ||b||=%.3e\n', ...
        toc(t0), size(A,1), nnz(A), norm(b));

if ~issparse(A)
    fprintf('A is dense — converting to sparse ...\n');
    A = sparse(A);
    fprintf('  nnz(A)=%d\n', nnz(A));
end

% symmetry check (cheap on a few entries)
sym_err = nnz(A - A') / max(1,nnz(A));
fprintf('  rel sym mismatch = %.2e\n', sym_err);



opts.threshold  =   0.25; %0.50;  
opts.min_size   =   512;%50000;  %
opts.max_levels =         25; %8%
opts.nu         = 2;
tol     = 1e-6;
maxiter = 500;

amg = [];
t_setup = NaN;
fprintf('\nBuilding AMG hierarchy + ichol ...\n');
t0 = tic;
try
    amg = amg_setup(A, opts);
    t_setup = toc(t0);
    fprintf('Setup time: %.1fs\n', t_setup);
catch ME
    fprintf('AMG setup FAILED after %.1fs: %s\n', toc(t0), ME.message);
    fprintf('Falling back to plain ichol(0)-PCG.\n');
end

t0 = tic;
if ~isempty(amg)
    fprintf('\nPCG with AMG+ichol (tol=%.0e, maxit=%d) ...\n', tol, maxiter);
    precond = @(r) amg_precondition(amg, r);
    [x, flag, relres, iter, resvec] = pcg(A, b, tol, maxiter, precond);
else
    fprintf('\nPCG with ichol(0) only (tol=%.0e, maxit=%d) ...\n', tol, maxiter);
    L = ichol(A, struct('type','nofill'));
    [x, flag, relres, iter, resvec] = pcg(A, b, tol, maxiter, L, L');
end
t_solve = toc(t0);
fprintf('flag=%d  iter=%d  relres=%.3e  time=%.1fs\n', flag, iter, relres, t_solve);

save('problem_case_solution.mat','x','iter','relres','flag','resvec', ...
     't_setup','t_solve','-v7.3');
fprintf('Saved problem_case_solution.mat\n');
