% test_amg.m - Test AMG+ichol preconditioned PCG solver (pure MATLAB)
%
% Authors: Z. Liu and Ludmil Zikatanov
%
% Reads raw network data directly from four_large_Cartesian_networks/,
% builds the system matrix, and solves — no intermediate .dat file needed.
%
% Usage:
%   >> test_amg                                           % default (small + medium)
%   >> test_amg({'100x100x10'})                           % one specific size
%   >> test_amg({'100x100x10', '100x100x100', '100x100x1000'})

function test_amg(sizes)

dname = 'four_large_Cartesian_networks';

if nargin < 1
    sizes = {'100x100x10', '100x100x100'};
elseif ischar(sizes) || isstring(sizes)
    sizes = {char(sizes)};
end

% AMG parameters
opts.threshold  = 0.25;
opts.min_size   = 128;
opts.max_levels = 20;
opts.nu         = 2;
tol        = 1e-6;
maxiter    = 200;

for f = 1:length(sizes)
    sz = sizes{f};
    filename = fullfile(dname, ['nwk_Cartesian_' sz]);

    fprintf('\n############################################################\n');
    fprintf('  %s\n', sz);
    fprintf('############################################################\n');

    %% Build system matrix from raw network data
    fprintf('Reading network data from %s ...\n', filename);
    t0 = tic;
    [A, b, interior_nodes] = read_nwk(filename);
    A
    n = size(A, 1);
    fprintf('  n = %d, nnz = %d  (%.2f s)\n', n, nnz(A), toc(t0));

    bnorm = norm(b);

    %% AMG+ichol setup
    fprintf('Building AMG hierarchy + ichol ...\n');
    t0 = tic;
    amg = amg_setup(A, opts);
    t_setup = toc(t0);
    fprintf('  Setup time: %.2f s\n', t_setup);

    %% PCG with AMG+ichol preconditioner
    fprintf('PCG with AMG+ichol:\n');
    precond = @(r) amg_precondition(amg, r);
    t0 = tic;
    [x_amg, ~, ~, iter] = pcg(A, b, tol, maxiter, precond);
    t_solve = toc(t0);
    Ax_amg = A*x_amg;
    relres_amg = norm(b - Ax_amg) / bnorm;
    fprintf('  iter = %d, relres = %.2e, time = %.2f s\n', iter, relres_amg, t_solve);
    fprintf('  ||A*u|| = %.6e,  ||b|| = %.6e\n', norm(Ax_amg), bnorm);

    %% Save solution p_I to file (one value per line, ordered by interior node)
    solfile = sprintf('p_interior_%s.txt', sz);
    fid = fopen(solfile, 'w');
    fprintf(fid, '%.16e\n', x_amg);
    fclose(fid);
    fprintf('  Solution saved to %s\n', solfile);

    %% Backslash (skip for large problems to avoid out-of-memory)
    t_bs = NaN;
    relres_bs = NaN;
    if n <= 2000000
        fprintf('MATLAB backslash:\n');
        t0 = tic;
        x_bs = A \ b;
        t_bs = toc(t0);
        Ax_bs = A*x_bs;
        relres_bs = norm(b - Ax_bs) / bnorm;
        fprintf('  relres = %.2e, time = %.2f s\n', relres_bs, t_bs);
        fprintf('  ||A*u|| = %.6e,  ||b|| = %.6e\n', norm(Ax_bs), bnorm);
    else
        fprintf('MATLAB backslash: skipped (n = %d too large)\n', n);
    end

    %% Summary
    fprintf('\n  %-25s %6s %12s %12s %12s\n', ...
        'Solver', 'Iter', 'Setup(s)', 'Solve(s)', 'Relres');
    fprintf('  %-25s %6d %12.2f %12.2f %12.2e\n', ...
        'AMG+ichol', iter, t_setup, t_solve, relres_amg);
    if ~isnan(t_bs)
        fprintf('  %-25s %6s %12s %12.2f %12.2e\n', ...
            'backslash', '--', '--', t_bs, relres_bs);
    end

end

end
