% run_problem_case_instrumented.m
% Same as run_problem_case.m but prints per-stage timing + RSS memory,
% and instruments amg_setup stage-by-stage so we can see exactly where
% it crashes / where memory blows up.

function run_problem_case_instrumented()
cd('/Users/david/Desktop/agg_pcg_ba1x/blood_networks_new/C_amg_release/src');

pid  = feature('getpid');
mem0 = rss_mb(pid);
T0   = tic;
stage = @(label) fprintf('[%7.1fs | RSS %8.1f MB | +%8.1f MB] %s\n', ...
                         toc(T0), rss_mb(pid), rss_mb(pid)-mem0, label);

stage('start');

stage('loading problem_case.mat ...');
S = load('/Users/david/Desktop/agg_pcg_ba1x/blood_networks_new/problem_case.mat');
A = S.A; b = S.b; clear S;
stage(sprintf('loaded: n=%d nnz=%d ||b||=%.3e', size(A,1), nnz(A), norm(b)));

if ~issparse(A), A = sparse(A); stage('converted to sparse'); end

opts.threshold  = 0.25;
opts.min_size   = 512;
opts.max_levels = 25;
opts.nu         = 2;
tol     = 1e-6;
maxiter = 500;

% ----- inlined + instrumented amg_setup -----
threshold  = opts.threshold;
min_size   = opts.min_size;
max_levels = opts.max_levels;
nu         = opts.nu;

levels = {};
levels{1}.A = A;
nlev = 1;
Ac = A;

stage(sprintf('AMG begin: level 1 n=%d nnz=%d', size(Ac,1), nnz(Ac)));

while size(Ac,1) >= min_size && nlev < max_levels
    n_fine = size(Ac,1);

    stage(sprintf('  L%d rs_strength (n=%d nnz=%d) ...', nlev, n_fine, nnz(Ac)));
    [Sstr, max_lam] = rs_strength(Ac, threshold);
    stage(sprintf('  L%d rs_strength done: nnz(S)=%d, lam_max=%d', nlev, nnz(Sstr), max_lam));

    stage(sprintf('  L%d rs_coarsening ...', nlev));
    [cf, n_coarse] = rs_coarsening(Sstr, max_lam);
    stage(sprintf('  L%d rs_coarsening done: n_coarse=%d', nlev, n_coarse));

    if n_coarse == 0 || n_coarse >= n_fine
        stage(sprintf('  coarsening stalled at level %d', nlev));
        break;
    end

    cf_order = [find(cf==1); find(cf~=1)];
    levels{nlev}.cf_order = cf_order;

    stage(sprintf('  L%d rs_interpolation ...', nlev));
    P = rs_interpolation(Ac, Sstr, cf);
    stage(sprintf('  L%d rs_interpolation done: nnz(P)=%d', nlev, nnz(P)));
    levels{nlev}.P = P;

    clear Sstr;

    stage(sprintf('  L%d Galerkin Ac = P''*Ac*P ...', nlev));
    Ac_new = P' * Ac * P;
    stage(sprintf('  L%d Galerkin done: n=%d nnz=%d', nlev, size(Ac_new,1), nnz(Ac_new)));

    nlev = nlev + 1;
    levels{nlev}.A = Ac_new;
    Ac = Ac_new;
    clear Ac_new;
end

stage(sprintf('AMG hierarchy built: %d levels', nlev));

stage('ichol(A) ...');
try
    Lf = ichol(A, struct('type','nofill'));
    stage(sprintf('ichol done: nnz(L)=%d', nnz(Lf)));
catch ME
    stage(['ichol FAILED: ' ME.message]);
    Lf = [];
end

amg.levels = levels;
amg.num_levels = nlev;
amg.nu = nu;
amg.A = A;
amg.L = Lf;

stage('PCG begin ...');
precond = @(r) amg_precondition(amg, r);
[x, flag, relres, iter, resvec] = pcg(A, b, tol, maxiter, precond);
stage(sprintf('PCG done: flag=%d iter=%d relres=%.3e', flag, iter, relres));

save('problem_case_solution.mat','x','iter','relres','flag','resvec','-v7.3');
stage('saved solution');
end

function mb = rss_mb(p)
[~, s] = system(sprintf('ps -o rss= -p %d', p));
mb = str2double(strtrim(s)) / 1024;
end
