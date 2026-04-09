function amg = amg_setup(A, opts)
% AMG_SETUP  Build classical RS AMG hierarchy with fine-level ichol.
%
%   amg = amg_setup(A, opts)
%
%   opts fields (all optional):
%     threshold  - strength threshold (default 0.25)
%     min_size   - stop coarsening below this size (default 128)
%     max_levels - maximum number of levels (default 20)
%     nu         - number of GS sweeps per smoothing step (default 2)
%
%   Returns struct with:
%     levels    - cell array, each with fields A, P, cf_order
%     L         - fine-level incomplete Cholesky factor
%     nu        - smoothing sweeps
%     A         - fine-level matrix (same as levels{1}.A)
%     num_levels

if nargin < 2, opts = struct(); end
threshold  = getfield_default(opts, 'threshold',  0.25);
min_size   = getfield_default(opts, 'min_size',   128);
max_levels = getfield_default(opts, 'max_levels', 20);
nu         = getfield_default(opts, 'nu',         2);

fprintf('Building RS hierarchy (threshold=%.2f, min_size=%d)\n', ...
    threshold, min_size);
fprintf('  Level 1: n = %d, nnz = %d\n', size(A,1), nnz(A));

levels = {};
levels{1}.A = A;

nlev = 1;
Ac = A;

while size(Ac, 1) >= min_size && nlev < max_levels
    n_fine = size(Ac, 1);

    % Step 1: Strength of connection
    [S, max_lam] = rs_strength(Ac, threshold);

    % Step 2: Classical RS coarsening
    [cf, n_coarse] = rs_coarsening(S,max_lam);

    if n_coarse == 0 || n_coarse >= n_fine
        fprintf('  Coarsening stalled at level %d (n_coarse=%d)\n', ...
            nlev, n_coarse);
        break;
    end

    % Build CF ordering: C-points first, then F-points
    cf_order = [find(cf == 1); find(cf ~= 1)];
    levels{nlev}.cf_order = cf_order;

    % Step 3: Standard interpolation
    P = rs_interpolation(Ac, S, cf);
    levels{nlev}.P = P;

    % Step 4: Galerkin coarse-grid operator Ac = P' * A * P
    Ac_new = P' * Ac * P;

    % P'*A*P is symmetric in exact arithmetic, but sparse mat-mat in
    % floating point can leave ULP-level asymmetry. That perturbs the
    % strength-of-connection threshold on the next level and forks the
    % coarsening trajectory away from the c_old reference. Force symmetry.
    Ac_new = (Ac_new + Ac_new') / 2;

    fprintf('  Level %d: n = %d, nnz = %d (ratio %.2f)\n', ...
        nlev + 1, n_coarse, nnz(Ac_new), n_fine / n_coarse);

    nlev = nlev + 1;
    levels{nlev}.A = Ac_new;
    Ac = Ac_new;
end

amg.levels = levels;
amg.num_levels = nlev;
amg.nu = nu;
amg.A = levels{1}.A;

% Compute ichol on fine level
amg.L = ichol(A, struct('type', 'nofill'));

fprintf('RS hierarchy complete: %d levels\n\n', nlev);
end

function val = getfield_default(s, name, default)
if isfield(s, name)
    val = s.(name);
else
    val = default;
end
end
