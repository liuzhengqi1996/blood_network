function amg = amg_setup(A, opts)
% Build classical RS AMG hierarchy (PMIS coarsening, standard interpolation).
%   opts.threshold, opts.min_size, opts.max_levels, opts.nu
addpath(fileparts(mfilename('fullpath')));

threshold  = getfield_d(opts,'threshold',0.25);
min_size   = getfield_d(opts,'min_size',512);
max_levels = getfield_d(opts,'max_levels',20);

amg.nu = getfield_d(opts,'nu',2);
amg.threshold = threshold;

fprintf('Building PMIS hierarchy (theta=%.2f, min_size=%d)\n', threshold, min_size);

lev = 1;
amg.A{lev}    = A;
amg.l1{lev}   = full(sum(abs(A),2));
fprintf('  Level 1: n=%d nnz=%d\n', size(A,1), nnz(A));

while size(amg.A{lev},1) >= min_size && lev < max_levels
    Ac = amg.A{lev};
    n_fine = size(Ac,1);

    S  = rs_strength(Ac, threshold);
    cf = pmis_coarsening(S);
    nC = nnz(cf == 1);
    if nC == 0 || nC >= n_fine
        fprintf('  Coarsening stalled at level %d (nC=%d)\n', lev, nC);
        break;
    end

    P = rs_standard_interpolation(Ac, S, cf);
    amg.P{lev} = P;

    Anew = P.' * Ac * P;
    Anew = (Anew + Anew.') * 0.5;

    lev = lev + 1;
    amg.A{lev}  = Anew;
    amg.l1{lev} = full(sum(abs(Anew),2));
    fprintf('  Level %d: n=%d nnz=%d (ratio %.2f)\n', ...
        lev, size(Anew,1), nnz(Anew), n_fine/nC);
end
amg.nlev = lev;

% Coarsest-level ichol for the coarse solve
Ac = amg.A{lev};
try
    amg.Lc = ichol(Ac, struct('type','nofill','diagcomp',1e-6));
catch
    amg.Lc = ichol(Ac, struct('type','nofill','diagcomp',1e-3));
end
fprintf('Hierarchy: %d levels, coarsest n=%d\n\n', lev, size(Ac,1));
end

function v = getfield_d(s, f, d)
if isfield(s,f), v = s.(f); else, v = d; end
end
