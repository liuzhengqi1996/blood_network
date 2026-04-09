function P = rs_standard_interpolation(A, S, cf)
% Classical RS standard (direct) interpolation.
% A,S are n x n sparse symmetric, cf in {+1,-1}.  Returns P (n x nc).
n  = size(A,1);
isC = (cf == 1);
nc = nnz(isC);
cmap = zeros(n,1);
cmap(isC) = (1:nc).';

% Row-major access via transposes (A.' is CSC -> entries grouped by row of A)
[Ap, Aj, Av] = csr(A);
[Sp, Sj, ~ ] = csr(S);

% Reserve enough triplets: 1 per C-row + (#strong C-nbrs) per F-row
% Upper bound: nnz(S) + nc
maxnnz = nnz(S) + nc;
Ir = zeros(maxnnz,1);
Jc = zeros(maxnnz,1);
Vv = zeros(maxnnz,1);
ptr = 0;

is_strong = false(n,1);
pos = zeros(n,1);

for i = 1:n
    if isC(i)
        ptr = ptr + 1;
        Ir(ptr) = i; Jc(ptr) = cmap(i); Vv(ptr) = 1.0;
        continue;
    end

    s_nbrs = Sj(Sp(i)+1 : Sp(i+1));
    if isempty(s_nbrs), continue; end
    is_strong(s_nbrs) = true;

    Ci_mask = isC(s_nbrs);
    Ci = s_nbrs(Ci_mask);
    if isempty(Ci)
        is_strong(s_nbrs) = false;
        continue;
    end
    nCi = numel(Ci);
    pos(Ci) = 1:nCi;
    p_vals = zeros(nCi,1);

    a_cols = Aj(Ap(i)+1 : Ap(i+1));
    a_vals = Av(Ap(i)+1 : Ap(i+1));

    diag_v = 0; weak = 0;
    for kk = 1:numel(a_cols)
        j = a_cols(kk); a = a_vals(kk);
        if j == i
            diag_v = a;
        elseif is_strong(j) && isC(j)
            p_vals(pos(j)) = p_vals(pos(j)) + a;
        elseif is_strong(j) && ~isC(j)
            % strong F-neighbor, handled in next loop
        else
            weak = weak + a;
        end
    end
    diag_hat = diag_v + weak;

    for kk = 1:numel(a_cols)
        kf = a_cols(kk);
        if kf == i, continue; end
        if ~is_strong(kf) || isC(kf), continue; end
        a_ik = a_vals(kk);

        k_cols = Aj(Ap(kf)+1 : Ap(kf+1));
        k_vals = Av(Ap(kf)+1 : Ap(kf+1));
        sel = isC(k_cols) & is_strong(k_cols);
        denom = sum(k_vals(sel));
        if abs(denom) < 1e-15
            diag_hat = diag_hat + a_ik;
        else
            sc = k_cols(sel); sv = k_vals(sel);
            for tt = 1:numel(sc)
                p_vals(pos(sc(tt))) = p_vals(pos(sc(tt))) + a_ik*sv(tt)/denom;
            end
        end
    end

    if abs(diag_hat) > 1e-15
        p_vals = p_vals / (-diag_hat);
    end

    rng_ = ptr+1 : ptr+nCi;
    Ir(rng_) = i;
    Jc(rng_) = cmap(Ci);
    Vv(rng_) = p_vals;
    ptr = ptr + nCi;

    is_strong(s_nbrs) = false;
    pos(Ci) = 0;
end

P = sparse(Ir(1:ptr), Jc(1:ptr), Vv(1:ptr), n, nc);
end

function [rp, ci, vv] = csr(A)
% Build CSR row pointers/indices/values from a sparse MATLAB matrix.
% rp is 0-based, length n+1.
n = size(A,1);
At = A.';            % CSC of A^T == row-major view of A
[ci, jj, vv] = find(At);
% In find(At): row=ci is column index of A, col=jj is row index of A.
% Entries are returned in column-major order of At, i.e. grouped by row of A.
rp = [0; cumsum(accumarray(jj, 1, [n 1]))];
end
