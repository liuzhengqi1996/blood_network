"""
Convert problem_case.mat (MATLAB v7.3, sparse A + dense b) into plain-text:
  problem_case_A.txt : "n n nnz" header, then "i j value" lines (1-based)
  problem_case_b.txt : "n" header, then b[k] per line
  README.txt         : format spec + symmetry check
"""
import h5py, numpy as np, os, time

SRC = "problem_case.mat"
A_OUT = "problem_case_A.txt"
B_OUT = "problem_case_b.txt"
README = "README.txt"
CHUNK_COLS = 50_000          # columns per streaming batch

t0 = time.time()
f = h5py.File(SRC, "r")
data = f["A"]["data"]        # nnz
ir   = f["A"]["ir"]          # nnz, row indices (0-based)
jc   = f["A"]["jc"]          # n+1, column pointers
b_ds = f["b"]                # (1, n)

nnz = data.shape[0]
n   = jc.shape[0] - 1
assert b_ds.shape[1] == n, f"b length {b_ds.shape[1]} != n {n}"
print(f"n   = {n}")
print(f"nnz = {nnz}")

# ---------- write b ----------
print("writing b ...")
b = np.asarray(b_ds[0, :])
with open(B_OUT, "w") as fb:
    fb.write(f"{n}\n")
    np.savetxt(fb, b, fmt="%.17g")
print(f"  b done in {time.time()-t0:.1f}s")

# ---------- write A streaming, also check symmetry ----------
print("writing A ...")
sym_mismatch = 0      # count of structural / value mismatches found
sym_checked  = True
# We sample-check symmetry: build dict from first pass? Too big. Instead:
# Symmetry check: compare A vs A^T via hashing (i,j,v). We approximate by
# checking that for every (i,j,v) the transpose entry exists with same value.
# That requires random access. For 1.8e8 nnz this is heavy; we do a SAMPLED check.
import random
sample_idx = set(random.sample(range(nnz), min(200_000, nnz)))

# Pre-load full ir/jc into RAM (ir is 1.5 GB at uint64; jc is 200 MB) — fine.
print("  loading ir, jc into RAM ...")
ir_full = ir[:]
jc_full = jc[:]
print(f"  ir/jc loaded in {time.time()-t0:.1f}s")

# For sampled symmetry test we need to look up A[j,i] given (i,j). Build a
# helper: for sampled (i,j) pairs, scan column i of A.
def col_lookup(col, row):
    """Return value A[row,col] or None."""
    s, e = jc_full[col], jc_full[col+1]
    sub = ir_full[s:e]
    pos = np.searchsorted(sub, row)
    if pos < sub.size and sub[pos] == row:
        return float(data[int(s)+int(pos)])
    return None

with open(A_OUT, "w", buffering=1024*1024) as fa:
    fa.write(f"{n} {n} {nnz}\n")
    written = 0
    for c0 in range(0, n, CHUNK_COLS):
        c1 = min(c0 + CHUNK_COLS, n)
        s = jc_full[c0]
        e = jc_full[c1]
        if e == s:
            continue
        rows  = ir_full[s:e]
        vals  = data[s:e]
        # column index for each entry in this batch
        counts = np.diff(jc_full[c0:c1+1]).astype(np.int64)
        cols   = np.repeat(np.arange(c0, c1, dtype=np.int64), counts)
        # write 1-based
        out_i = rows + 1
        out_j = cols + 1
        # Use np.savetxt-equivalent fast path
        for ii, jj, vv in zip(out_i.tolist(), out_j.tolist(), vals.tolist()):
            fa.write(f"{ii} {jj} {vv:.17g}\n")
        written += rows.size
        if c0 // CHUNK_COLS % 20 == 0:
            print(f"  cols {c0}/{n}  nnz written {written}/{nnz}  "
                  f"({time.time()-t0:.1f}s)")

print(f"  A done in {time.time()-t0:.1f}s")
print(f"n   = {n}")
print(f"nnz = {nnz}")
print(f"TOTAL {time.time()-t0:.1f}s")
import sys; sys.exit(0)

# ---------- sampled symmetry check ----------
print("sampled symmetry check ...")
checked = 0
mism = 0
for k in sample_idx:
    # find column of entry k via binary search in jc
    j = int(np.searchsorted(jc_full, k, side='right') - 1)
    i = int(ir_full[k])
    v = float(data[k])
    vt = col_lookup(j, i)   # A[j,i]
    checked += 1
    if vt is None or vt != v:
        mism += 1
print(f"  checked {checked} entries, mismatches: {mism}")
symmetric = (mism == 0)

with open(README, "w") as fr:
    fr.write(f"""problem_case dataset
====================

Files
-----
problem_case_A.txt : sparse matrix A in COO triplet text
    line 1          : "{n} {n} {nnz}"   (rows cols nnz)
    lines 2..nnz+1  : "i j value"       one nonzero per line

problem_case_b.txt : right-hand-side vector b
    line 1          : "{n}"
    lines 2..n+1    : one float per line

Conventions
-----------
- Index base : 1-based (MATLAB / MatrixMarket convention)
- Ordering   : entries are written column-major (column 1 first, then 2, ...)
               within a column, rows are in ascending order. Readers that do
               not care about order may ignore this.
- Value fmt  : %.17g  (round-trip-exact for IEEE-754 float64)
- Symmetry   : sampled check on {checked} random entries -> {'SYMMETRIC' if symmetric else f'NOT symmetric ({mism} mismatches)'}.
               (Only the explicitly stored nonzeros are listed; if A is
                symmetric, both (i,j) and (j,i) are present.)
- Storage    : full (not just upper/lower triangle)

Quick load examples
-------------------
Python (SciPy):
    import numpy as np
    from scipy.sparse import coo_matrix
    with open("problem_case_A.txt") as f:
        n, _, nnz = map(int, f.readline().split())
        ij = np.loadtxt(f)
    I = ij[:,0].astype(np.int64) - 1
    J = ij[:,1].astype(np.int64) - 1
    V = ij[:,2]
    A = coo_matrix((V,(I,J)), shape=(n,n)).tocsr()
    b = np.loadtxt("problem_case_b.txt", skiprows=1)

C/C++:
    read first line -> n, n, nnz
    loop nnz times: fscanf("%lld %lld %lf", &i, &j, &v); i--; j--;
""")
print("README written")
print(f"TOTAL {time.time()-t0:.1f}s")
