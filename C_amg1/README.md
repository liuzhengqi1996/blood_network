# Classical Ruge-Stuben AMG

C implementation of the classical algebraic multigrid (AMG) method of
Ruge and Stuben [1], including strength of connection, C/F splitting,
and standard (direct) interpolation.  Integrates with the existing
aggregation-based AMG framework via the `dcsr_mat` CSR data structure
from `agg_lib.h`.

## Project layout

```
jl_and_C_amg/
├── include/            Headers (agg_lib.h, amg.h)
├── src/                Library sources (10 .c files)
├── lib/                Shared library output (libamg.so)
├── test_amg.c          C test driver
├── solve_amg.jl        Julia interface script
├── Makefile
└── matrixA.csr   Sample matrix (COO format)
└── matrixP.csr   Sample matrix (COO format)
```

## Build

```
make            # builds test_amg and lib/libamg.so
make clean      # removes all objects, test_amg, and lib/libamg.so
```

Both `test_amg` and `lib/libamg.so` are built by default.  The executable
links dynamically against the shared library via `-Llib -lamg` with an
`$ORIGIN/lib` rpath, so it finds `libamg.so` at runtime without setting
`LD_LIBRARY_PATH`.

## Running the C test

```
./test_amg <matrix_file> [threshold] [nu] [min_size]
```

| Argument | Default | Description |
|---|---|---|
| `matrix_file` | (required) | Path to a COO-format matrix file |
| `threshold` | 0.25 | Strength-of-connection threshold |
| `nu` | 2 | Number of smoothing sweeps per V-cycle |
| `min_size` | 4096 | Coarsening stops below this size, for small examples use 32 |

If `rhs.txt` exists (format: first line is `n`, followed by `n` values),
it is loaded as the right-hand side; otherwise a random RHS is used.

### Examples

**matrixA.csr** (n=137, nnz=1665):

```
$ ./test_amg matrixA.csr 0.125 2 25
  Level 1: n = 137, nnz = 1665
  Level 2: n = 41, nnz = 449 (ratio 3.34)
  Level 3: n = 18, nnz = 146 (ratio 2.28)
RS hierarchy complete: 3 levels

  AMG+ichol:    4 iterations, relres = 6.95e-12
  AMG V-cycle:  7 iterations, relres = 2.17e-11
```

**matrixP.csr** (n=137, nnz=1081):

```
$ ./test_amg matrixP.csr 0.125 2 25
  Level 1: n = 137, nnz = 1081
  Level 2: n = 41, nnz = 377 (ratio 3.34)
  Level 3: n = 14, nnz = 100 (ratio 2.93)
RS hierarchy complete: 3 levels

  AMG+ichol:    3 iterations, relres = 4.95e-11
  AMG V-cycle:  5 iterations, relres = 2.27e-11
```

## Julia interface

### Requirements

Julia (>= 1.6) with the standard library modules `SparseArrays`,
`LinearAlgebra`, `Random`, and `Printf` (all ship with Julia, no
packages to install).

### Usage

```
julia solve_amg.jl <matrix_file> [threshold] [nu] [min_size]
```

| Argument | Default | Description |
|---|---|---|
| `matrix_file` | (required) | Path to a COO-format matrix file |
| `threshold` | 0.25 | Strength-of-connection threshold |
| `nu` | 2 | Number of smoothing sweeps per V-cycle |
| `min_size` | 4096 | Coarsening stops below this size |

Same `rhs.txt` convention as the C driver.

### Examples

**matrixA.csr**:

```
$ julia solve_amg.jl matrixA.csr 0.125 2 25
  AMG levels: 3
  PCG with AMG+ichol: 4 iterations, relres = 1.91e-12
```

**matrixP.csr**:

```
$ julia solve_amg.jl matrixP.csr 0.125 2 25
  AMG levels: 3
  PCG with AMG+ichol: 3 iterations, relres = 3.73e-11
```

### How the Julia-C binding works

The script uses Julia's built-in `ccall` to call functions exported by
`lib/libamg.so` directly, with no wrapper packages needed.  The shared
library exposes an opaque-handle API through five C functions:

| C function | Signature | Purpose |
|---|---|---|
| `amg_solver_create` | `(n, nnz, ia, ja, val, threshold, min_size, max_levels, nu) → void*` | Build AMG hierarchy + ichol factor from 0-based CSR arrays |
| `amg_solver_solve` | `(handle, b, x, tol, maxiter) → int` | PCG solve with AMG+ichol preconditioner; returns iteration count |
| `amg_solver_size` | `(handle) → int` | Matrix dimension |
| `amg_solver_num_levels` | `(handle) → int` | Number of AMG levels |
| `amg_solver_free` | `(handle) → void` | Release all resources |

The Julia side performs the following steps:

1. **Read COO file** into a Julia `SparseMatrixCSC` via `sparse()`.
2. **Symmetrize** with `A = (A + A') / 2`.
3. **Convert to 0-based CSR** by transposing the CSC structure
   (`SparseMatrixCSC(A')` gives CSC of $A^T$ = CSR of $A$) and
   subtracting 1 from `colptr` and `rowval`.
4. **`ccall(:amg_solver_create, ...)`** passes the CSR arrays as
   `Ptr{Cint}` / `Ptr{Cdouble}` — Julia pins the arrays for the
   duration of the call and the C side copies them internally.
5. **`ccall(:amg_solver_solve, ...)`** fills the solution vector
   in-place.
6. **Verify** the residual in Julia: `norm(b - A*x) / norm(b)`.
7. **`ccall(:amg_solver_free, ...)`** releases the C-side memory.

The library path is resolved at compile time as a `const`:

```julia
const LIBAMG = joinpath(@__DIR__, "lib", "libamg.so")
```

This is required because `ccall` needs the library to be a constant
expression.  Since the script lives alongside the Makefile, `@__DIR__`
resolves to the project root and `lib/libamg.so` is found automatically.

## Source files

| File | Location | Description |
|---|---|---|
| `rs_classical.c` | `src/` | Strength of connection, RS coarsening, standard interpolation |
| `build_hierarchy.c` | `src/` | AMG hierarchy construction |
| `amg_cycle.c` | `src/` | V-cycle (forward/backward Gauss-Seidel) |
| `ichol.c` | `src/` | Incomplete Cholesky factorization |
| `amg_ichol_precond.c` | `src/` | Combined AMG + ichol preconditioner |
| `amg_pcg.c` | `src/` | Preconditioned conjugate gradient solver |
| `amg_wrapper.c` | `src/` | Shared library wrapper (opaque-handle API) |
| `vector_ops.c` | `src/` | Vector operations (copy, dot, norm, axpy) |
| `sparse_ops.c` | `src/` | CSR matrix operations (SpMV, RAP) |
| `io.c` | `src/` | COO file reader, COO-to-CSR conversion |
| `agg_lib.h` | `include/` | Core types (`dcsr_mat`, `dcoo_mat`) and inline utilities |
| `amg.h` | `include/` | AMG types and function prototypes |

## Algorithms

### 1. Strength of Connection

**Function:** `rs_strength(A, theta, S)`

Given a sparse matrix $A$ and a threshold $\theta \in (0,1)$ (typically
$\theta = 0.25$), point $j$ is a **strong connection** of point $i$ if

$$
-a_{ij} \;\ge\; \theta \cdot \max_{k \ne i}\bigl(-a_{ik}\bigr).
$$

This identifies couplings whose magnitude is within a factor $\theta$ of
the largest off-diagonal coupling in the row.  The output $S$ is a sparse
matrix containing only the strong off-diagonal entries of $A$.

**Rationale.**  For an M-matrix (positive diagonal, non-positive
off-diagonal), the off-diagonal entries $a_{ij} \le 0$ represent
couplings between unknowns.  Algebraically smooth error (error not
reduced by relaxation) varies slowly along strong connections and can
change rapidly only across weak ones [1, Section 3].  The strength
threshold separates connections that matter for interpolation from
those that can be safely ignored.

**Complexity:** $O(\text{nnz}(A))$ — three passes over the matrix
(compute row maxima, count strong entries, fill).

### 2. C/F Splitting (Classical RS Coarsening)

**Function:** `rs_coarsening(S, cf, n_coarse)`

Partitions unknowns into **coarse points** (C) that form the next-level
grid and **fine points** (F) that are interpolated from C-points.  The
splitting must satisfy two heuristic principles [1]:

- **(H1)** For each F-point $i$, every strong connection of $i$ is either
  a C-point or is itself strongly connected to a C-point in $C_i$.
  (Ensures interpolation can reach all important couplings.)

- **(H2)** The set of C-points forms a maximal independent set with
  respect to strong connections.
  (Keeps the coarse grid small.)

**Algorithm (first pass).**  Uses the *influence measure*:

$$
\lambda_i = \bigl|\{j : i \in S_j\}\bigr|
$$

which counts how many points are strongly influenced by $i$.  This
requires the transpose $S^T$.  The greedy procedure is:

1. Select the undecided point $i$ with largest $\lambda_i$.
2. Mark $i$ as C.
3. For each undecided point $j$ that $i$ strongly influences ($j \in S_i^T$):
   - Mark $j$ as F.
   - For each undecided point $k$ that $j$ strongly depends on ($k \in S_j$, $k \ne i$):
     increment $\lambda_k$ (making $k$ more attractive as a future C-point,
     since the newly created F-point $j$ needs C-neighbors for interpolation).
4. Repeat until no undecided points remain.

When $j$ becomes an F-point, incrementing $\lambda$ for its
strong dependencies ensures that nearby undecided points are more
likely to be selected as C-points in subsequent iterations.  This is
the mechanism that enforces heuristic (H1).

**Second pass.**  Any F-point that ended up without a strong C-neighbor
is promoted to C.  This can happen in rare configurations where the
greedy first pass leaves an F-point isolated from all C-points.

**Complexity:** $O(n^2)$ in the current implementation due to the
linear scan for the maximum $\lambda$.  Production codes use a bucket
or heap structure to achieve $O(n + \text{nnz}(S))$; see [3, Appendix A].

### 3. Standard (Direct) Interpolation

**Function:** `rs_standard_interpolation(A, S, cf, P)`

Constructs the prolongation operator $P$ ($n \times n_c$) that maps
coarse-grid corrections to the fine grid.

**C-points.**  Row $i$ of $P$ has a single entry $P_{i,\hat\imath} = 1$
where $\hat\imath$ is the coarse-grid index of $i$ (identity injection).

**F-points.**  For each F-point $i$, classify its neighbors:

| Set | Definition | Role |
|---|---|---|
| $C_i$ | $\{j \in S_i : j \text{ is C}\}$ | Strong C-neighbors (interpolation sources) |
| $D_i^s$ | $\{k \in S_i : k \text{ is F}\}$ | Strong F-neighbors (redistributed) |
| $D_i^w$ | $\{m \in N_i \setminus S_i : m \ne i\}$ | Weak connections (lumped to diagonal) |

**Derivation.**  Assume the error $e$ is algebraically smooth,
satisfying $Ae \approx 0$ pointwise.  Row $i$ gives

$$
a_{ii}\,e_i + \sum_{j \in N_i} a_{ij}\,e_j \;\approx\; 0.
$$

Partition $N_i$ into $C_i$, $D_i^s$, and $D_i^w$:

**Step 1 — Lump weak connections.**
For $m \in D_i^w$, approximate $e_m \approx e_i$ (weak connections
imply the error varies little):

$$
\hat{a}_{ii}\,e_i + \sum_{j \in C_i} a_{ij}\,e_j
  + \sum_{k \in D_i^s} a_{ik}\,e_k \;\approx\; 0
$$

where $\hat{a}_{ii} = a_{ii} + \sum_{m \in D_i^w} a_{im}$.

**Step 2 — Redistribute strong F-neighbors.**
For $k \in D_i^s$, approximate $e_k$ by interpolation from $C_i$:

$$
e_k \;\approx\; \sum_{j \in C_i}
  \frac{a_{kj}}{\displaystyle\sum_{\ell \in C_i} a_{k\ell}}\;e_j.
$$

This distributes the influence of F-neighbor $k$ to the C-points
in $C_i$, weighted by how strongly $k$ is connected to each.  If
$k$ has no connection to any point in $C_i$ (the denominator vanishes),
then $a_{ik}$ is added to $\hat{a}_{ii}$ instead.

**Step 3 — Final weights.**
Substituting into the smoothness equation and solving for $e_i$:

$$
e_i \;\approx\; \sum_{j \in C_i} w_{ij}\,e_j
$$

where

$$
\boxed{
w_{ij} = -\frac{1}{\hat{a}_{ii}}
  \left(a_{ij} + \sum_{k \in D_i^s}
    \frac{a_{ik}\,a_{kj}}{\displaystyle\sum_{\ell \in C_i} a_{k\ell}}
  \right).
}
$$

**Properties.**

- **Constant preservation:** If row $i$ of $A$ has zero row sum
  ($\sum_j a_{ij} = 0$, typical for interior nodes of a Laplacian),
  then $\sum_{j \in C_i} w_{ij} = 1$, so the constant vector is
  interpolated exactly.

- **Non-negative weights:** For M-matrices ($a_{ii} > 0$,
  $a_{ij} \le 0$ for $j \ne i$), all interpolation weights satisfy
  $w_{ij} \ge 0$.

**Complexity:** $O(\text{nnz}(A) \cdot \max\_\text{row\_nnz})$ in
the worst case due to the inner loop over row $k$ for each strong
F-neighbor.  In practice this is close to $O(\text{nnz}(A))$ for
matrices with bounded stencil width.

**Implementation notes.**  The code uses sentinel-based marker arrays
(`is_strong[j] == i` marks $j \in S_i$) to avoid per-row clearing,
and a position array (`pos[j]` gives the offset in $P$'s value array)
for $O(1)$ accumulation into interpolation weights.

## References

[1] J.W. Ruge and K. Stuben, "Algebraic Multigrid," in *Multigrid
    Methods* (S.F. McCormick, ed.), Frontiers in Applied Mathematics,
    Vol. 3, SIAM, Philadelphia, 1987, pp. 73-130.

[2] K. Stuben, "A Review of Algebraic Multigrid," *Journal of
    Computational and Applied Mathematics*, 128(1-2):281-309, 2001.
    (Section 3: strength of connection; Section 3.3: standard
    interpolation formula; Appendix A: coarsening algorithms.)

[3] U. Trottenberg, C.W. Oosterlee, and A. Schuller, *Multigrid*,
    Academic Press, 2001.  (Chapter 8: Algebraic Multigrid.)

[4] Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed.,
    SIAM, 2003.  (Chapter 13: Multigrid Methods.)

[5] Van Emden Henson and Ulrike Meier Yang, "BoomerAMG: A Parallel
    Algebraic Multigrid Solver and Preconditioner," *Applied Numerical
    Mathematics*, 41(1):155-177, 2002.
    (Reference implementation in hypre; discusses parallel RS coarsening
    variants including CLJP and PMIS.)
