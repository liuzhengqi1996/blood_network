/*
 * amg_wrapper.c - Shared library wrapper for AMG + ichol solver
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Exposes an opaque-handle API for FFI (Julia, Python, etc.):
 *   amg_solver_create      - build hierarchy + ichol from 0-based CSR
 *   amg_solver_set_matrix  - replace solve matrix (A and ichol), keep AMG hierarchy
 *   amg_solver_set_smoother- change smoother type and parameters
 *   amg_solver_solve       - PCG solve with AMG+ichol preconditioner
 *   amg_solver_size        - return matrix dimension
 *   amg_solver_num_levels  - return number of AMG levels
 *   amg_solver_free        - cleanup
 */
#include "amg.h"

typedef struct {
  dcsr_mat     A;
  amg_data     amg;
  ichol_factor L;
} amg_solver;

/* Preconditioner callback matching amg_pcg's signature */
typedef struct {
  const amg_data*     amg;
  const dcsr_mat*     A;
  const ichol_factor* L;
} precond_ctx;

static void precond_wrap(const void* data, const REAL* g, REAL* Bg) {
  const precond_ctx* ctx = (const precond_ctx*)data;
  amg_ichol_precond(ctx->amg, ctx->A, ctx->L, g, Bg);
}

/*
 * amg_solver_create - Build AMG hierarchy and ichol factor
 *
 * Parameters (all 0-based CSR):
 *   n, nnz       - matrix dimensions
 *   ia[n+1]      - row pointers
 *   ja[nnz]      - column indices
 *   val[nnz]     - values
 *   threshold    - strength threshold (e.g. 0.25)
 *   min_size     - coarsening stops below this size
 *   max_levels   - maximum AMG levels
 *   nu           - smoothing sweeps
 *
 * Returns opaque handle (NULL on failure).
 */
void* amg_solver_create(int n, int nnz,
                         const int* ia, const int* ja, const double* val,
                         double threshold, int min_size, int max_levels,
                         int nu) {
  amg_solver* s = (amg_solver*)calloc(1, sizeof(amg_solver));
  if (!s) return NULL;

  /* Copy CSR input into a dcsr_mat */
  dcsr_mat Ain;
  dcsr_mat_alloc(&Ain, n, n, nnz);
  memcpy(Ain.ia, ia, (n + 1) * sizeof(INT));
  memcpy(Ain.ja, ja, nnz * sizeof(INT));
  memcpy(Ain.val, val, nnz * sizeof(REAL));

  /* Symmetrize: A = (Ain + Ain') / 2 */
  dcsr_mat At;
  dcsr_mat_alloc(&At, n, n, nnz);
  dcsr_transz(&Ain, NULL, &At);

  dcsr_mat Asum;
  dcsr_mat_add(&Ain, &At, &Asum);
  for (INT k = 0; k < Asum.nnz; k++) Asum.val[k] *= 0.5;

  dcsr_mat_free(&Ain);
  dcsr_mat_free(&At);
  s->A = Asum;

  /* Build AMG hierarchy */
  memset(&s->amg, 0, sizeof(amg_data));
  s->amg.threshold  = threshold;
  s->amg.min_size   = min_size;
  s->amg.max_levels = (max_levels > 0 && max_levels <= AMG_MAX_LEVELS)
                        ? max_levels : AMG_MAX_LEVELS;
  s->amg.maxit_prol = 0;
  s->amg.lump       = 0;
  s->amg.nu         = (nu > 0) ? nu : 1;
  s->amg.smoother     = SMOOTHER_CF_GS;
  s->amg.omega_jacobi = 0.5;
  s->amg.coarsening   = COARSEN_PMIS;

  amg_build_hierarchy(&s->amg, &s->A);

  /* Incomplete Cholesky */
  ichol_compute(&s->A, &s->L);

  return (void*)s;
}

/*
 * amg_solver_update - Update solver for a spectrally equivalent matrix
 *
 * Reuses the existing AMG hierarchy structure (CF splitting, P sparsity)
 * but recomputes P values, coarse operators (RAP), and ichol factors.
 *
 * Parameters (all 0-based CSR, same sparsity as original):
 *   handle       - from amg_solver_create
 *   n, nnz       - matrix dimensions (must match original)
 *   ia[n+1]      - row pointers
 *   ja[nnz]      - column indices
 *   val[nnz]     - new values
 *
 * Returns 0 on success.
 */
int amg_solver_update(void* handle, int n, int nnz,
                       const int* ia, const int* ja, const double* val) {
  amg_solver* s = (amg_solver*)handle;

  /* Build new CSR matrix */
  dcsr_mat A_new;
  dcsr_mat_alloc(&A_new, n, n, nnz);
  memcpy(A_new.ia, ia, (n + 1) * sizeof(INT));
  memcpy(A_new.ja, ja, nnz * sizeof(INT));
  memcpy(A_new.val, val, nnz * sizeof(REAL));

  /* Symmetrize: A = (A_new + A_new') / 2 */
  dcsr_mat At;
  dcsr_mat_alloc(&At, n, n, nnz);
  dcsr_transz(&A_new, NULL, &At);

  dcsr_mat Asum;
  dcsr_mat_add(&A_new, &At, &Asum);
  for (INT k = 0; k < Asum.nnz; k++) Asum.val[k] *= 0.5;

  dcsr_mat_free(&A_new);
  dcsr_mat_free(&At);

  /* Free old fine-level ichol and matrix */
  ichol_free(&s->L);
  dcsr_mat_free(&s->A);
  s->A = Asum;

  /* Rebuild AMG hierarchy values */
  amg_rebuild_values(&s->amg, &s->A);

  /* Recompute fine-level ichol */
  ichol_compute(&s->A, &s->L);

  return 0;
}

/*
 * amg_solver_set_matrix - Replace solve matrix A and ichol, keep AMG hierarchy
 *
 * Use this for two-matrix preconditioning: build AMG from one matrix (Ap),
 * then set a different matrix (A) for the actual solve.  The AMG hierarchy
 * (CF splitting, interpolation, coarse operators) is NOT modified.
 * The two matrices only need the same dimensions, not the same sparsity.
 *
 * Parameters (all 0-based CSR):
 *   handle       - from amg_solver_create
 *   n, nnz       - matrix dimensions
 *   ia[n+1]      - row pointers
 *   ja[nnz]      - column indices
 *   val[nnz]     - values
 *
 * Returns 0 on success.
 */
int amg_solver_set_matrix(void* handle, int n, int nnz,
                           const int* ia, const int* ja, const double* val) {
  amg_solver* s = (amg_solver*)handle;

  /* Build new CSR matrix */
  dcsr_mat A_new;
  dcsr_mat_alloc(&A_new, n, n, nnz);
  memcpy(A_new.ia, ia, (n + 1) * sizeof(INT));
  memcpy(A_new.ja, ja, nnz * sizeof(INT));
  memcpy(A_new.val, val, nnz * sizeof(REAL));

  /* Symmetrize: A = (A_new + A_new') / 2 */
  dcsr_mat At;
  dcsr_mat_alloc(&At, n, n, nnz);
  dcsr_transz(&A_new, NULL, &At);

  dcsr_mat Asum;
  dcsr_mat_add(&A_new, &At, &Asum);
  for (INT k = 0; k < Asum.nnz; k++) Asum.val[k] *= 0.5;

  dcsr_mat_free(&A_new);
  dcsr_mat_free(&At);

  /* Replace A and ichol, leave AMG hierarchy untouched */
  ichol_free(&s->L);
  dcsr_mat_free(&s->A);
  s->A = Asum;

  ichol_compute(&s->A, &s->L);

  return 0;
}

/*
 * amg_solver_solve - Solve A*x = b using PCG with AMG+ichol
 *
 * Parameters:
 *   handle       - from amg_solver_create
 *   b[n]         - right-hand side
 *   x[n]         - solution (zeroed on entry recommended)
 *   tol          - relative residual tolerance
 *   maxiter      - maximum CG iterations
 *   print_level  - 0 = silent, 1 = per-iteration residuals
 *
 * Returns iteration count (negative if not converged).
 */
int amg_solver_solve(void* handle,
                     const double* b, double* x,
                     double tol, int maxiter, int print_level) {
  amg_solver* s = (amg_solver*)handle;

  precond_ctx ctx;
  ctx.amg = &s->amg;
  ctx.A   = &s->A;
  ctx.L   = &s->L;

  return amg_pcg(&s->A, b, x, precond_wrap, &ctx, tol, maxiter, print_level);
}

/*
 * amg_solver_set_smoother - Change smoother type and parameters
 *
 * Smoother types:
 *   0 = L1-Jacobi (parameter-free, robust default)
 *   1 = lexicographic Gauss-Seidel
 *   2 = CF-ordered Gauss-Seidel (default)
 *   3 = damped Jacobi (uses omega)
 *
 * Parameters:
 *   handle       - from amg_solver_create
 *   smoother     - smoother type (0-3)
 *   omega        - Jacobi damping factor (only used when smoother=3)
 */
void amg_solver_set_smoother(void* handle, int smoother, double omega) {
  amg_solver* s = (amg_solver*)handle;
  s->amg.smoother = smoother;
  if (smoother == SMOOTHER_JACOBI)
    s->amg.omega_jacobi = omega;
}

/* Return matrix dimension */
int amg_solver_size(void* handle) {
  amg_solver* s = (amg_solver*)handle;
  return s->A.row;
}

/* Return number of AMG levels */
int amg_solver_num_levels(void* handle) {
  amg_solver* s = (amg_solver*)handle;
  return s->amg.num_levels;
}

/* Free all resources */
void amg_solver_free(void* handle) {
  if (!handle) return;
  amg_solver* s = (amg_solver*)handle;
  ichol_free(&s->L);
  amg_free(&s->amg);
  dcsr_mat_free(&s->A);
  free(s);
}
