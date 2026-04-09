/*
 * AMG + Incomplete Cholesky Preconditioner - Header
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Classical Ruge-Stuben AMG hierarchy with CF-ordered Gauss-Seidel smoothing.
 *
 * Uses dcsr_mat from agg_lib.h (0-based CSR).
 */

#ifndef AMG_H
#define AMG_H

#include "agg_lib.h"
#include <time.h>

#define AMG_MAX_LEVELS 20

/* ---------- Incomplete Cholesky factor (lower triangular CSR) ---------- */
typedef struct {
  INT  n;
  INT* ia;    /* row pointers, size n+1 */
  INT* ja;    /* column indices */
  REAL* val;  /* values */
  INT  nnz;
} ichol_factor;

/* ---------- AMG level structure ---------- */
typedef struct {
  dcsr_mat A;          /* system matrix at this level */
  dcsr_mat P;          /* prolongation (n_fine x n_coarse) */
  dcsr_mat A_filtered; /* strength matrix S */
  INT* mis;            /* mis[i]=1 if C-point, 0 if F-point */
  INT* isolated;       /* isolated[i]=1 if isolated node */
  INT* cf_order;       /* CF ordering: C-points first, then F-points */
  INT n;               /* number of unknowns */
  INT n_coarse;        /* number of coarse points */
  ichol_factor L_ichol;/* ichol factor (coarsest level only) */
} amg_level;

/* ---------- AMG hierarchy ---------- */
typedef struct {
  amg_level levels[AMG_MAX_LEVELS];
  INT num_levels;
  /* options */
  REAL threshold;
  INT  min_size;
  INT  max_levels;
  INT  maxit_prol;
  INT  lump;
  INT  nu;  /* smoothing sweeps */
} amg_data;

/* ======================================================================
 * Function prototypes
 * ====================================================================== */

/* --- hierarchy construction (build_hierarchy.c) --- */
void amg_build_hierarchy(amg_data* amg, const dcsr_mat* A);
void amg_rebuild_values(amg_data* amg, const dcsr_mat* A_new);
void amg_free(amg_data* amg);

/* --- V-cycle (amg_cycle.c) --- */
void amg_backslash(const amg_data* amg, INT lev,
                   const REAL* b, REAL* x, INT nu);
void amg_fwdslash(const amg_data* amg, INT lev,
                  const REAL* b, REAL* x, INT nu);
void amg_vcycle_precond(const amg_data* amg,
                        const REAL* g, REAL* x);

/* --- ichol (ichol.c) --- */
void ichol_compute(const dcsr_mat* A, ichol_factor* L);
void ichol_solve(const ichol_factor* L, const REAL* b, REAL* x);
void ichol_free(ichol_factor* L);

/* --- combined preconditioner (amg_ichol_precond.c) --- */
void amg_ichol_precond(const amg_data* amg, const dcsr_mat* A,
                       const ichol_factor* L,
                       const REAL* g, REAL* x);

/* --- classical RS AMG (rs_classical.c) --- */
void rs_strength(const dcsr_mat* A, REAL theta, dcsr_mat* S, INT* max_lam);
void rs_coarsening(const dcsr_mat* S, INT max_lam, INT* cf, INT* n_coarse);
void rs_standard_interpolation_sparsity(const dcsr_mat* S,
                                        const INT* cf, dcsr_mat* P);
void rs_standard_interpolation_values(const dcsr_mat* A, const dcsr_mat* S,
                                      const INT* cf, dcsr_mat* P);
void rs_standard_interpolation(const dcsr_mat* A, const dcsr_mat* S,
                               const INT* cf, dcsr_mat* P);

/* --- PCG solver (amg_pcg.c) --- */
INT amg_pcg(const dcsr_mat* A,
            const REAL* b, REAL* x,
            void (*precond)(const void* data, const REAL* g, REAL* Bg),
            const void* precond_data,
            REAL tol, INT maxiter, INT print_level);

#endif /* AMG_H */
