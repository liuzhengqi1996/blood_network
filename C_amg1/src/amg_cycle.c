/*
 * amg_cycle.c - AMG backslash/fwdslash cycles and V-cycle preconditioner
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Smoothers: L1-Jacobi (default), Gauss-Seidel, CF-ordered GS, damped Jacobi.
 * Coarsest level solved by ichol-PCG.
 */
#include "amg.h"

/* ichol preconditioner wrapper for amg_pcg */
static void ichol_precond_coarse(const void* data, const REAL* g, REAL* Bg) {
  const ichol_factor* L = (const ichol_factor*)data;
  ichol_solve(L, g, Bg);
}

/* ======================================================================
 *  Lexicographic Gauss-Seidel
 * ====================================================================== */
static void fwd_gs(const dcsr_mat* A, INT n,
                   const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT i = 0; i < n; i++) {
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        if (A->ja[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->ja[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

static void bwd_gs(const dcsr_mat* A, INT n,
                   const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT i = n - 1; i >= 0; i--) {
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        if (A->ja[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->ja[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

/* ======================================================================
 *  CF-ordered Gauss-Seidel
 * ====================================================================== */
static void cf_fwd_gs(const dcsr_mat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = 0; ii < n; ii++) {
      INT i = cf_order[ii];
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        if (A->ja[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->ja[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

static void cf_bwd_gs(const dcsr_mat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = n - 1; ii >= 0; ii--) {
      INT i = cf_order[ii];
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        if (A->ja[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->ja[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

/* ======================================================================
 *  Damped Jacobi:  x = x + omega * D^{-1} (b - Ax)
 * ====================================================================== */
static void jacobi_smooth(const dcsr_mat* A, INT n,
                          const REAL* b, REAL* x, INT nu, REAL omega) {
  REAL* r = (REAL*)malloc(n * sizeof(REAL));
  for (INT sweep = 0; sweep < nu; sweep++) {
    /* r_i = (b_i - A_i x) / a_ii */
    for (INT i = 0; i < n; i++) {
      REAL s = b[i], diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        if (A->ja[k] == i) diag = A->val[k];
        s -= A->val[k] * x[A->ja[k]];
      }
      r[i] = s / diag;
    }
    for (INT i = 0; i < n; i++)
      x[i] += omega * r[i];
  }
  free(r);
}

/* ======================================================================
 *  L1-Jacobi smoother (Baker, Falgout, Kolev, Yang 2011)
 *
 *    x = x + D_l1^{-1} (b - Ax)
 *    (D_l1)_ii = sum_j |a_ij|  (l1 row norm)
 *
 *  Automatically well-damped: effective omega_i = a_ii / ||a_i||_1.
 *  Convergent for SPD A without tuning parameters.
 * ====================================================================== */
static void l1_smooth(const dcsr_mat* A, const REAL* l1_diag, INT n,
                      const REAL* b, REAL* x, INT nu) {
  REAL* r = (REAL*)malloc(n * sizeof(REAL));
  for (INT sweep = 0; sweep < nu; sweep++) {
    /* r = b - A*x */
    for (INT i = 0; i < n; i++) {
      REAL s = b[i];
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++)
        s -= A->val[k] * x[A->ja[k]];
      r[i] = s;
    }
    /* x = x + D_l1^{-1} * r */
    for (INT i = 0; i < n; i++)
      x[i] += r[i] / l1_diag[i];
  }
  free(r);
}

/* ======================================================================
 *  Smoother dispatch: forward (pre-smooth) and backward (post-smooth)
 *
 *  For GS variants, forward/backward are different orderings.
 *  For Jacobi/L1, both directions apply the same symmetric update.
 * ====================================================================== */
static void smooth_fwd(const amg_data* amg, INT lev,
                       const REAL* b, REAL* x, INT nu) {
  const amg_level* L = &amg->levels[lev];
  switch (amg->smoother) {
    case SMOOTHER_GS:
      fwd_gs(&L->A, L->n, b, x, nu);
      break;
    case SMOOTHER_CF_GS:
      cf_fwd_gs(&L->A, L->cf_order, L->n, b, x, nu);
      break;
    case SMOOTHER_JACOBI:
      jacobi_smooth(&L->A, L->n, b, x, nu, amg->omega_jacobi);
      break;
    case SMOOTHER_L1:
    default:
      l1_smooth(&L->A, L->l1_diag, L->n, b, x, nu);
      break;
  }
}

static void smooth_bwd(const amg_data* amg, INT lev,
                       const REAL* b, REAL* x, INT nu) {
  const amg_level* L = &amg->levels[lev];
  switch (amg->smoother) {
    case SMOOTHER_GS:
      bwd_gs(&L->A, L->n, b, x, nu);
      break;
    case SMOOTHER_CF_GS:
      cf_bwd_gs(&L->A, L->cf_order, L->n, b, x, nu);
      break;
    case SMOOTHER_JACOBI:
      jacobi_smooth(&L->A, L->n, b, x, nu, amg->omega_jacobi);
      break;
    case SMOOTHER_L1:
    default:
      l1_smooth(&L->A, L->l1_diag, L->n, b, x, nu);
      break;
  }
}

/* ======================================================================
 *  Backslash recursive: pre-smooth + coarse correct
 * ====================================================================== */
static void backslash_rec(const amg_data* amg, INT lev,
                          const REAL* b, REAL* x, INT nu) {
  INT nlev = amg->num_levels;
  const amg_level* L = &amg->levels[lev];
  INT n = L->n;

  if (lev == nlev - 1) {
    /* Coarsest level: ichol-preconditioned PCG */
    amg_pcg(&L->A, b, x, ichol_precond_coarse, &L->L_ichol, 1e-12, n, 0);
    return;
  }

  /* Pre-smoothing */
  smooth_fwd(amg, lev, b, x, nu);

  /* Residual: r = b - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  abyvam(b, L->A.ia, L->A.ja, L->A.val, x, n, r);

  /* Restrict: rc = P' * r */
  INT nc = L->n_coarse;
  REAL* rc = (REAL*)calloc(nc, sizeof(REAL));
  vbya(L->P.ia, L->P.ja, L->P.val, r, n, nc, rc);

  /* Solve on coarse grid */
  REAL* ec = (REAL*)calloc(nc, sizeof(REAL));
  backslash_rec(amg, lev + 1, rc, ec, nu);

  /* Prolongate: x = x + P * ec */
  abyvcs(L->P.ia, L->P.ja, L->P.val, ec, n, x);

  free(r); free(rc); free(ec);
}

void amg_backslash(const amg_data* amg, INT lev,
                   const REAL* b, REAL* x, INT nu) {
  backslash_rec(amg, lev, b, x, nu);
}

/* ======================================================================
 *  Fwdslash recursive: coarse correct + post-smooth
 * ====================================================================== */
static void fwdslash_rec(const amg_data* amg, INT lev,
                         const REAL* b, REAL* x, INT nu) {
  INT nlev = amg->num_levels;
  const amg_level* L = &amg->levels[lev];
  INT n = L->n;

  if (lev == nlev - 1) {
    amg_pcg(&L->A, b, x, ichol_precond_coarse, &L->L_ichol, 1e-12, n, 0);
    return;
  }

  /* Residual: r = b - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  abyvam(b, L->A.ia, L->A.ja, L->A.val, x, n, r);

  /* Restrict: rc = P' * r */
  INT nc = L->n_coarse;
  REAL* rc = (REAL*)calloc(nc, sizeof(REAL));
  vbya(L->P.ia, L->P.ja, L->P.val, r, n, nc, rc);

  /* Solve on coarse grid */
  REAL* ec = (REAL*)calloc(nc, sizeof(REAL));
  fwdslash_rec(amg, lev + 1, rc, ec, nu);

  /* Prolongate: x = x + P * ec */
  abyvcs(L->P.ia, L->P.ja, L->P.val, ec, n, x);

  /* Post-smoothing */
  smooth_bwd(amg, lev, b, x, nu);

  free(r); free(rc); free(ec);
}

void amg_fwdslash(const amg_data* amg, INT lev,
                  const REAL* b, REAL* x, INT nu) {
  fwdslash_rec(amg, lev, b, x, nu);
}

/* ======================================================================
 *  Symmetric V-cycle: backslash + fwdslash
 * ====================================================================== */
void amg_vcycle_precond(const amg_data* amg,
                        const REAL* g, REAL* x) {
  INT n = amg->levels[0].n;
  INT nu = amg->nu;

  nullv(x, n);
  amg_backslash(amg, 0, g, x, nu);

  /* r = g - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  abyvam(g, amg->levels[0].A.ia, amg->levels[0].A.ja,
         amg->levels[0].A.val, x, n, r);

  /* fwdslash on residual */
  REAL* corr = (REAL*)calloc(n, sizeof(REAL));
  amg_fwdslash(amg, 0, r, corr, nu);

  /* x = x + corr */
  uuplmv(x, corr, n, 1.0);

  free(r); free(corr);
}
