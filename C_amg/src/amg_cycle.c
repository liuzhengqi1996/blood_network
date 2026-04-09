/*
 * amg_cycle.c - AMG backslash/fwdslash cycles and V-cycle preconditioner
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Uses CF-ordered GS smoothing, coarsest level solved by ichol-PCG.
 */
#include "amg.h"

/* ichol preconditioner wrapper for amg_pcg */
static void ichol_precond_coarse(const void* data, const REAL* g, REAL* Bg) {
  const ichol_factor* L = (const ichol_factor*)data;
  ichol_solve(L, g, Bg);
}

/* Forward GS with CF ordering on A*x = b */
static void cf_fwd_gs(const dcsr_mat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = 0; ii < n; ii++) {
      INT i = cf_order[ii];
      REAL rhs = b[i];
      REAL diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        INT j = A->ja[k];
        if (j == i)
          diag = A->val[k];
        else
          rhs -= A->val[k] * x[j];
      }
      x[i] = rhs / diag;
    }
  }
}

/* Backward GS with CF ordering on A*x = b */
static void cf_bwd_gs(const dcsr_mat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = n - 1; ii >= 0; ii--) {
      INT i = cf_order[ii];
      REAL rhs = b[i];
      REAL diag = 1.0;
      for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
        INT j = A->ja[k];
        if (j == i)
          diag = A->val[k];
        else
          rhs -= A->val[k] * x[j];
      }
      x[i] = rhs / diag;
    }
  }
}

/* ====================================================================== */
/*  Backslash recursive: pre-smooth + coarse correct                      */
/* ====================================================================== */
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

  /* Pre-smoothing: forward GS with CF ordering */
  cf_fwd_gs(&L->A, L->cf_order, n, b, x, nu);

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

/* ====================================================================== */
/*  Fwdslash recursive: coarse correct + post-smooth                      */
/* ====================================================================== */
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

  /* Post-smoothing: backward GS with CF ordering */
  cf_bwd_gs(&L->A, L->cf_order, n, b, x, nu);

  free(r); free(rc); free(ec);
}

void amg_fwdslash(const amg_data* amg, INT lev,
                  const REAL* b, REAL* x, INT nu) {
  fwdslash_rec(amg, lev, b, x, nu);
}

/* ====================================================================== */
/*  Symmetric V-cycle: backslash + fwdslash                               */
/* ====================================================================== */
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
