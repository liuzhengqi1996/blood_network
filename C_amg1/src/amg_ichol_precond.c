/*
 * amg_ichol_precond.c - Combined AMG + ichol preconditioner
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Translated from amg_ichol_precond.m
 *
 *   (1) u1 = S * g              (backslash AMG cycle)
 *   (2) u2 = u1 + L\(L'\(g - A*u1))  (ichol correction)
 *   (3) x  = u2 + S^T * (g - A*u2)   (fwdslash AMG cycle)
 */
#include "amg.h"

/* Data bundle passed through void* to PCG */
typedef struct {
  const amg_data* amg;
  const dcsr_mat* A;
  const ichol_factor* L;
} amg_ichol_data;

static void amg_ichol_apply(const void* data, const REAL* g, REAL* x) {
  const amg_ichol_data* d = (const amg_ichol_data*)data;
  INT n = d->amg->levels[0].n;
  INT nu = d->amg->nu;

  REAL* u1 = (REAL*)calloc(n, sizeof(REAL));
  REAL* r1 = (REAL*)calloc(n, sizeof(REAL));
  REAL* ichol_corr = (REAL*)calloc(n, sizeof(REAL));
  REAL* r2 = (REAL*)calloc(n, sizeof(REAL));
  REAL* corr = (REAL*)calloc(n, sizeof(REAL));

  /* Step 1: u1 = backslash(g) */
  amg_backslash(d->amg, 0, g, u1, nu);

  /* Step 2: r1 = g - A*u1; u2 = u1 + ichol_solve(r1) */
  abyvam(g, d->A->ia, d->A->ja, d->A->val, u1, n, r1);
  ichol_solve(d->L, r1, ichol_corr);
  /* u2 stored in x: x = u1 + ichol_corr */
  copyv(u1, x, n);
  uuplmv(x, ichol_corr, n, 1.0);

  /* Step 3: r2 = g - A*u2; x = u2 + fwdslash(r2) */
  abyvam(g, d->A->ia, d->A->ja, d->A->val, x, n, r2);
  amg_fwdslash(d->amg, 0, r2, corr, nu);
  uuplmv(x, corr, n, 1.0);

  free(u1); free(r1); free(ichol_corr); free(r2); free(corr);
}

void amg_ichol_precond(const amg_data* amg, const dcsr_mat* A,
                       const ichol_factor* L,
                       const REAL* g, REAL* x) {
  amg_ichol_data d;
  d.amg = amg;
  d.A = A;
  d.L = L;
  amg_ichol_apply(&d, g, x);
}
