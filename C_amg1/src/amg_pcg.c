/*
 * amg_pcg.c - Preconditioned Conjugate Gradient solver
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Generic PCG that takes a function pointer for the preconditioner.
 * Returns number of iterations (negative if not converged).
 */
#include "amg.h"

INT amg_pcg(const dcsr_mat* A,
            const REAL* b, REAL* x,
            void (*precond)(const void* data, const REAL* g, REAL* Bg),
            const void* precond_data,
            REAL tol, INT maxiter, INT print_level) {
  INT n = A->row;

  REAL* r   = (REAL*)calloc(n, sizeof(REAL));
  REAL* z   = (REAL*)calloc(n, sizeof(REAL));
  REAL* p   = (REAL*)calloc(n, sizeof(REAL));
  REAL* Ap  = (REAL*)calloc(n, sizeof(REAL));

  /* r = b - A*x */
  abyvam(b, A->ia, A->ja, A->val, x, n, r);

  REAL bnorm;
  l2norm(b, &bnorm, n);
  if (bnorm == 0.0) bnorm = 1.0;

  /* z = B*r */
  precond(precond_data, r, z);

  /* p = z */
  copyv(z, p, n);

  REAL rz_old;
  scpro(r, z, &rz_old, n);
  REAL rz0 = rz_old;  /* initial (Br,r) for relative B-norm residual */

  REAL rnorm;
  l2norm(r, &rnorm, n);

  if (print_level >= 1) {
    fprintf(stderr, "  %4s  %8s  %8s  %6s\n",
            "iter", "||r||_B", "rel", "rho");
    fprintf(stderr, "  %4d  %8.2e  %8.2e\n",
            0, sqrt(rz0), 1.0);
  }

  INT iter;
  for (iter = 1; iter <= maxiter; iter++) {
    /* Ap = A*p */
    abyvg(A->ia, A->ja, A->val, p, n, Ap);

    /* alpha = rz / (p'*Ap) */
    REAL pAp;
    scpro(p, Ap, &pAp, n);
    REAL alpha = rz_old / pAp;

    /* x = x + alpha*p */
    uuplmv(x, p, n, alpha);

    /* r = r - alpha*Ap */
    uuplmv(r, Ap, n, -alpha);

    l2norm(r, &rnorm, n);

    /* z = B*r */
    precond(precond_data, r, z);

    REAL rz_new;
    scpro(r, z, &rz_new, n);

    REAL Bnorm_rel = sqrt(rz_new / rz0);
    REAL reduction = rz_new / rz_old;

    if (print_level >= 1) {
      fprintf(stderr, "  %4d  %8.2e  %8.2e  %6.4f\n",
              iter, sqrt(rz_new), Bnorm_rel, reduction);
    }

    if (Bnorm_rel < tol) break;

    REAL beta = rz_new / rz_old;
    rz_old = rz_new;

    /* p = z + beta*p */
    umuplv(p, z, n, beta);
  }

  REAL relres = rnorm / bnorm;
  if (iter > maxiter) {
    fprintf(stderr, "  PCG: no convergence after %d iterations, relres = %.2e\n",
            maxiter, relres);
    iter = -maxiter;
  }

  free(r); free(z); free(p); free(Ap);
  return iter;
}
