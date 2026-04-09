/*
 * ichol.c - Incomplete Cholesky factorization (ichol(0))
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Computes L such that A ≈ L*L', keeping the sparsity pattern of
 * the lower triangle of A.
 */
#include "amg.h"

void ichol_compute(const dcsr_mat* A, ichol_factor* L) {
  INT n = A->row;
  L->n = n;

  /* Count lower-triangular entries (including diagonal) */
  INT nnz_lower = 0;
  for (INT i = 0; i < n; i++)
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++)
      if (A->ja[k] <= i) nnz_lower++;

  L->nnz = nnz_lower;
  L->ia = (INT*)calloc(n + 1, sizeof(INT));
  L->ja = (INT*)calloc(nnz_lower, sizeof(INT));
  L->val = (REAL*)calloc(nnz_lower, sizeof(REAL));

  /* Extract lower triangle */
  INT pos = 0;
  for (INT i = 0; i < n; i++) {
    L->ia[i] = pos;
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      INT j = A->ja[k];
      if (j <= i) {
        L->ja[pos] = j;
        L->val[pos] = A->val[k];
        pos++;
      }
    }
  }
  L->ia[n] = pos;

  /* Build column-to-position lookup for fast access:
   * For each row i, find the position of column j in L's row i. */
  /* We need L stored by rows. For ichol(0) we process row by row. */

  /* ichol(0) algorithm:
   * For i = 0..n-1:
   *   For each k in lower-triangle row i (k < i):
   *     L(i,k) = (A(i,k) - sum_{j<k} L(i,j)*L(k,j)) / L(k,k)
   *   L(i,i) = sqrt(A(i,i) - sum_{j<i} L(i,j)^2)
   *
   * But we only update entries in the sparsity pattern.
   */

  /* We need to access L(k, j) for j < k, which means accessing row k.
   * Use a dense work vector indexed by column. */
  REAL* w = (REAL*)calloc(n, sizeof(REAL));
  INT* w_marker = (INT*)calloc(n, sizeof(INT));
  for (INT i = 0; i < n; i++) w_marker[i] = -1;

  /* For diagonal positions */
  INT* diag_pos = (INT*)malloc(n * sizeof(INT));
  for (INT i = 0; i < n; i++) {
    for (INT k = L->ia[i]; k < L->ia[i + 1]; k++) {
      if (L->ja[k] == i) { diag_pos[i] = k; break; }
    }
  }

  for (INT i = 0; i < n; i++) {
    /* Load row i of L into w */
    for (INT k = L->ia[i]; k < L->ia[i + 1]; k++) {
      w[L->ja[k]] = L->val[k];
      w_marker[L->ja[k]] = i;
    }

    /* For each off-diagonal entry L(i,k) with k < i */
    for (INT kk = L->ia[i]; kk < L->ia[i + 1]; kk++) {
      INT k = L->ja[kk];
      if (k >= i) continue;

      /* L(i,k) = (w[k] - sum_{j<k} L(i,j)*L(k,j)) / L(k,k) */
      REAL s = w[k];
      /* Iterate over row k of L for j < k */
      for (INT jj = L->ia[k]; jj < L->ia[k + 1]; jj++) {
        INT j = L->ja[jj];
        if (j >= k) break;
        if (w_marker[j] == i) {
          s -= w[j] * L->val[jj];
        }
      }
      REAL Lkk = L->val[diag_pos[k]];
      w[k] = s / Lkk;
    }

    /* Diagonal: L(i,i) = sqrt(w[i] - sum_{j<i} w[j]^2) */
    REAL diag_val = w[i];
    for (INT kk = L->ia[i]; kk < L->ia[i + 1]; kk++) {
      INT j = L->ja[kk];
      if (j >= i) break;
      diag_val -= w[j] * w[j];
    }
    if (diag_val <= 0.0) diag_val = 1e-14; /* safeguard */
    w[i] = sqrt(diag_val);

    /* Write back to L */
    for (INT kk = L->ia[i]; kk < L->ia[i + 1]; kk++) {
      L->val[kk] = w[L->ja[kk]];
    }

    /* Clear w */
    for (INT kk = L->ia[i]; kk < L->ia[i + 1]; kk++) {
      w[L->ja[kk]] = 0.0;
      w_marker[L->ja[kk]] = -1;
    }
  }

  free(w); free(w_marker); free(diag_pos);
}

/* Solve L * L' * x = b  =>  L*y = b (fwd), L'*x = y (bwd) */
void ichol_solve(const ichol_factor* L, const REAL* b, REAL* x) {
  INT n = L->n;
  REAL* y = (REAL*)calloc(n, sizeof(REAL));

  /* Forward solve: L * y = b */
  for (INT i = 0; i < n; i++) {
    REAL s = b[i];
    REAL diag = 1.0;
    for (INT k = L->ia[i]; k < L->ia[i + 1]; k++) {
      INT j = L->ja[k];
      if (j < i)
        s -= L->val[k] * y[j];
      else if (j == i)
        diag = L->val[k];
    }
    y[i] = s / diag;
  }

  /* Backward solve: L' * x = y
   * L' is upper triangular; process rows in reverse.
   * L'(j,i) = L(i,j) for j >= i.
   * For row j of L': x[j] = (y[j] - sum_{i>j} L(i,j)*x[i]) / L(j,j)
   * We accumulate contributions from L stored by rows.
   */
  copyv(y, x, n);

  /* Find diagonal values */
  REAL* Ldiag = (REAL*)calloc(n, sizeof(REAL));
  for (INT i = 0; i < n; i++) {
    for (INT k = L->ia[i]; k < L->ia[i + 1]; k++) {
      if (L->ja[k] == i) { Ldiag[i] = L->val[k]; break; }
    }
  }

  for (INT i = n - 1; i >= 0; i--) {
    x[i] /= Ldiag[i];
    /* Subtract contribution to earlier rows */
    for (INT k = L->ia[i]; k < L->ia[i + 1]; k++) {
      INT j = L->ja[k];
      if (j < i) x[j] -= L->val[k] * x[i];
    }
  }

  free(y); free(Ldiag);
}

void ichol_free(ichol_factor* L) {
  if (L->ia) { free(L->ia); L->ia = NULL; }
  if (L->ja) { free(L->ja); L->ja = NULL; }
  if (L->val) { free(L->val); L->val = NULL; }
  L->n = 0; L->nnz = 0;
}
