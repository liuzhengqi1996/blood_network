/*
 * Matrix I/O Functions
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Translated from Fortran readija.f and solver.f
 *
 * Uses 0-based indexing (C style).
 */

#include "agg_lib.h"

/*
 * dcoo2dcsr - Convert COO to CSR format
 *
 * Converts from coordinate format to CSR format.
 * All arrays use 0-based indexing.
 */
void dcoo2dcsr(dcoo_mat* acoo, dcsr_mat* a_out) {
  const INT m = acoo->row, nnz = acoo->nnz;

  /* Allocate output arrays */
  a_out->row = acoo->row;
  a_out->col = acoo->col;
  a_out->nnz = acoo->nnz;
  a_out->ia = (INT*)calloc(a_out->row + 1, sizeof(INT));
  a_out->ja = (INT*)calloc(a_out->nnz, sizeof(INT));
  a_out->val = (REAL*)calloc(a_out->nnz, sizeof(REAL));

  INT* ia = a_out->ia;
  INT* ja = a_out->ja;
  REAL* a = a_out->val;
  INT* row_idx = acoo->irows;
  INT* col_idx = acoo->icols;
  REAL* aij = acoo->val;

  INT* ind = (INT*)calloc(m + 1, sizeof(INT));

  /* Count number of nonzeros in each row */
  for (INT i = 0; i < nnz; ++i) ind[row_idx[i] + 1]++;

  /* Set row pointer */
  ia[0] = 0;
  for (INT i = 1; i <= m; ++i) {
    ia[i] = ia[i - 1] + ind[i];
    ind[i] = ia[i];
  }

  /* Set column index and values */
  for (INT i = 0; i < nnz; ++i) {
    INT iind = row_idx[i];
    INT jind = ind[iind];
    ja[jind] = col_idx[i];
    a[jind] = aij[i];
    ind[iind] = ++jind;
  }

  free(ind);
}

/*
 * dcsr_mat_dcoo - Read a matrix from file in COO format and convert to CSR
 *
 * File format:
 *   nrows ncols nnz
 *   row col value    (repeated nnz times)
 *
 * Indices can be 0-based or 1-based (auto-detected and converted to 0-based)
 */
INT dcsr_mat_dcoo(FILE* fp, dcsr_mat* a_out) {
  INT nnz;

  /* Read COO matrix header */
  dcoo_mat acoo;
  if (fscanf(fp, "%d %d %d", &acoo.row, &acoo.col, &acoo.nnz) != 3) {
    fprintf(stderr, "Error reading matrix dimensions\n");
    return 1;
  }
  nnz = acoo.nnz;
  fprintf(stdout,"%%%%Degrees of freedom: %d, Nonzeros: %d\n", acoo.row,
         acoo.nnz);
  fflush(stdout);

  acoo.irows = (INT*)calloc(nnz, sizeof(INT));
  acoo.icols = (INT*)calloc(nnz, sizeof(INT));
  acoo.val = (REAL*)calloc(nnz, sizeof(REAL));

  if (acoo.irows == NULL || acoo.icols == NULL || acoo.val == NULL) {
    fprintf(stderr, "Error allocating COO arrays\n");
    return 1;
  }

  for (INT k = 0; k < nnz; k++) {
    if (fscanf(fp, "%d %d %lf", &acoo.irows[k], &acoo.icols[k], &acoo.val[k]) !=
        3) {
      fprintf(stderr, "Error reading matrix entry %d\n", k);
      dcoo_mat_free(&acoo);
      return 1;
    }
  }

  /* Check if indices need shifting (1-based to 0-based) */
  INT min_idx = acoo.irows[0];
  for (INT k = 0; k < nnz; k++) {
    if (acoo.irows[k] < min_idx) min_idx = acoo.irows[k];
    if (acoo.icols[k] < min_idx) min_idx = acoo.icols[k];
  }

  if (min_idx == 1) {
    for (INT k = 0; k < nnz; k++) {
      acoo.irows[k]--;
      acoo.icols[k]--;
    }
  }

  /* Convert COO to CSR */
  dcoo2dcsr(&acoo, a_out);

  /* Free COO matrix */
  dcoo_mat_free(&acoo);

  return 0;
}
