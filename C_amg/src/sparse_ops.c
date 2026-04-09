/*
 * Sparse Matrix Operations
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Translated from Fortran srclib.f
 *
 * Uses 0-based indexing (C style).
 * CSR format: ia[0..n], ja[0..nnz-1], a[0..nnz-1]
 * Row i has entries from ja[ia[i]] to ja[ia[i+1]-1]
 */

#include "agg_lib.h"

/*
 * abyvg - Sparse matrix-vector multiplication: y = A * x
 */
void abyvg(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
           REAL* y) {
  for (INT i = 0; i < n; i++) {
    REAL hold0 = 0.0;
    for (INT ji = ia[i]; ji < ia[i + 1]; ji++) {
      hold0 += a[ji] * x[ja[ji]];
    }
    y[i] = hold0;
  }
}

/*
 * abyvcs - Sparse matrix-vector multiply-add: y = y + A * x
 */
void abyvcs(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
            REAL* y) {
  for (INT i = 0; i < n; i++) {
    REAL hold0 = y[i];
    for (INT ji = ia[i]; ji < ia[i + 1]; ji++) {
      hold0 += a[ji] * x[ja[ji]];
    }
    y[i] = hold0;
  }
}

/*
 * abyvam - Residual calculation: z = x - A * y
 */
void abyvam(const REAL* x, const INT* ia, const INT* ja, const REAL* a,
            const REAL* y, INT n, REAL* z) {
  for (INT i = 0; i < n; i++) {
    REAL hold0 = x[i];
    for (INT ji = ia[i]; ji < ia[i + 1]; ji++) {
      hold0 -= a[ji] * y[ja[ji]];
    }
    z[i] = hold0;
  }
}

/*
 * vbya - Transpose matrix-vector multiplication: y = A^T * x
 */
void vbya(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
          INT m_dim, REAL* y) {
  /* Initialize y to zero */
  for (INT i = 0; i < m_dim; i++) {
    y[i] = 0.0;
  }

  /* Compute y = A^T * x */
  for (INT i = 0; i < n; i++) {
    REAL hold0 = x[i];
    for (INT ji = ia[i]; ji < ia[i + 1]; ji++) {
      y[ja[ji]] += a[ji] * hold0;
    }
  }
}

/*
 * dcsr_cp - Copy a CSR matrix (dst must be pre-allocated)
 */
void dcsr_cp(const dcsr_mat* src, dcsr_mat* dst) {
  dst->row = src->row;
  dst->col = src->col;
  dst->nnz = src->nnz;
  memcpy(dst->ia, src->ia, (src->row + 1) * sizeof(INT));
  memcpy(dst->ja, src->ja, src->nnz * sizeof(INT));
  memcpy(dst->val, src->val, src->nnz * sizeof(REAL));
}

/***********************************************************************************************/
/*!
 * \fn void dcsr_rap(dcsr_mat *r_in, dcsr_mat *a_in, dcsr_mat *p_in, dcsr_mat
 * *rap_out)
 *
 * \brief Triple sparse matrix multiplication rap_out = r_in * a_in * p_in
 * (HAZmath)
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication
 * Package. Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 * \note Index starts at 0!!!
 *
 */
void dcsr_rap(dcsr_mat* r_in, dcsr_mat* a_in, dcsr_mat* p_in,
              dcsr_mat* rap_out) {
  INT nc = r_in->row;
  INT *ria = r_in->ia, *rja = r_in->ja;
  REAL* rval = r_in->val;

  INT nf = a_in->row;
  INT *aia = a_in->ia, *aja = a_in->ja;
  REAL* aval = a_in->val;

  INT *pia = p_in->ia, *pja = p_in->ja;
  REAL* pval = p_in->val;

  INT rap_size, *rap_ia = NULL, *rap_ja = NULL;
  REAL* rap_val = NULL;
  INT *ps_marker = NULL, *as_marker = NULL;

  INT ic, i1, i2, i3, jj1, jj2, jj3, jj_counter, jj_row_beg;
  REAL r_entry, ra_prod, rap_prod;

  INT nc_x_1 = nc, nf_x_1 = nf;
  INT minus_one_len = nc_x_1 + nf_x_1;
  INT total_alloc = minus_one_len + nc + 1 + 1;

  ps_marker = (INT*)calloc(total_alloc, sizeof(INT));
  as_marker = ps_marker + nc_x_1;

  /* First Pass: Determine size of rap and set up rap_ia */
  rap_ia = (INT*)calloc(nc + 1, sizeof(INT));
  for (INT ii = 0; ii < minus_one_len; ++ii) ps_marker[ii] = -1;

  jj_counter = 0;
  for (ic = 0; ic < nc; ic++) {
    ps_marker[ic] = jj_counter;
    jj_row_beg = jj_counter;
    jj_counter++;

    for (jj1 = ria[ic]; jj1 < ria[ic + 1]; jj1++) {
      i1 = rja[jj1];
      for (jj2 = aia[i1]; jj2 < aia[i1 + 1]; jj2++) {
        i2 = aja[jj2];
        if (as_marker[i2] != ic) {
          as_marker[i2] = ic;
          for (jj3 = pia[i2]; jj3 < pia[i2 + 1]; jj3++) {
            i3 = pja[jj3];
            if (ps_marker[i3] < jj_row_beg) {
              ps_marker[i3] = jj_counter;
              jj_counter++;
            }
          }
        }
      }
    }
    rap_ia[ic] = jj_row_beg;
  }

  rap_ia[nc] = jj_counter;
  rap_size = jj_counter;

  rap_ja = (INT*)calloc(rap_size, sizeof(INT));
  rap_val = (REAL*)calloc(rap_size, sizeof(REAL));

  for (INT ii = 0; ii < minus_one_len; ++ii) ps_marker[ii] = -1;

  jj_counter = 0;
  for (ic = 0; ic < nc; ic++) {
    ps_marker[ic] = jj_counter;
    jj_row_beg = jj_counter;
    rap_ja[jj_counter] = ic;
    rap_val[jj_counter] = 0.0;
    jj_counter++;

    for (jj1 = ria[ic]; jj1 < ria[ic + 1]; jj1++) {
      r_entry = rval[jj1];
      i1 = rja[jj1];
      for (jj2 = aia[i1]; jj2 < aia[i1 + 1]; jj2++) {
        ra_prod = r_entry * aval[jj2];
        i2 = aja[jj2];
        if (as_marker[i2] != ic) {
          as_marker[i2] = ic;
          for (jj3 = pia[i2]; jj3 < pia[i2 + 1]; jj3++) {
            rap_prod = ra_prod * pval[jj3];
            i3 = pja[jj3];
            if (ps_marker[i3] < jj_row_beg) {
              ps_marker[i3] = jj_counter;
              rap_val[jj_counter] = rap_prod;
              rap_ja[jj_counter] = i3;
              jj_counter++;
            } else {
              rap_val[ps_marker[i3]] += rap_prod;
            }
          }
        } else {
          for (jj3 = pia[i2]; jj3 < pia[i2 + 1]; jj3++) {
            i3 = pja[jj3];
            rap_prod = ra_prod * pval[jj3];
            rap_val[ps_marker[i3]] += rap_prod;
          }
        }
      }
    }
  }

  rap_out->row = nc;
  rap_out->col = nc;
  rap_out->nnz = rap_size;
  rap_out->ia = rap_ia;
  rap_out->ja = rap_ja;
  rap_out->val = rap_val;

  free(ps_marker);
}

/*********************************EOF***********************************/
