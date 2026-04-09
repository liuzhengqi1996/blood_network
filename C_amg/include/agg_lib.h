/*
 * Aggregation-based AMG Library - Header File
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Translated from Fortran agg2c library
 *
 * This library implements an aggregation-based algebraic multigrid (AMG)
 * solver for sparse linear systems in CSR format.
 *
 * Uses 0-based indexing (C style).
 * CSR format: Row i has entries at ja[ia[i]] to ja[ia[i+1]-1] (exclusive upper
 * bound)
 */

#ifndef AGG_LIB_H
#define AGG_LIB_H

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ========================================================================
 * Basic Types
 * ======================================================================== */
#ifndef INT
#define INT int
#endif
#ifndef REAL
#define REAL double
#endif

/* Sparse matrix in CSR format (double values) */
typedef struct {
  INT row;   /* number of rows */
  INT col;   /* number of columns */
  INT nnz;   /* number of nonzero entries */
  INT* ia;   /* row pointers, size row+1 */
  INT* ja;   /* column indices, size nnz */
  REAL* val; /* nonzero values, size nnz */
} dcsr_mat;

/* Sparse matrix in COO format */
typedef struct {
  INT row;    /* number of rows */
  INT col;    /* number of columns */
  INT nnz;    /* number of nonzero entries */
  INT* icols; /* columns indices size nnz */
  INT* irows; /* row indices indices, size nnz */
  REAL* val;  /* nonzero values, size nnz */
} dcoo_mat;

/* Free dcsr_mat structure */
static inline void dcsr_mat_free(dcsr_mat* a) {
  if (a == NULL) return;
  if (a->ia) {
    free(a->ia);
    a->ia = NULL;
  }
  if (a->ja) {
    free(a->ja);
    a->ja = NULL;
  }
  if (a->val) {
    free(a->val);
    a->val = NULL;
  }
  a->row = 0;
  a->col = 0;
  a->nnz = 0;
}

/* Free dcoo_mat structure */
static inline void dcoo_mat_free(dcoo_mat* a) {
  if (a == NULL) return;
  if (a->irows) {
    free(a->irows);
    a->irows = NULL;
  }
  if (a->icols) {
    free(a->icols);
    a->icols = NULL;
  }
  if (a->val) {
    free(a->val);
    a->val = NULL;
  }
  a->row = 0;
  a->col = 0;
  a->nnz = 0;
}

/* Allocate a dcsr_mat (inline version for generate_sparse_mat.c) */
static inline void dcsr_mat_alloc(dcsr_mat* a, INT nrows, INT ncols, INT nnz) {
  a->row = nrows;
  a->col = ncols;
  a->nnz = nnz;
  a->ia = (INT*)calloc(nrows + 1, sizeof(INT));
  a->ja = (INT*)calloc(nnz, sizeof(INT));
  a->val = (REAL*)calloc(nnz, sizeof(REAL));
}

/*
 * dcsr_transz - Transpose a CSR matrix with optional permutation and optional
 * values being NULL;
 *
 * Computes AT = A^T if p is NULL, or AT = (A^T)*P if p is not NULL.
 * if A_in->val is null then AT also has at->val=NULL.
 * The storage for AT (ia, ja, val) must be pre-allocated.
 * If A->val is NULL, only transposes the structure.
 *
 * Reference: Fred G. Gustavson. Two fast algorithms for sparse
 *            matrices: multiplication and permuted transposition.
 *            ACM Trans. Math. Software, 4(3):250-269, 1978.
 */

static inline void dcsr_transz(const dcsr_mat* a_in, INT* p, dcsr_mat* at_out) {
  const INT n = a_in->row, m = a_in->col, nnz = a_in->nnz;
  at_out->row = m;
  at_out->col = n;
  at_out->nnz = nnz;

  /* Initialize row pointers to zero */
  for (INT i = 0; i <= m; i++) at_out->ia[i] = 0;

  /* Count entries per column of A */
  for (INT k = 0; k < nnz; k++) at_out->ia[a_in->ja[k] + 1]++;

  /* Cumulative sum */
  for (INT i = 1; i <= m; i++) at_out->ia[i] += at_out->ia[i - 1];

  /* Fill in transpose */
  INT* count = (INT*)calloc(m, sizeof(INT));
  if (p == NULL) {
    if (a_in->val != NULL && at_out->val != NULL) {
      for (INT i = 0; i < n; i++) {
        for (INT k = a_in->ia[i]; k < a_in->ia[i + 1]; k++) {
          INT j = a_in->ja[k];
          INT pos = at_out->ia[j] + count[j];
          at_out->ja[pos] = i;
          at_out->val[pos] = a_in->val[k];
          count[j]++;
        }
      }
    } else {
      for (INT i = 0; i < n; i++) {
        for (INT k = a_in->ia[i]; k < a_in->ia[i + 1]; k++) {
          INT j = a_in->ja[k];
          INT pos = at_out->ia[j] + count[j];
          at_out->ja[pos] = i;
          count[j]++;
        }
      }
    }
  } else {
    if (a_in->val != NULL && at_out->val != NULL) {
      for (INT i = 0; i < n; i++) {
        INT pi = p[i];
        for (INT k = a_in->ia[pi]; k < a_in->ia[pi + 1]; k++) {
          INT j = a_in->ja[k];
          INT pos = at_out->ia[j] + count[j];
          at_out->ja[pos] = i;
          at_out->val[pos] = a_in->val[k];
          count[j]++;
        }
      }
    } else {
      for (INT i = 0; i < n; i++) {
        INT pi = p[i];
        for (INT k = a_in->ia[pi]; k < a_in->ia[pi + 1]; k++) {
          INT j = a_in->ja[k];
          INT pos = at_out->ia[j] + count[j];
          at_out->ja[pos] = i;
          count[j]++;
        }
      }
    }
  }
  free(count);
}

/* Sparse matrix addition: C = A + B (inline) */
static inline INT dcsr_mat_add(dcsr_mat * a, dcsr_mat * b, dcsr_mat * c) {
  if (a->row != b->row || a->col != b->col) return -1;

  INT n = a->row, m = a->col;
  INT* marker = (INT*)calloc(m, sizeof(INT));
  for (INT i = 0; i < m; i++) marker[i] = -1;

  /* First pass: count nnz */
  INT nnz = 0;
  for (INT i = 0; i < n; i++) {
    for (INT j = a->ia[i]; j < a->ia[i + 1]; j++) {
      marker[a->ja[j]] = i;
      nnz++;
    }
    for (INT j = b->ia[i]; j < b->ia[i + 1]; j++) {
      if (marker[b->ja[j]] != i) {
        marker[b->ja[j]] = i;
        nnz++;
      }
    }
  }

  dcsr_mat_alloc(c, n, m, nnz);
  for (INT i = 0; i < m; i++) marker[i] = -1;

  REAL* w = (REAL*)calloc(m, sizeof(REAL));
  INT pos = 0;
  for (INT i = 0; i < n; i++) {
    c->ia[i] = pos;
    for (INT j = a->ia[i]; j < a->ia[i + 1]; j++) {
      INT col = a->ja[j];
      marker[col] = pos;
      c->ja[pos] = col;
      w[col] = a->val[j];
      pos++;
    }
    for (INT j = b->ia[i]; j < b->ia[i + 1]; j++) {
      INT col = b->ja[j];
      if (marker[col] >= c->ia[i]) {
        w[col] += b->val[j];
      } else {
        marker[col] = pos;
        c->ja[pos] = col;
        w[col] = b->val[j];
        pos++;
      }
    }
    for (INT j = c->ia[i]; j < pos; j++) {
      c->val[j] = w[c->ja[j]];
      w[c->ja[j]] = 0.0;
    }
  }
  c->ia[n] = pos;
  c->nnz = pos;
  free(marker);
  free(w);
  return 0;
}

/* ========================================================================
 * Vector Operations (from srclib.f)
 * All operations use 0-based indexing: arrays indexed from 0 to n-1
 * ======================================================================== */

/* Zero a REAL precision vector: u[0..n-1] = 0 */
void nullv(REAL * u, INT n);

/* Copy a REAL precision vector: v[0..n-1] = u[0..n-1] */
void copyv(const REAL* u, REAL* v, INT n);

/* Inner product: scpr = u[0..n-1]^T * v[0..n-1] */
void scpro(const REAL* u, const REAL* v, REAL* scpr, INT n);

/* L2 norm: l2nr = ||u[0..n-1]||_2 */
void l2norm(const REAL* u, REAL* l2nr, INT n);

/* u[0..n-1] = u[0..n-1] + smult * v[0..n-1] */
void uuplmv(REAL * u, const REAL* v, INT n, REAL smult);

/* u[0..n-1] = smult * u[0..n-1] + v[0..n-1] */
void umuplv(REAL * u, const REAL* v, INT n, REAL smult);

/* Generate random vector with entries in [0, 1] */
void randv(REAL * u, INT n);

/* ========================================================================
 * CSR Matrix Helper Functions
 * ======================================================================== */

/* Copy a CSR matrix (dst must be pre-allocated) */
void dcsr_cp(const dcsr_mat* src, dcsr_mat* dst);

/* ========================================================================
 * Sparse Matrix Operations (from srclib.f)
 * CSR format: Row i has entries at ja[ia[i]] to ja[ia[i+1]-1]
 * ======================================================================== */

/* Sparse matrix-vector multiplication: y = A * x */
void abyvg(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
           REAL* y);

/* Sparse matrix-vector multiply-add: y = y + A * x */
void abyvcs(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
            REAL* y);

/* Residual calculation: z = x - A * y */
void abyvam(const REAL* x, const INT* ia, const INT* ja, const REAL* a,
            const REAL* y, INT n, REAL* z);

/* Transpose matrix-vector multiplication: y = x^T * A */
void vbya(const INT* ia, const INT* ja, const REAL* a, const REAL* x, INT n,
          INT m_dim, REAL* y);

/* Triple sparse matrix multiplication: rap = r * a * p */
void dcsr_rap(dcsr_mat *r_in, dcsr_mat *a_in, dcsr_mat *p_in, dcsr_mat *rap_out);

/* ========================================================================
 * I/O Functions (from solver.f and readija.f)
 * ======================================================================== */

/* Convert COO to CSR format */
void dcoo2dcsr(dcoo_mat * acoo, dcsr_mat * a);

/* Read matrix from file in COO format and convert to CSR */
INT dcsr_mat_dcoo(FILE * fp, dcsr_mat * a_out);

#endif /* AGG_LIB_H */
