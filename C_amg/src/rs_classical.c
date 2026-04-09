/*
 * rs_classical.c — Classical Ruge-Stuben AMG Components
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 *   1. Strength of Connection  (RS criterion)
 *   2. C/F Splitting           (classical RS coarsening, first + second pass)
 *   3. Standard Interpolation  (direct interpolation)
 *
 * Uses dcsr_mat from agg_lib.h (0-based CSR).
 */
#include "amg.h"

/* C/F markers */
#define RS_C_PT    1
#define RS_F_PT   (-1)
#define RS_UNDECIDED 0

/* ======================================================================
 *  1. Strength of Connection
 *
 *  For row i of A, off-diagonal entry j is "strong" if
 *
 *      -a_{ij}  >=  theta * max_{k != i} ( -a_{ik} )
 *
 *  Output S has the same dimension as A and stores only strong
 *  off-diagonal entries with their original values from A.
 * ====================================================================== */
void rs_strength(const dcsr_mat *A, REAL theta, dcsr_mat *S, INT *max_lam)
{
  INT n = A->row;

  /* --- pass 1: row-wise max of (-a_{ij}) for j != i --- */
  REAL *row_max = (REAL *)calloc(n, sizeof(REAL));
  INT ml = 0;
  for (INT i = 0; i < n; i++) {
    REAL mx = 0.0;
    INT lambda_i = (INT)(A->ia[i + 1] - A->ia[i]);
    if (lambda_i > ml) ml = lambda_i;
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      if (A->ja[k] == i) continue;
      REAL v = -(A->val[k]);
      if (v > mx) mx = v;
    }
    row_max[i] = mx;
  }
  if (ml < 1) ml = 1;
  *max_lam = ml;

  /* --- pass 2: count strong connections per row --- */
  INT *s_ia = (INT *)calloc(n + 1, sizeof(INT));
  for (INT i = 0; i < n; i++) {
    REAL thr = theta * row_max[i];
    INT cnt = 0;
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      INT j = A->ja[k];
      if (j != i && -(A->val[k]) >= thr)
        cnt++;
    }
    s_ia[i + 1] = cnt;
  }
  for (INT i = 0; i < n; i++)
    s_ia[i + 1] += s_ia[i];

  INT s_nnz = s_ia[n];
  INT  *s_ja  = (INT  *)malloc(s_nnz * sizeof(INT));
  REAL *s_val = (REAL *)malloc(s_nnz * sizeof(REAL));

  /* --- pass 3: fill --- */
  for (INT i = 0; i < n; i++) {
    REAL thr = theta * row_max[i];
    INT p = s_ia[i];
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      INT j = A->ja[k];
      if (j != i && -(A->val[k]) >= thr) {
        s_ja[p]  = j;
        s_val[p] = A->val[k];
        p++;
      }
    }
  }

  free(row_max);

  S->row = n;  S->col = n;  S->nnz = s_nnz;
  S->ia  = s_ia;  S->ja = s_ja;  S->val = s_val;
}

/* ======================================================================
 *  2. C/F Splitting — Classical RS Coarsening
 *
 *  First pass (greedy independent-set style):
 *    lambda_i = |{ j undecided : i in S_j }|  (influence measure)
 *    Repeat:
 *      pick i with max lambda  ->  C-point
 *      for each undecided j that i influences (j in S_i^T):
 *        mark j as F
 *        for each undecided k in S_j (k != i): lambda_k++
 *
 *  Second pass:
 *    Any F-point without a strong C-neighbor is promoted to C.
 *
 *  Output:
 *    cf[i] = RS_C_PT (1) or RS_F_PT (-1), length n.
 *    *n_coarse = number of C-points.
 * ====================================================================== */
void rs_coarsening(const dcsr_mat *S, INT max_lam, INT *cf, INT *n_coarse)
{
  INT n = S->row;

  /* --- transpose S: ST[i] = { j : i in S_j } = points influenced by i --- */
  dcsr_mat ST;
  dcsr_mat_alloc(&ST, n, n, S->nnz);
  dcsr_transz(S, NULL, &ST);

  /* --- initialize lambda and cf --- */
  INT *lambda = (INT *)calloc(n, sizeof(INT));
  cf[0]=RS_UNDECIDED;
  for (INT i = 1; i < n; i++) {
    lambda[i] = ST.ia[i + 1] - ST.ia[i];
    cf[i] = RS_UNDECIDED;
  }

  /* --- bucket structure for O(1) max-lambda lookup --- */
  INT num_buckets = 2 * max_lam + 1;
  if (num_buckets < 1) num_buckets = 1;

  INT *bucket_head = (INT *)malloc(num_buckets * sizeof(INT));
  INT *next = (INT *)malloc(n * sizeof(INT));
  INT *prev = (INT *)malloc(n * sizeof(INT));

  for (INT k = 0; k < num_buckets; k++)
    bucket_head[k] = -1;

  /* Insert all points into their buckets */
  for (INT i = n - 1; i >= 0; i--) {
    INT b = lambda[i];
    next[i] = bucket_head[b];
    prev[i] = -1;
    if (bucket_head[b] >= 0)
      prev[bucket_head[b]] = i;
    bucket_head[b] = i;
  }

  INT top_bucket = max_lam;

  /* --- first pass: greedy selection using buckets, O(n + nnz(S)) --- */
  for (;;) {
    /* Find highest non-empty bucket */
    while (top_bucket >= 0 && bucket_head[top_bucket] < 0)
      top_bucket--;
    if (top_bucket < 0) break;

    /* Remove best from its bucket */
    INT best = bucket_head[top_bucket];
    {
      INT b = lambda[best];
      if (prev[best] >= 0) next[prev[best]] = next[best];
      else                  bucket_head[b] = next[best];
      if (next[best] >= 0) prev[next[best]] = prev[best];
    }
    cf[best] = RS_C_PT;

    /* For each undecided j influenced by best (row best of S^T) */
    for (INT kk = ST.ia[best]; kk < ST.ia[best + 1]; kk++) {
      INT j = ST.ja[kk];
      if (cf[j] != RS_UNDECIDED) continue;

      /* Remove j from its bucket */
      {
        INT b = lambda[j];
        if (prev[j] >= 0) next[prev[j]] = next[j];
        else               bucket_head[b] = next[j];
        if (next[j] >= 0) prev[next[j]] = prev[j];
      }
      cf[j] = RS_F_PT;

      /* j just became F => undecided strong dependencies of j gain value */
      for (INT ll = S->ia[j]; ll < S->ia[j + 1]; ll++) {
        INT k = S->ja[ll];
        if (cf[k] != RS_UNDECIDED) continue;

        /* Remove k from its current bucket */
        {
          INT b = lambda[k];
          if (prev[k] >= 0) next[prev[k]] = next[k];
          else               bucket_head[b] = next[k];
          if (next[k] >= 0) prev[next[k]] = prev[k];
        }

        lambda[k]++;
        INT new_b = lambda[k];

        /* Insert k into its new bucket */
        next[k] = bucket_head[new_b];
        prev[k] = -1;
        if (bucket_head[new_b] >= 0)
          prev[bucket_head[new_b]] = k;
        bucket_head[new_b] = k;

        if (new_b > top_bucket) top_bucket = new_b;
      }
    }
  }

  free(bucket_head);
  free(next);
  free(prev);

  /* --- second pass: every F-point must have >= 1 strong C-neighbor --- */
  for (INT i = 0; i < n; i++) {
    if (cf[i] != RS_F_PT) continue;
    INT has_c = 0;
    for (INT k = S->ia[i]; k < S->ia[i + 1]; k++) {
      if (cf[S->ja[k]] == RS_C_PT) { has_c = 1; break; }
    }
    if (!has_c) cf[i] = RS_C_PT;
  }

  /* count coarse points */
  INT nc = 0;
  for (INT i = 0; i < n; i++)
    if (cf[i] == RS_C_PT) nc++;
  *n_coarse = nc;

  free(lambda);
  dcsr_mat_free(&ST);
}

/* ======================================================================
 *  3. Standard (Direct) Interpolation
 *
 *  For C-point i:  row i of P has a single 1 at coarse index of i.
 *
 *  For F-point i:
 *    C_i   = { j in S_i : cf[j] = C }        strong C-neighbors
 *    D_i^s = { k in S_i : cf[k] = F }        strong F-neighbors
 *    D_i^w = { m in N_i \ S_i : m != i }     weak connections
 *
 *    a_hat_ii = a_{ii} + sum_{m in D_i^w} a_{im}
 *
 *    For each j in C_i:
 *
 *        w_{ij} = -1 / a_hat_ii * ( a_{ij}
 *                   + sum_{k in D_i^s}  a_{ik} * a_{kj}
 *                                      / sum_{l in C_i} a_{kl} )
 *
 *    If a strong F-neighbor k has no connection to any l in C_i,
 *    a_{ik} is lumped into a_hat_ii instead.
 *
 *  Output P has dimension (n x nc) in CSR format.
 * ====================================================================== */
/* ----------------------------------------------------------------------
 *  3a. Build sparsity pattern of P
 *
 *  Allocates and fills P->ia and P->ja (CSR structure).
 *  P->val is allocated and zeroed but not filled with weights.
 *  Only needs S and cf (not A).
 * ---------------------------------------------------------------------- */
void rs_standard_interpolation_sparsity(const dcsr_mat *S,
                                        const INT      *cf,
                                        dcsr_mat       *P)
{
  INT n = S->row;

  /* --- coarse-grid index map --- */
  INT nc = 0;
  INT *cmap = (INT *)malloc(n * sizeof(INT));
  for (INT i = 0; i < n; i++)
    cmap[i] = (cf[i] == RS_C_PT) ? nc++ : -1;

  /* --- count P entries per row --- */
  INT *P_ia = (INT *)calloc(n + 1, sizeof(INT));
  for (INT i = 0; i < n; i++) {
    if (cf[i] == RS_C_PT) {
      P_ia[i + 1] = 1;
    } else {
      INT cnt = 0;
      for (INT k = S->ia[i]; k < S->ia[i + 1]; k++)
        if (cf[S->ja[k]] == RS_C_PT) cnt++;
      P_ia[i + 1] = cnt;
    }
  }
  for (INT i = 0; i < n; i++)
    P_ia[i + 1] += P_ia[i];

  INT P_nnz  = P_ia[n];
  INT  *P_ja = (INT  *)malloc(P_nnz * sizeof(INT));
  REAL *P_v  = (REAL *)calloc(P_nnz, sizeof(REAL));

  /* --- fill column indices --- */
  for (INT i = 0; i < n; i++) {
    if (cf[i] == RS_C_PT) {
      P_ja[P_ia[i]] = cmap[i];
    } else {
      INT p = P_ia[i];
      for (INT k = S->ia[i]; k < S->ia[i + 1]; k++) {
        INT j = S->ja[k];
        if (cf[j] == RS_C_PT) {
          P_ja[p] = cmap[j];
          p++;
        }
      }
    }
  }

  free(cmap);

  P->row = n;  P->col = nc;  P->nnz = P_nnz;
  P->ia  = P_ia;  P->ja = P_ja;  P->val = P_v;
}

/* ----------------------------------------------------------------------
 *  3b. Compute interpolation weights (fill P->val)
 *
 *  Assumes P already has its sparsity pattern set (ia, ja allocated).
 *  For F-point rows, column order in P->ja must match the order of
 *  C-points in S_i (as produced by rs_standard_interpolation_sparsity).
 * ---------------------------------------------------------------------- */
void rs_standard_interpolation_values(const dcsr_mat *A,
                                      const dcsr_mat *S,
                                      const INT      *cf,
                                      dcsr_mat       *P)
{
  INT n = A->row;

  /* --- workspace arrays (size n, reused per row) --- */
  /* is_strong[j] == i  iff  j is in S_i  (strong connection of row i) */
  INT *is_strong = (INT *)malloc(n * sizeof(INT));
  /* pos[j] = offset into P's value array for the C_i entry of j */
  INT *pos = (INT *)malloc(n * sizeof(INT));
  memset(is_strong, -1, n * sizeof(INT));
  memset(pos,       -1, n * sizeof(INT));

  for (INT i = 0; i < n; i++) {

    /* ---- C-point: identity injection ---- */
    if (cf[i] == RS_C_PT) {
      P->val[P->ia[i]] = 1.0;
      continue;
    }

    /* ---- F-point: standard interpolation ---- */

    /* Mark strong connections of row i */
    for (INT k = S->ia[i]; k < S->ia[i + 1]; k++)
      is_strong[S->ja[k]] = i;

    /* Record positions of C_i entries in P */
    {
      INT p = P->ia[i];
      for (INT k = S->ia[i]; k < S->ia[i + 1]; k++) {
        INT j = S->ja[k];
        if (cf[j] == RS_C_PT) {
          P->val[p] = 0.0;
          pos[j]    = p;
          p++;
        }
      }
    }

    /* Walk row i of A to classify each entry */
    REAL diag = 0.0, weak_sum = 0.0;
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      INT  j = A->ja[k];
      REAL a = A->val[k];
      if (j == i) {
        diag = a;                       /* diagonal */
      } else if (is_strong[j] == i && cf[j] == RS_C_PT) {
        P->val[pos[j]] += a;           /* direct a_{ij} for j in C_i */
      } else if (is_strong[j] == i && cf[j] == RS_F_PT) {
        /* strong F-neighbor: handled below */
      } else {
        weak_sum += a;                  /* weak connection -> lump */
      }
    }
    REAL diag_hat = diag + weak_sum;

    /* Distribute strong F-neighbor contributions */
    for (INT k = A->ia[i]; k < A->ia[i + 1]; k++) {
      INT kk = A->ja[k];
      if (kk == i || is_strong[kk] != i || cf[kk] != RS_F_PT)
        continue;

      REAL a_ik = A->val[k];

      /* denominator: sum_{l in C_i} a_{kk,l}  (row kk of A, cols in C_i) */
      REAL denom = 0.0;
      for (INT ll = A->ia[kk]; ll < A->ia[kk + 1]; ll++) {
        INT m = A->ja[ll];
        if (cf[m] == RS_C_PT && is_strong[m] == i)
          denom += A->val[ll];
      }

      if (fabs(denom) < 1e-15) {
        diag_hat += a_ik;               /* no common C-point -> diagonal */
      } else {
        /* distribute a_ik to each l in C_i proportional to a_{kk,l} */
        for (INT ll = A->ia[kk]; ll < A->ia[kk + 1]; ll++) {
          INT m = A->ja[ll];
          if (cf[m] == RS_C_PT && is_strong[m] == i)
            P->val[pos[m]] += a_ik * A->val[ll] / denom;
        }
      }
    }

    /* Scale: w_{ij} = numerator / (-diag_hat) */
    if (fabs(diag_hat) > 1e-15) {
      for (INT pp = P->ia[i]; pp < P->ia[i + 1]; pp++)
        P->val[pp] /= -diag_hat;
    }

    /* Clean up pos[] (only reset entries we set for this row) */
    for (INT k = S->ia[i]; k < S->ia[i + 1]; k++) {
      INT j = S->ja[k];
      if (cf[j] == RS_C_PT) pos[j] = -1;
    }
  }

  free(is_strong);
  free(pos);
}

/* ----------------------------------------------------------------------
 *  3c. Standard (Direct) Interpolation — convenience wrapper
 *
 *  Builds both the sparsity pattern and the numerical values of P.
 * ---------------------------------------------------------------------- */
void rs_standard_interpolation(const dcsr_mat *A,
                               const dcsr_mat *S,
                               const INT      *cf,
                               dcsr_mat       *P)
{
  rs_standard_interpolation_sparsity(S, cf, P);
  rs_standard_interpolation_values(A, S, cf, P);
}
