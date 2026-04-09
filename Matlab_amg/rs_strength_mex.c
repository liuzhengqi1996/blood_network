/*
 * rs_strength_mex.c — Strength of connection, Ruge-Stuben criterion (MEX)
 *
 * Authors: Z. Liu and Ludmil Zikatanov
 *
 * Usage:
 *   [S, lam_max] = rs_strength_mex(A, theta)
 *
 * Inputs:
 *   A     : sparse n x n symmetric SPD matrix
 *   theta : strength threshold (scalar, e.g. 0.25)
 *
 * Outputs:
 *   S       : sparse n x n matrix containing only strong off-diagonal entries
 *             with their original values from A
 *   lam_max : scalar, max column count of A (used for bucket allocation)
 *
 * For each row i, off-diagonal j is "strong" if
 *     max(0, -a_{ij})  >=  theta * max_{k != i} max(0, -a_{ik})
 *
 * Because A is symmetric, MATLAB's CSC storage doubles as CSR:
 *   column i = row i.  The output S is built directly in CSC by
 *   iterating columns of A and filtering strong entries.
 */
#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize  n, i;
    mwIndex k;

    /* ---- check arguments ---- */
    if (nrhs != 2)
        mexErrMsgIdAndTxt("rs_str:nrhs",
            "Two inputs required: A, theta");
    if (!mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("rs_str:notSparse", "A must be sparse.");

    /* ---- unpack A (symmetric: CSC = CSR) ---- */
    n = mxGetM(prhs[0]);
    mwIndex *a_jc  = mxGetJc(prhs[0]);
    mwIndex *a_ir  = mxGetIr(prhs[0]);
    double  *a_val = mxGetPr(prhs[0]);
    double   theta = mxGetScalar(prhs[1]);

    /* ---- pass 1: row-wise max of max(0, -a_{ij}) for j != i ---- */
    double *row_max = (double *)mxCalloc(n, sizeof(double));

    int lambda_i;
    int max_lam = 0;

    for (i = 0; i < n; i++) {
        double mx = 0.0;
        for (k = a_jc[i]; k < a_jc[i + 1]; k++) {
            if (a_ir[k] == (mwIndex)i) continue;       /* skip diagonal */
            double v = -a_val[k];
            if (v < 0.0) v = 0.0;                      /* clamp */
            if (v > mx) mx = v;
        }
        row_max[i] = mx;
        lambda_i = (int)(a_jc[i + 1] - a_jc[i]);
        if (lambda_i > max_lam) max_lam = lambda_i;
    }

    if (max_lam < 1) max_lam = 1;




    /* ---- pass 2: count strong entries per column ---- */
    mwIndex *col_cnt = (mwIndex *)mxCalloc(n, sizeof(mwIndex));

    for (i = 0; i < n; i++) {                          /* column i of A */
        for (k = a_jc[i]; k < a_jc[i + 1]; k++) {
            mwIndex row = a_ir[k];
            if (row == (mwIndex)i) continue;            /* skip diagonal */
            double v = -a_val[k];
            if (v < 0.0) v = 0.0;
            if (v >= theta * row_max[row])
                col_cnt[i]++;
        }
    }

    /* ---- build column pointers ---- */
    mwSize s_nnz = 0;
    for (i = 0; i < n; i++) s_nnz += col_cnt[i];

    plhs[0] = mxCreateSparse(n, n, (s_nnz > 0 ? s_nnz : 1), mxREAL);
    mwIndex *s_jc  = mxGetJc(plhs[0]);
    mwIndex *s_ir  = mxGetIr(plhs[0]);
    double  *s_val = mxGetPr(plhs[0]);
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        *(int *)mxGetData(plhs[1]) = max_lam;
    }

    s_jc[0] = 0;
    for (i = 0; i < n; i++)
        s_jc[i + 1] = s_jc[i] + col_cnt[i];

    /* ---- pass 3: fill S in CSC ---- */
    memset(col_cnt, 0, n * sizeof(mwIndex));

    for (i = 0; i < n; i++) {                          /* column i of A */
        for (k = a_jc[i]; k < a_jc[i + 1]; k++) {
            mwIndex row = a_ir[k];
            if (row == (mwIndex)i) continue;
            double v = -a_val[k];
            if (v < 0.0) v = 0.0;
            if (v >= theta * row_max[row]) {
                mwIndex pos = s_jc[i] + col_cnt[i];
                s_ir[pos]  = row;
                s_val[pos] = a_val[k];                 /* original value */
                col_cnt[i]++;
            }
        }
    }

    /* ---- cleanup ---- */
    mxFree(row_max);
    mxFree(col_cnt);
}
