/*
 * cf_gauss_seidel_mex.c — CF-ordered Gauss-Seidel smoother (MEX)
 *
 * Authors: Z. Liu and Ludmil Zikatanov
 *
 * Usage:
 *   x = cf_gauss_seidel_mex(A, b, x0, cf_order, nu, direction)
 *
 * Inputs:
 *   A         : sparse n x n symmetric SPD matrix
 *   b         : n x 1 RHS vector
 *   x0        : n x 1 initial guess
 *   cf_order  : n x 1 permutation vector (1-based, MATLAB convention)
 *   nu        : number of GS sweeps (scalar)
 *   direction : 0 = forward, 1 = backward
 *
 * Output:
 *   x : n x 1 updated solution after nu GS sweeps
 *
 * Key insight: because A is symmetric, MATLAB's CSC storage is equivalent
 * to CSR — column i contains the same entries as row i.  So jc[] serves
 * as "row pointers" and ir[] as "column indices" for row-wise access.
 */
#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize   n, ii;
    mwIndex  k;
    int      nu, backward, sweep;
    mwIndex *jc, *ir;
    double  *av, *x, *diag;
    const double *b, *x0, *cf_dbl;

    /* ---- check arguments ---- */
    if (nrhs != 6)
        mexErrMsgIdAndTxt("cf_gs:nrhs",
            "Six inputs required: A, b, x, cf_order, nu, direction");
    if (!mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("cf_gs:notSparse", "A must be sparse.");

    /* ---- unpack inputs ---- */
    jc  = mxGetJc(prhs[0]);
    ir  = mxGetIr(prhs[0]);
    av  = mxGetPr(prhs[0]);
    n   = mxGetM(prhs[0]);

    b       = mxGetPr(prhs[1]);
    x0      = mxGetPr(prhs[2]);
    cf_dbl  = mxGetPr(prhs[3]);
    nu      = (int)mxGetScalar(prhs[4]);
    backward = (int)mxGetScalar(prhs[5]);

    /* ---- allocate output (copy of x0) ---- */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x = mxGetPr(plhs[0]);
    memcpy(x, x0, n * sizeof(double));

    /* ---- pre-extract diagonal ---- */
    diag = (double *)mxMalloc(n * sizeof(double));
    for (ii = 0; ii < n; ii++) {
        diag[ii] = 1.0;                  /* safe fallback */
        for (k = jc[ii]; k < jc[ii + 1]; k++) {
            if (ir[k] == (mwIndex)ii) {
                diag[ii] = av[k];
                break;
            }
        }
    }

    /* ---- Gauss-Seidel sweeps ---- */
    if (!backward) {
        /* forward: visit cf_order(1), cf_order(2), ... , cf_order(n) */
        for (sweep = 0; sweep < nu; sweep++) {
            for (ii = 0; ii < n; ii++) {
                mwSize i = (mwSize)(cf_dbl[ii] - 1);   /* 1-based -> 0-based */
                double s = b[i];
                for (k = jc[i]; k < jc[i + 1]; k++) {
                    mwIndex j = ir[k];
                    if (j != (mwIndex)i)
                        s -= av[k] * x[j];
                }
                x[i] = s / diag[i];
            }
        }
    } else {
        /* backward: visit cf_order(n), cf_order(n-1), ... , cf_order(1) */
        for (sweep = 0; sweep < nu; sweep++) {
            for (ii = n; ii > 0; ii--) {
                mwSize i = (mwSize)(cf_dbl[ii - 1] - 1);
                double s = b[i];
                for (k = jc[i]; k < jc[i + 1]; k++) {
                    mwIndex j = ir[k];
                    if (j != (mwIndex)i)
                        s -= av[k] * x[j];
                }
                x[i] = s / diag[i];
            }
        }
    }

    mxFree(diag);
}
