/*
 * rs_interpolation_mex.c — Standard (direct) RS interpolation (MEX)
 *
 * Authors: Z. Liu and Ludmil Zikatanov
 *
 * Usage:
 *   P = rs_interpolation_mex(A, S, cf)
 *
 * Inputs:
 *   A  : sparse n x n symmetric SPD matrix
 *   S  : sparse n x n strength-of-connection matrix
 *   cf : n x 1 double vector (1 = C-point, -1 = F-point)
 *
 * Output:
 *   P  : sparse n x nc interpolation matrix
 *
 * Because A is symmetric, MATLAB's CSC = CSR for row access.
 * S is not symmetric; CSR is built by transposing the CSC data.
 *
 * The algorithm follows the standard RS direct interpolation:
 *   C-points: injection  P(i, cmap(i)) = 1
 *   F-points: weighted interpolation from strong C-neighbors,
 *             with strong-F-neighbor redistribution.
 */
#include "mex.h"
#include <string.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize  n, i;
    mwIndex k;

    /* ---- check arguments ---- */
    if (nrhs != 3)
        mexErrMsgIdAndTxt("rs_interp:nrhs",
            "Three inputs required: A, S, cf");
    if (!mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("rs_interp:A", "A must be sparse.");
    if (!mxIsSparse(prhs[1]))
        mexErrMsgIdAndTxt("rs_interp:S", "S must be sparse.");

    /* ---- unpack A (symmetric: CSC = CSR) ---- */
    n = mxGetM(prhs[0]);
    mwIndex *a_jc  = mxGetJc(prhs[0]);
    mwIndex *a_ir  = mxGetIr(prhs[0]);
    double  *a_val = mxGetPr(prhs[0]);

    /* ---- unpack S (CSC) and build CSR ---- */
    mwIndex *s_jc  = mxGetJc(prhs[1]);
    mwIndex *s_ir  = mxGetIr(prhs[1]);
    mwSize   s_nnz = s_jc[n];

    int *s_ia = (int *)mxCalloc(n + 1, sizeof(int));
    int *s_ja = (int *)mxMalloc((s_nnz > 0 ? s_nnz : 1) * sizeof(int));

    for (k = 0; k < (mwIndex)s_nnz; k++)
        s_ia[s_ir[k] + 1]++;
    for (i = 0; i < n; i++)
        s_ia[i + 1] += s_ia[i];

    {
        int *cnt = (int *)mxCalloc(n, sizeof(int));
        for (i = 0; i < n; i++) {
            for (k = s_jc[i]; k < s_jc[i + 1]; k++) {
                int row = (int)s_ir[k];
                s_ja[s_ia[row] + cnt[row]] = (int)i;
                cnt[row]++;
            }
        }
        mxFree(cnt);
    }

    /* ---- unpack cf ---- */
    const double *cf_dbl = mxGetPr(prhs[2]);
    int *cf = (int *)mxMalloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        cf[i] = (int)cf_dbl[i];

    /* ---- coarse-grid index map (0-based) ---- */
    int *cmap = (int *)mxMalloc(n * sizeof(int));
    int  nc = 0;
    for (i = 0; i < n; i++) {
        if (cf[i] == 1) cmap[i] = nc++;
        else            cmap[i] = -1;
    }

    /* ---- workspace (reused per row) ---- */
    int    *is_strong = (int *)mxMalloc(n * sizeof(int));
    int    *pos       = (int *)mxMalloc(n * sizeof(int));
    int    *C_i       = (int *)mxMalloc(n * sizeof(int));
    double *w         = (double *)mxMalloc(n * sizeof(double));
    memset(is_strong, -1, n * sizeof(int));
    memset(pos,       -1, n * sizeof(int));

    /* ---- COO triplets for P ---- */
    mwSize max_nnz_P = s_nnz + (mwSize)nc + 1;
    int    *p_row = (int *)   mxMalloc(max_nnz_P * sizeof(int));
    int    *p_col = (int *)   mxMalloc(max_nnz_P * sizeof(int));
    double *p_val = (double *)mxMalloc(max_nnz_P * sizeof(double));
    mwSize  p_nnz = 0;

    /* ---- main loop: build P row by row ---- */
    for (i = 0; i < n; i++) {

        /* ---- C-point: injection ---- */
        if (cf[i] == 1) {
            p_row[p_nnz] = (int)i;
            p_col[p_nnz] = cmap[i];
            p_val[p_nnz] = 1.0;
            p_nnz++;
            continue;
        }

        /* ---- F-point: standard interpolation ---- */

        /* mark strong connections of row i and collect C_i */
        int n_ci = 0;
        for (int kk = s_ia[i]; kk < s_ia[i + 1]; kk++) {
            int j = s_ja[kk];
            is_strong[j] = (int)i;          /* mark j as strong for row i */
            if (cf[j] == 1) {
                pos[j]    = n_ci;
                C_i[n_ci] = j;
                w[n_ci]   = 0.0;
                n_ci++;
            }
        }

        if (n_ci == 0) {
            /* no strong C-neighbors — clear marks and skip */
            for (int kk = s_ia[i]; kk < s_ia[i + 1]; kk++)
                is_strong[s_ja[kk]] = -1;
            continue;
        }

        /* walk row i of A  (= column i of CSC, since A is symmetric) */
        double diag_val = 0.0, weak_sum = 0.0;

        for (k = a_jc[i]; k < a_jc[i + 1]; k++) {
            int j   = (int)a_ir[k];
            double a = a_val[k];

            if (j == (int)i) {
                diag_val = a;                           /* diagonal */
            } else if (is_strong[j] == (int)i && cf[j] == 1) {
                w[pos[j]] += a;                         /* strong C-neighbor */
            } else if (is_strong[j] == (int)i && cf[j] == -1) {
                /* strong F-neighbor — handled below */
            } else {
                weak_sum += a;                          /* weak → lump */
            }
        }
        double diag_hat = diag_val + weak_sum;

        /* distribute strong-F-neighbor contributions */
        for (k = a_jc[i]; k < a_jc[i + 1]; k++) {
            int kk = (int)a_ir[k];
            if (kk == (int)i || is_strong[kk] != (int)i || cf[kk] != -1)
                continue;

            double a_ik = a_val[k];

            /* denominator: sum_{l in C_i} a_{kk,l}  (row kk of A) */
            double denom = 0.0;
            for (mwIndex ll = a_jc[kk]; ll < a_jc[kk + 1]; ll++) {
                int m = (int)a_ir[ll];
                if (cf[m] == 1 && is_strong[m] == (int)i)
                    denom += a_val[ll];
            }

            if (fabs(denom) < 1e-15) {
                diag_hat += a_ik;                       /* no common C → lump */
            } else {
                for (mwIndex ll = a_jc[kk]; ll < a_jc[kk + 1]; ll++) {
                    int m = (int)a_ir[ll];
                    if (cf[m] == 1 && is_strong[m] == (int)i)
                        w[pos[m]] += a_ik * a_val[ll] / denom;
                }
            }
        }

        /* scale: w_{ij} = numerator / (-diag_hat) */
        if (fabs(diag_hat) > 1e-15) {
            for (int idx = 0; idx < n_ci; idx++)
                w[idx] /= -diag_hat;
        }

        /* store triplets */
        for (int idx = 0; idx < n_ci; idx++) {
            p_row[p_nnz] = (int)i;
            p_col[p_nnz] = cmap[C_i[idx]];
            p_val[p_nnz] = w[idx];
            p_nnz++;
        }

        /* clean up workspace for this row */
        for (int kk = s_ia[i]; kk < s_ia[i + 1]; kk++) {
            int j = s_ja[kk];
            is_strong[j] = -1;
            if (cf[j] == 1) pos[j] = -1;
        }
    }

    /* ---- convert COO → CSC for mxCreateSparse ---- */
    plhs[0] = mxCreateSparse(n, (mwSize)nc, (p_nnz > 0 ? p_nnz : 1), mxREAL);
    mwIndex *pjc = mxGetJc(plhs[0]);
    mwIndex *pir = mxGetIr(plhs[0]);
    double  *ppr = mxGetPr(plhs[0]);

    /* count entries per column of P */
    for (k = 0; k < (mwIndex)p_nnz; k++)
        pjc[p_col[k] + 1]++;
    for (int j = 0; j < nc; j++)
        pjc[j + 1] += pjc[j];

    /* fill CSC (rows arrive in ascending order per column because
     * the main loop iterates i = 0 .. n-1 sequentially) */
    {
        int *cnt = (int *)mxCalloc(nc, sizeof(int));
        for (k = 0; k < (mwIndex)p_nnz; k++) {
            int col     = p_col[k];
            mwIndex dst = pjc[col] + cnt[col];
            pir[dst]    = (mwIndex)p_row[k];
            ppr[dst]    = p_val[k];
            cnt[col]++;
        }
        mxFree(cnt);
    }

    /* ---- cleanup ---- */
    mxFree(cf);   mxFree(cmap);
    mxFree(s_ia); mxFree(s_ja);
    mxFree(is_strong); mxFree(pos);
    mxFree(C_i);  mxFree(w);
    mxFree(p_row); mxFree(p_col); mxFree(p_val);
}
