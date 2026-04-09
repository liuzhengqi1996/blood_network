/*
 * rs_coarsening_mex.c — Classical Ruge-Stuben C/F splitting (MEX)
 *
 * Authors: Z. Liu and Ludmil Zikatanov
 *
 * Usage:
 *   [cf, nc] = rs_coarsening_mex(S, lam_max)
 *
 * Inputs:
 *   S       : sparse n x n strength-of-connection matrix (from rs_strength)
 *   lam_max : scalar, max column count of A (from rs_strength)
 *
 * Outputs:
 *   cf : n x 1 double vector (1 = C-point, -1 = F-point)
 *   nc : scalar, number of C-points
 */
#include "mex.h"
#include <string.h>

#define C_PT       1
#define F_PT     (-1)
#define UNDECIDED  0

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize   n, i;
    mwIndex  k;

    /* ---- check arguments ---- */
    if (nrhs != 2)
        mexErrMsgIdAndTxt("rs_coarsen:nrhs", "Two inputs required: S, lam_max");
    if (!mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("rs_coarsen:notSparse", "S must be sparse.");

    /* lam_max must be an int32 scalar */
    if (!mxIsInt32(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("rs_coarsen:badLamMax",
            "lam_max must be a scalar int32.");
    int lam_max = *(int *)mxGetData(prhs[1]);

    /* ---- unpack S (CSC) ---- */
    n = mxGetM(prhs[0]);
    mwIndex *s_jc = mxGetJc(prhs[0]);   /* column pointers          */
    mwIndex *s_ir = mxGetIr(prhs[0]);   /* row indices               */
    mwSize   s_nnz = s_jc[n];

   

    /* ---- build CSR for S ---- */
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

    /* ---- initialize lambda and cf ---- */
    int *lambda = (int *)mxCalloc(n, sizeof(int));
    int *cf     = (int *)mxCalloc(n, sizeof(int));   /* UNDECIDED = 0 */

    for (i = 0; i < n; i++)
        lambda[i] = (int)(s_jc[i + 1] - s_jc[i]);   /* |S^T row i| */

    int num_buckets = 2 * lam_max + 1;
    if (num_buckets < 1) num_buckets = 1;

    int *bucket_head = (int *)mxMalloc(num_buckets * sizeof(int));
    int *nxt = (int *)mxMalloc(n * sizeof(int));
    int *prv = (int *)mxMalloc(n * sizeof(int));

    for (int b = 0; b < num_buckets; b++)
        bucket_head[b] = -1;

    /* insert all points (reverse order → lowest index first in bucket) */
    for (int ii = (int)n - 1; ii >= 0; ii--) {
        int b = lambda[ii];
        nxt[ii] = bucket_head[b];
        prv[ii] = -1;
        if (bucket_head[b] >= 0)
            prv[bucket_head[b]] = ii;
        bucket_head[b] = ii;
    }

    int top_bucket = lam_max;

    /* ---- first pass: greedy selection ---- */
    for (;;) {
        /* find highest non-empty bucket */
        while (top_bucket >= 0 && bucket_head[top_bucket] < 0)
            top_bucket--;
        if (top_bucket < 0) break;

        int best = bucket_head[top_bucket];

        /* remove best from its bucket */
        {
            int b = lambda[best];
            if (prv[best] >= 0) nxt[prv[best]] = nxt[best];
            else                bucket_head[b]  = nxt[best];
            if (nxt[best] >= 0) prv[nxt[best]]  = prv[best];
        }
        cf[best] = C_PT;

        /* for each undecided j influenced by best (S^T row best) */
        for (k = s_jc[best]; k < s_jc[best + 1]; k++) {
            int j = (int)s_ir[k];
            if (cf[j] != UNDECIDED) continue;

            /* remove j from its bucket */
            {
                int b = lambda[j];
                if (prv[j] >= 0) nxt[prv[j]] = nxt[j];
                else             bucket_head[b] = nxt[j];
                if (nxt[j] >= 0) prv[nxt[j]]   = prv[j];
            }
            cf[j] = F_PT;

            /* j became F ⇒ increment lambda for undecided strong deps of j */
            for (int ll = s_ia[j]; ll < s_ia[j + 1]; ll++) {
                int kk = s_ja[ll];
                if (cf[kk] != UNDECIDED) continue;

                /* remove kk from its current bucket */
                {
                    int b = lambda[kk];
                    if (prv[kk] >= 0) nxt[prv[kk]] = nxt[kk];
                    else              bucket_head[b] = nxt[kk];
                    if (nxt[kk] >= 0) prv[nxt[kk]]  = prv[kk];
                }

                lambda[kk]++;
                int new_b = lambda[kk];

                /* insert kk into new bucket */
                nxt[kk] = bucket_head[new_b];
                prv[kk] = -1;
                if (bucket_head[new_b] >= 0)
                    prv[bucket_head[new_b]] = kk;
                bucket_head[new_b] = kk;

                if (new_b > top_bucket)
                    top_bucket = new_b;
            }
        }
    }

    /* ---- second pass: F-points with no strong C-neighbor → C ---- */
    for (i = 0; i < n; i++) {
        if (cf[i] != F_PT) continue;
        int has_c = 0;
        for (int kk = s_ia[i]; kk < s_ia[i + 1]; kk++) {
            if (cf[s_ja[kk]] == C_PT) { has_c = 1; break; }
        }
        if (!has_c) cf[i] = C_PT;
    }

    /* ---- count coarse points ---- */
    int nc = 0;
    for (i = 0; i < n; i++)
        if (cf[i] == C_PT) nc++;

    /* ---- outputs ---- */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *cf_out = mxGetPr(plhs[0]);
    for (i = 0; i < n; i++)
        cf_out[i] = (double)cf[i];

    if (nlhs > 1)
        plhs[1] = mxCreateDoubleScalar((double)nc);

    /* ---- cleanup ---- */
    mxFree(s_ia);  mxFree(s_ja);
    mxFree(lambda); mxFree(cf);
    mxFree(bucket_head); mxFree(nxt); mxFree(prv);
}
