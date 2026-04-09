/*
 * Vector Operations
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Translated from Fortran srclib.f
 *
 * Uses 0-based indexing (C style).
 */

#include "agg_lib.h"

/*
 * nullv - Zero a REAL precision vector
 */
void nullv(REAL *u, INT n)
{
    for (INT i = 0; i < n; i++) {
        u[i] = 0.0;
    }
}

/*
 * copyv - Copy a REAL precision vector: v = u
 */
void copyv(const REAL *u, REAL *v, INT n)
{
    for (INT i = 0; i < n; i++) {
        v[i] = u[i];
    }
}

/*
 * scpro - Inner product: scpr = u^T * v
 */
void scpro(const REAL *u, const REAL *v, REAL *scpr, INT n)
{
    *scpr = 0.0;
    for (INT i = 0; i < n; i++) {
        *scpr += u[i] * v[i];
    }
}

/*
 * l2norm - L2 norm: l2nr = ||u||_2
 */
void l2norm(const REAL *u, REAL *l2nr, INT n)
{
    scpro(u, u, l2nr, n);
    *l2nr = sqrt(*l2nr);
}

/*
 * uuplmv - u = u + smult * v
 */
void uuplmv(REAL *u, const REAL *v, INT n, REAL smult)
{
    for (INT i = 0; i < n; i++) {
        u[i] = u[i] + v[i] * smult;
    }
}

/*
 * umuplv - u = smult * u + v
 */
void umuplv(REAL *u, const REAL *v, INT n, REAL smult)
{
    for (INT i = 0; i < n; i++) {
        u[i] = smult * u[i] + v[i];
    }
}

/*
 * randv - Generate random vector with entries in [0, 1]
 */
void randv(REAL *u, INT n)
{
    for (INT i = 0; i < n; i++) {
        u[i] = (REAL)rand() / (REAL)RAND_MAX;
    }
}
