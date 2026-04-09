/*
 * build_hierarchy.c - Build multilevel AMG hierarchy
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Classical Ruge-Stuben coarsening using rs_strength, rs_coarsening,
 * rs_standard_interpolation from rs_classical.c.
 */
#include "amg.h"

/* ======================================================================
 *  Classical Ruge-Stuben hierarchy
 * ====================================================================== */
void amg_build_hierarchy(amg_data* amg, const dcsr_mat* A) {
  REAL threshold = amg->threshold;
  INT min_size   = amg->min_size;
  INT max_levels = amg->max_levels;

  if (max_levels > AMG_MAX_LEVELS) max_levels = AMG_MAX_LEVELS;

  /* Copy A to level 0 */
  amg_level* lev = &amg->levels[0];
  dcsr_mat_alloc(&lev->A, A->row, A->col, A->nnz);
  dcsr_cp(A, &lev->A);
  lev->n = A->row;
  lev->mis = NULL;
  lev->isolated = NULL;
  lev->cf_order = NULL;

  fprintf(stderr,
    "Building RS hierarchy (threshold=%.2f, min_size=%d)\n",
    threshold, min_size);
  fprintf(stderr, "  Level 1: n = %d, nnz = %d\n", A->row, A->nnz);

  INT nlev = 1;
  const dcsr_mat* Ac = &amg->levels[0].A;

  while (Ac->row >= min_size && nlev < max_levels) {
    INT n_fine = Ac->row;
    amg_level* cur = &amg->levels[nlev - 1];

    /* Step 1: Strength of connection (stored in A_filtered) */
    INT max_lam = 0;
    rs_strength(Ac, threshold, &cur->A_filtered, &max_lam);

    /* Step 2: Classical RS coarsening */
    INT *cf = (INT *)calloc(n_fine, sizeof(INT));
    INT n_coarse = 0;
    rs_coarsening(&cur->A_filtered, max_lam, cf, &n_coarse);
    cur->n_coarse = n_coarse;

    if (n_coarse == 0 || n_coarse >= n_fine) {
      fprintf(stderr, "  Coarsening stalled at level %d (n_coarse=%d)\n",
              nlev, n_coarse);
      free(cf);
      break;
    }

    /* Populate mis[] and isolated[] for compatibility with V-cycle */
    cur->mis = (INT *)calloc(n_fine, sizeof(INT));
    cur->isolated = (INT *)calloc(n_fine, sizeof(INT));
    for (INT i = 0; i < n_fine; i++)
      cur->mis[i] = (cf[i] == 1) ? 1 : 0;   /* C-point = 1 */

    /* Build CF ordering: C-points first, then F-points */
    cur->cf_order = (INT *)malloc(n_fine * sizeof(INT));
    INT idx = 0;
    for (INT i = 0; i < n_fine; i++)
      if (cf[i] == 1) cur->cf_order[idx++] = i;
    for (INT i = 0; i < n_fine; i++)
      if (cf[i] != 1) cur->cf_order[idx++] = i;

    /* Step 3: Standard interpolation */
    rs_standard_interpolation(&cur->A, &cur->A_filtered, cf, &cur->P);

    /* Step 3b: Smooth P (near-harmonic iteration) */
    /* smooth_P(&cur->A_filtered, &cur->P, cur->mis, cur->isolated,
             amg->maxit_prol); */

    free(cf);

    /* Step 4: Galerkin coarse-grid operator Ac = P^T A P */
    dcsr_mat Pt;
    dcsr_mat_alloc(&Pt, cur->P.col, cur->P.row, cur->P.nnz);
    dcsr_transz(&cur->P, NULL, &Pt);

    dcsr_mat Ac_new;
    dcsr_rap(&Pt, &cur->A, &cur->P, &Ac_new);
    dcsr_mat_free(&Pt);

    /* Symmetrize: Ac = (Ac + Ac') / 2 */
    {
      dcsr_mat Ac_t;
      dcsr_mat_alloc(&Ac_t, Ac_new.col, Ac_new.row, Ac_new.nnz);
      dcsr_transz(&Ac_new, NULL, &Ac_t);

      dcsr_mat Ac_sum;
      dcsr_mat_add(&Ac_new, &Ac_t, &Ac_sum);
      for (INT k = 0; k < Ac_sum.nnz; k++) Ac_sum.val[k] *= 0.5;

      dcsr_mat_free(&Ac_new);
      dcsr_mat_free(&Ac_t);
      Ac_new = Ac_sum;
    }

    fprintf(stderr, "  Level %d: n = %d, nnz = %d (ratio %.2f)\n",
            nlev + 1, n_coarse, Ac_new.nnz, (REAL)n_fine / n_coarse);

    /* Store coarse level */
    amg_level* next = &amg->levels[nlev];
    next->A = Ac_new;
    next->n = n_coarse;
    next->mis = NULL;
    next->isolated = NULL;
    next->cf_order = NULL;

    nlev++;
    Ac = &amg->levels[nlev - 1].A;
  }

  amg->num_levels = nlev;

  /* Compute ichol on the coarsest level for the coarse solver */
  amg_level* coarsest = &amg->levels[nlev - 1];
  ichol_compute(&coarsest->A, &coarsest->L_ichol);
  fprintf(stderr, "RS hierarchy complete: %d levels (coarsest ichol nnz=%d)\n\n",
          nlev, coarsest->L_ichol.nnz);
}

/* ======================================================================
 *  Rebuild hierarchy values for a spectrally equivalent matrix
 *
 *  Reuses the existing CF splitting, P sparsity pattern, and strength
 *  graph structure from a previously built hierarchy.  Only recomputes:
 *    - A at each level (level 0 from A_new, coarser via RAP)
 *    - P values via rs_standard_interpolation_values
 *    - ichol on the coarsest level
 *
 *  Precondition: A_new must have the same sparsity pattern as the
 *  matrix used to build the original hierarchy.
 * ====================================================================== */
void amg_rebuild_values(amg_data* amg, const dcsr_mat* A_new) {
  INT nlev = amg->num_levels;

  /* Free old ichol on coarsest level */
  ichol_free(&amg->levels[nlev - 1].L_ichol);

  /* Replace level 0 matrix */
  amg_level* lev0 = &amg->levels[0];
  dcsr_mat_free(&lev0->A);
  dcsr_mat_alloc(&lev0->A, A_new->row, A_new->col, A_new->nnz);
  dcsr_cp(A_new, &lev0->A);

  fprintf(stderr,
    "Rebuilding AMG values (reusing hierarchy structure)\n");
  fprintf(stderr, "  Level 1: n = %d, nnz = %d\n", A_new->row, A_new->nnz);

  for (INT k = 0; k < nlev - 1; k++) {
    amg_level* cur = &amg->levels[k];
    INT n = cur->n;

    /* Reconstruct cf array from mis (mis[i]=1 for C, 0 for F) */
    INT *cf = (INT *)malloc(n * sizeof(INT));
    for (INT i = 0; i < n; i++)
      cf[i] = cur->mis[i] ? 1 : -1;  /* RS_C_PT=1, RS_F_PT=-1 */

    /* Recompute P values (sparsity pattern is preserved) */
    rs_standard_interpolation_values(&cur->A, &cur->A_filtered, cf, &cur->P);

    free(cf);

    /* Recompute coarse-grid operator: Ac = P^T A P */
    dcsr_mat Pt;
    dcsr_mat_alloc(&Pt, cur->P.col, cur->P.row, cur->P.nnz);
    dcsr_transz(&cur->P, NULL, &Pt);

    dcsr_mat Ac_new;
    dcsr_rap(&Pt, &cur->A, &cur->P, &Ac_new);
    dcsr_mat_free(&Pt);

    /* Symmetrize: Ac = (Ac + Ac') / 2 */
    {
      dcsr_mat Ac_t;
      dcsr_mat_alloc(&Ac_t, Ac_new.col, Ac_new.row, Ac_new.nnz);
      dcsr_transz(&Ac_new, NULL, &Ac_t);

      dcsr_mat Ac_sum;
      dcsr_mat_add(&Ac_new, &Ac_t, &Ac_sum);
      for (INT j = 0; j < Ac_sum.nnz; j++) Ac_sum.val[j] *= 0.5;

      dcsr_mat_free(&Ac_new);
      dcsr_mat_free(&Ac_t);
      Ac_new = Ac_sum;
    }

    fprintf(stderr, "  Level %d: n = %d, nnz = %d\n",
            k + 2, Ac_new.row, Ac_new.nnz);

    /* Replace coarse level matrix */
    amg_level* next = &amg->levels[k + 1];
    dcsr_mat_free(&next->A);
    next->A = Ac_new;
  }

  /* Recompute ichol on the coarsest level */
  amg_level* coarsest = &amg->levels[nlev - 1];
  ichol_compute(&coarsest->A, &coarsest->L_ichol);
  fprintf(stderr, "AMG values rebuilt: %d levels (coarsest ichol nnz=%d)\n\n",
          nlev, coarsest->L_ichol.nnz);
}

/* ====================================================================== */
void amg_free(amg_data* amg) {
  for (INT k = 0; k < amg->num_levels; k++) {
    amg_level* lev = &amg->levels[k];
    dcsr_mat_free(&lev->A);
    if (k < amg->num_levels - 1) {
      dcsr_mat_free(&lev->P);
      dcsr_mat_free(&lev->A_filtered);
    }
    if (lev->mis) { free(lev->mis); lev->mis = NULL; }
    if (lev->isolated) { free(lev->isolated); lev->isolated = NULL; }
    if (lev->cf_order) { free(lev->cf_order); lev->cf_order = NULL; }
    if (k == amg->num_levels - 1) ichol_free(&lev->L_ichol);
  }
  amg->num_levels = 0;
}
