/*
 * test_amg.c - Test AMG V-cycle and AMG+ichol preconditioners
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Usage: ./test_amg <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]
 *
 * One matrix:  ./test_amg matrixA.csr 0.125 2 25
 *   Build AMG from A, ichol from A, solve Ax=b.
 *
 * Two matrices: ./test_amg matrixA.csr matrixP.csr 0.125 2 25
 *   Build AMG+ichol from Ap (matrixP.csr), solve Ax=b.
 *
 * Matrix file format (1-indexed COO, from matrix_save.m):
 *   nr nc nnz
 *   i j value   (repeated nnz times)
 */
#include "amg.h"
#include <sys/resource.h>

/* Preconditioner wrappers for amg_pcg function pointer */

typedef struct {
  const amg_data* amg;
  const dcsr_mat* A;
  const ichol_factor* L;
} pcg_precond_data;

/* static void vcycle_precond_wrap(const void* data, const REAL* g, REAL* Bg) {
  const amg_data* amg = (const amg_data*)data;
  amg_vcycle_precond(amg, g, Bg);
} */

static void ichol_precond_wrap(const void* data, const REAL* g, REAL* Bg) {
  const pcg_precond_data* d = (const pcg_precond_data*)data;
  amg_ichol_precond(d->amg, d->A, d->L, g, Bg);
}

/* Load and symmetrize a matrix from a COO file */
static int load_matrix(const char* filename, dcsr_mat* A) {
  FILE* fp = fopen(filename, "r");
  if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return 1; }

  if (dcsr_mat_dcoo(fp, A) != 0) {
    fprintf(stderr, "Error reading matrix from %s\n", filename);
    fclose(fp);
    return 1;
  }
  fclose(fp);

  /* Symmetrize: A = (A + A') / 2 */
  dcsr_mat At;
  dcsr_mat_alloc(&At, A->col, A->row, A->nnz);
  dcsr_transz(A, NULL, &At);

  dcsr_mat Asum;
  dcsr_mat_add(A, &At, &Asum);
  for (INT k = 0; k < Asum.nnz; k++) Asum.val[k] *= 0.5;

  dcsr_mat_free(A);
  dcsr_mat_free(&At);
  *A = Asum;

  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]\n", argv[0]);
    return 1;
  }

  /* Parse arguments: detect whether arg2 is a filename or a number */
  const char* solve_file = argv[1];
  const char* amg_file   = NULL;
  int arg_offset = 2;  /* index of first numeric argument */

  if (argc >= 3) {
    /* If arg2 starts with a digit or '.', treat as threshold (one-matrix mode) */
    char c = argv[2][0];
    if (c != '.' && (c < '0' || c > '9')) {
      amg_file = argv[2];
      arg_offset = 3;
    }
  }

  REAL threshold = (argc > arg_offset)     ? atof(argv[arg_offset])     : 0.25;
  INT  nu        = (argc > arg_offset + 1) ? atoi(argv[arg_offset + 1]) : 2;
  INT  min_size  = (argc > arg_offset + 2) ? atoi(argv[arg_offset + 2]) : 4096;
  INT  coarsen   = (argc > arg_offset + 3) ? atoi(argv[arg_offset + 3]) : COARSEN_PMIS;

  int two_matrix = (amg_file != NULL);

  /* Load solve matrix A */
  fprintf(stderr, "Loading solve matrix from %s ...\n", solve_file);
  clock_t t0 = clock();

  dcsr_mat A;
  if (load_matrix(solve_file, &A) != 0) return 1;

  INT n = A.row;
  double t_load = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  n = %d, nnz = %d  (%.2f s)\n", n, A.nnz, t_load);

  /* Load AMG matrix Ap (or alias to A) */
  dcsr_mat Ap_storage;
  dcsr_mat* Ap;
  if (two_matrix) {
    fprintf(stderr, "Loading AMG matrix from %s ...\n", amg_file);
    t0 = clock();
    if (load_matrix(amg_file, &Ap_storage) != 0) return 1;
    double t_load2 = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "  n = %d, nnz = %d  (%.2f s)\n",
            Ap_storage.row, Ap_storage.nnz, t_load2);
    if (Ap_storage.row != n) {
      fprintf(stderr, "ERROR: matrix dimensions do not match (%d vs %d)\n",
              Ap_storage.row, n);
      return 1;
    }
    Ap = &Ap_storage;
  } else {
    Ap = &A;
  }

  /* RHS: try rhs.txt, fall back to random */
  REAL* b = (REAL*)malloc(n * sizeof(REAL));
  REAL* x = (REAL*)calloc(n, sizeof(REAL));
  {
    FILE* frhs = fopen("rhs.txt", "r");
    if (frhs) {
      INT nrhs;
      if (fscanf(frhs, "%d", &nrhs) == 1 && nrhs == n) {
        for (INT i = 0; i < n; i++) {
          if (fscanf(frhs, "%lf", &b[i]) != 1) {
            fprintf(stderr, "Error reading rhs.txt at entry %d, using random RHS\n", i);
            randv(b, n);
            break;
          }
        }
        fprintf(stderr, "RHS loaded from rhs.txt\n");
      } else {
        fprintf(stderr, "rhs.txt dimension mismatch (expected %d), using random RHS\n", n);
        randv(b, n);
      }
      fclose(frhs);
    } else {
      fprintf(stderr, "rhs.txt not found, using random RHS\n");
      randv(b, n);
    }
  }

  /* Build AMG hierarchy from Ap */
  amg_data amg;
  memset(&amg, 0, sizeof(amg));
  amg.threshold  = threshold;
  amg.min_size   = min_size;
  amg.max_levels = AMG_MAX_LEVELS;
  amg.maxit_prol = 0;
  amg.lump       = 0;
  amg.nu         = nu;
  amg.smoother     = SMOOTHER_CF_GS;
  amg.omega_jacobi = 0.5;
  amg.coarsening   = coarsen;

  if (two_matrix)
    fprintf(stderr, "\nBuilding AMG hierarchy from %s ...\n", amg_file);

  t0 = clock();
  amg_build_hierarchy(&amg, Ap);
  double t_hier = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "Hierarchy built in %.2f s\n", t_hier);

  /* Report nnz per F-row in P at each level */
  for (INT lev = 0; lev < amg.num_levels - 1; lev++) {
    const amg_level* L_lev = &amg.levels[lev];
    const dcsr_mat* P = &L_lev->P;
    const INT* mis = L_lev->mis;
    INT nf = 0, min_nnz = P->row, max_nnz = 0;
    long total_nnz = 0;
    for (INT i = 0; i < P->row; i++) {
      if (mis[i]) continue; /* skip C-points */
      nf++;
      INT rnnz = P->ia[i + 1] - P->ia[i];
      total_nnz += rnnz;
      if (rnnz < min_nnz) min_nnz = rnnz;
      if (rnnz > max_nnz) max_nnz = rnnz;
    }
    fprintf(stderr, "  Level %d P: %d F-rows, nnz/F-row: min=%d, max=%d, avg=%.1f\n",
            lev + 1, nf, min_nnz, max_nnz, nf > 0 ? (double)total_nnz / nf : 0.0);
  }

  /* Incomplete Cholesky: from Ap in two-matrix mode, from A otherwise */
  ichol_factor L;
  t0 = clock();
  if (two_matrix) {
    fprintf(stderr, "Computing ichol of AMG matrix %s ...\n", amg_file);
    ichol_compute(Ap, &L);
  } else {
    ichol_compute(&A, &L);
  }
  double t_ichol = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "ichol computed in %.2f s (nnz(L) = %d)\n", t_ichol, L.nnz);

  REAL tol = 1e-6;
  INT maxiter = 200;

  /* --- PCG with AMG+ichol --- */
  pcg_precond_data pdata;
  pdata.amg = &amg;
  pdata.A = two_matrix ? Ap : &A;
  pdata.L = &L;

  nullv(x, n);
  if (two_matrix)
    fprintf(stderr, "\nPCG solve of %s with AMG(%s)+ichol(%s):\n",
            solve_file, amg_file, amg_file);
  else
    fprintf(stderr, "\nPCG with AMG+ichol preconditioner:\n");
  t0 = clock();
  INT iter1 = amg_pcg(&A, b, x, ichol_precond_wrap, &pdata, tol, maxiter, 1);
  double t1 = (double)(clock() - t0) / CLOCKS_PER_SEC;

  /* Compute final relres */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  abyvam(b, A.ia, A.ja, A.val, x, n, r);
  REAL rnorm, bnorm;
  l2norm(r, &rnorm, n);
  l2norm(b, &bnorm, n);
  fprintf(stderr, "  iterations = %d, relres = %.2e, time = %.2f s\n",
          iter1, rnorm / bnorm, t1);

  /* --- PCG with AMG V-cycle only (commented out) --- */
  /*
  nullv(x, n);
  if (two_matrix)
    fprintf(stderr, "\nPCG solve of %s with AMG(%s) V-cycle only:\n",
            solve_file, amg_file);
  else
    fprintf(stderr, "\nPCG with AMG V-cycle preconditioner:\n");
  t0 = clock();
  INT iter2 = amg_pcg(&A, b, x, vcycle_precond_wrap, &amg, tol, maxiter, 1);
  double t2 = (double)(clock() - t0) / CLOCKS_PER_SEC;

  abyvam(b, A.ia, A.ja, A.val, x, n, r);
  l2norm(r, &rnorm, n);
  fprintf(stderr, "  iterations = %d, relres = %.2e, time = %.2f s\n",
          iter2, rnorm / bnorm, t2);
  */

  /* Summary */
  fprintf(stderr, "\n========================================\n");
  if (two_matrix)
    fprintf(stderr, "  AMG from %s, solve %s\n", amg_file, solve_file);
  fprintf(stderr, "  %-25s %6s %12s %12s\n", "Preconditioner", "Iter", "Setup(s)", "Solve(s)");
  fprintf(stderr, "  %-25s %6d %12.2f %12.2f\n", "AMG+ichol", iter1, t_hier + t_ichol, t1);
  /* fprintf(stderr, "  %-25s %6d %12.2f %12.2f\n", "AMG V-cycle only", iter2, t_hier, t2); */
  fprintf(stderr, "========================================\n");

  /* Report peak memory usage */
  {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    double mem_mb = usage.ru_maxrss / 1024.0;  /* ru_maxrss is in KB on Linux */
    if (mem_mb >= 1024.0)
      fprintf(stderr, "  Peak memory: %.2f GB\n", mem_mb / 1024.0);
    else
      fprintf(stderr, "  Peak memory: %.1f MB\n", mem_mb);
    fprintf(stderr, "========================================\n");
  }

  /* Cleanup */
  if (two_matrix) dcsr_mat_free(&Ap_storage);
  free(b); free(x); free(r);
  ichol_free(&L);
  amg_free(&amg);
  dcsr_mat_free(&A);

  return 0;
}
