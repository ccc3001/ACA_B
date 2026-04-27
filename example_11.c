#include "basic.h"
#include "cluster.h"
#include "interpolation.h"
#include "supermatrix.h"
#include "blas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*****************************************************
main function
******************************************************/
int
main(int argc, char **argv) {
  int n, d, i, j, k, nmin, rank, nd, seed, col;
  double eta, eps, h, total_t, exact, error, oo_error;
  double *nodes_x, *u, *v, *xn, *yn, *e_col, *error_col;
  clock_t start_t, end_t;
  size_t size_supermatrix, size_fullmatrix;
  pclustertree ct;
  psupermatrix s;
  FILE *fp;

  if (argv[1] == NULL) {
    printf("****************************************************************\n");
    printf("Execute ./example_11 n d rank nmin eta eps\nn: number of points ");
    printf("per spatial dimension\nd: spatial dimension\nrank: ");
    printf("rank for ACA\nnmin: stopping criterion for cluster bisection");
    printf("\neta: admissibility condition\n");
    printf("eps: accuracy for recompression of supermatrix\n");
    printf("****************************************************************\n");
    exit(EXIT_FAILURE);
  }
  else if (argc < 7) {
    printf("******************************************************\n");
    printf("Too few arguments.\n");
    printf("******************************************************\n");
    exit(EXIT_FAILURE);
  }
  else {
    n = atoi(argv[1]);  
    d = atoi(argv[2]);
    rank = atoi(argv[3]);
    nmin = atoi(argv[4]);
    eta = atof(argv[5]);
    eps = atof(argv[6]);
    printf("******************************************************\n");
    printf("Used values:\nn = %d\nd = %d\nrank = %d\n", n, d, rank);
    printf("nmin = %d\neta = %3.1f\neps = %g\n", nmin, eta, eps);
    printf("******************************************************\n");
    if (nmin < 2 * rank - 1) {
      printf("******************************************************\n");
      printf("nmin >= 2 * rank - 1 isn't fulfilled.\n");
      printf("So the rank is too large for the compressed QR decomposition.\n");
      printf("******************************************************\n");
      exit(EXIT_FAILURE);
    }
  }
  
  nd = (int) pow((double) n, (double) d);
  printf("Create %d^%d = %d points.\n", n, d, nd);
  seed = 4;
  nodes_x = allocate_doubles(d * nd);
  /* Create n^d random points in [0, 1]. */
  /* fill_random_entries(d * nd, nodes_x, seed); */
  /* Create n^d equidistant points in (-1, 1)^d. */
  h = 2.0 / (n + 1);
  if (d == 1) { 
    for (i = 0; i < n; i++) {
      nodes_x[i] = -1 + (i + 1) * h;
    }
  }
  else if (d == 2) { 
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        nodes_x[i*n + j] = -1 + (i+1) * h;
        nodes_x[n*n + i*n + j]= -1 + (j+1) * h;
      }
    }
  }
  else if (d == 3) { 
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) { 
        for (k = 0; k < n; k++) {
          nodes_x[i*n*n + j*n + k] = -1 + (i+1) * h;
          nodes_x[n*n*n + i*n*n + j*n + k] = -1 + (j+1) * h;
          nodes_x[2*n*n*n + i*n*n + j*n + k] = -1 + (k+1) * h;
          /* printf("node %d: (%f, %f, %f)\n", i*n*n + j*n + k,
                 -1 + (i+1)*h, -1 + (j+1)*h, -1 + (k+1)*h); */
        }
      }
    }
  }
  else {
    printf("Point creation not implemented for d = %d.\n", d);
    exit(EXIT_FAILURE);
  }   
  
  start_t = clock();
  ct = create_geometric_clustertree(nd, d, nodes_x, nmin);
  reorder_nodes(nd, d, nodes_x, ct->permute);
  s = build_supermatrix_from_cluster(ct->root, ct->root, eta);
  end_t = clock();
  total_t = (double) (end_t - start_t) / (double) CLOCKS_PER_SEC;
  printf("\nTotal time taken by CPU for clustertree and H-matrix creation: %.6g sec\n\n", total_t);

  start_t = clock();
  fill_kernelaca_supermatrix(s, d, rank, nodes_x, test_function_log, 0, 0, nd); 
  /* fill_kernelaca_supermatrix(s, d, rank, nodes_x, test_function_gaussian, 0, 0, nd); */
  /* fill_kernelaca_delta_supermatrix(s, d, rank, eps, nodes_x, test_function_log, 0, 0, nd); */
  /* fill_kernelaca_delta_supermatrix(s, d, rank, eps, nodes_x, test_function_gaussian, 0, 0, nd); */
  end_t = clock();
  total_t = (double) (end_t - start_t) / (double) CLOCKS_PER_SEC;
  printf("Total time taken by CPU to fill H-matrix: %.6g sec\n\n", total_t);

  size_supermatrix = getsize_supermatrix(s);
  printf("Size of storage in megabytes to save this %d x %d H-matrix: %.4g\n",
         nd, nd, convertb2_mb(size_supermatrix));
  size_fullmatrix = getsize_fullmatrix(nd, nd);
  printf("Size of storage in megabytes to save a %d x %d fullmatrix: %.4g\n",
         nd, nd, convertb2_mb(size_fullmatrix));
    
  if (nd < 50) {
    printf("\n");
    output_clustertree(ct); 
    write_supermatrix(s, 0, 0);
  }
  if (nd < 1000) {
    output_supermatrix(s, "s.ps");
  }

  fp = fopen("storage.txt", "a");
  if (fp == NULL) {
    printf("Opening file failed.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "%d %.4g\n", nd, convertb2_mb(size_supermatrix));
  fclose(fp);

  col = nd / 2;
  printf("\nAccuracy of the single column %d:\n", col); 
  v = allocate_doubles(nd);  
  u = allocate_doubles(nd); 
  e_col = allocate_doubles(nd);  
  error_col = allocate_doubles(nd); 
  xn = allocate_doubles(d);  
  yn = allocate_doubles(d);  
  clear_entries(nd, e_col);
  e_col[col] = 1.0;

  eval_supermatrix(s, e_col, v);
  for (j = 0; j < d; j++) { 
    yn[j] = nodes_x[j*nd + col];
  }
  for (i = 0; i < nd; i++) {
    for (j = 0; j < d; j++) { 
      xn[j] = nodes_x[j*nd + i]; 
    }
    u[i] = test_function_log(d, xn, yn); 
    /* u[i] = test_function_gaussian(d, xn, yn); */
    /* printf("%6.4g %6.4g %6.4g\n", v[i], u[i], fabs(v[i] - u[i])); */
  }
  
  exact = sqrt( ddot_(&nd, u, eins_, u, eins_) ); /* 2-norm of exact matrix column */
  dcopy_(&nd, u, eins_, error_col, eins_); /* error_col = u */
  daxpy_(&nd, dmeins_, v, eins_, error_col, eins_); /* error_col = error_col - v */
  error = sqrt( ddot_(&nd, error_col, eins_, error_col, eins_) );
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2                         = %.4g\n", col, col, error);
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2 / nd                    = %.4g\n", col, col, error / nd);
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2 / || A_full(:, %d) ||_2 = %.4g\n", col, col, col, error / exact);
  oo_error = 0.0;
  for (j = 0; j < nd; j++) {
    oo_error = max_double(fabs(error_col[j]), oo_error);
  }
  printf("|| A_full(:, %d) - A_h(:, %d) ||_oo                        = %.4g\n", col, col, oo_error);
    
  fp = fopen("error.txt", "a");
  if (fp == NULL) {
    printf("Opening file failed.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "%d %.4g %.4g %.4g %.4g\n", nd, error, error / nd, error / exact, oo_error);
  fclose(fp);
  
  start_t = clock();  
  recompress_supermatrix(s, eps);
  end_t = clock();
  total_t = (double) (end_t - start_t) / (double) CLOCKS_PER_SEC;
  printf("\nTotal time taken by CPU for H-matrix recompression: %.6g sec\n\n", total_t);
  
  size_supermatrix = getsize_supermatrix(s);
  printf("Size of storage in megabytes to save this %d x %d recompressed H-matrix: %.4g\n",
         nd, nd, convertb2_mb(size_supermatrix));
    
  if (nd < 1000) {
    output_supermatrix(s, "s_recompressed.ps");
  }

  fp = fopen("storage_recompressed.txt", "a");
  if (fp == NULL) {
    printf("Opening file failed.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "%d %.4g\n", nd, convertb2_mb(size_supermatrix));
  fclose(fp);
  
  printf("\nAccuracy of the single column %d:\n", col); 
  eval_supermatrix(s, e_col, v);
  
  dcopy_(&nd, u, eins_, error_col, eins_); /* error_col = u */
  daxpy_(&nd, dmeins_, v, eins_, error_col, eins_); /* error_col = error_col - v */
  error = sqrt( ddot_(&nd, error_col, eins_, error_col, eins_) );
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2                         = %.4g\n", col, col, error);
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2 / nd                    = %.4g\n", col, col, error / nd);
  printf("|| A_full(:, %d) - A_h(:, %d) ||_2 / || A_full(:, %d) ||_2 = %.4g\n", col, col, col, error / exact);
  oo_error = 0.0;
  for (j = 0; j < nd; j++) {
    oo_error = max_double(fabs(error_col[j]), oo_error);
  }
  printf("|| A_full(:, %d) - A_h(:, %d) ||_oo                        = %.4g\n", col, col, oo_error);
    
  fp = fopen("error_recompressed.txt", "a");
  if (fp == NULL) {
    printf("Opening file failed.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "%d %.4g %.4g %.4g %.4g\n", nd, error, error / nd, error / exact, oo_error);
  fclose(fp);

  del_doubles(nodes_x);
  del_doubles(u);
  del_doubles(v);
  del_doubles(e_col);
  del_doubles(error_col);
  del_doubles(xn);
  del_doubles(yn);
  del_clustertree(ct);
  del_supermatrix(s);
  
  return 0;
}
