#include "aca.h"
#include "basic.h"
#include "blas.h"

#include <math.h>
#include <assert.h>
#include <stdio.h>

double 
compute_entry_aca(pcrkmatrix r, int d, int k, int row, int col, 
                  const double *x, const double *y, 
                  double (*test_function)(int d, const double *x, const double *y)) {
  int i;
  double value;

  assert(r != NULL);
  assert(r->k >= k);  
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  assert(row >= 0);
  assert(col >= 0);

  /*************/
  /* exercise */
  /*************/
  value = test_function(d, x, y); 
  for (i = 0; i < k; i++) { 
    value -= r->a[i*r->rows + row] * r->b[i*r->cols + col];
  }

  return value;
}

prkmatrix
aca_rkmatrix(int d, int k, int rows, int cols, 
             const double *nodes_x, const double *nodes_y, 
             double (*test_function)(int d, const double *x, const double *y)) {
  int i, j, piv_row, piv_col;
  double max_col_value;
  double *x, *y;
  prkmatrix r;

  assert(nodes_x != NULL || d == 0 || rows == 0);
  assert(nodes_y != NULL || d == 0 || cols == 0);
  
  k = min(rows, k);
  k = min(cols, k);

  r = new_zero_rkmatrix(k, rows, cols);

  x = allocate_doubles(d * rows);
  y = allocate_doubles(d * cols);
  
  /* The rows d-dimensional nodes x_i are saved in nodes_x and the 
     cols d-dimensional nodes y_i are saved in nodes_y. 
     nodes_x and nodes_y begin with all first components, then all second
     components etc. We reorder the entries of nodes_x and nodes_y to begin
     with all d components of the first node, then all d components of the 
     second one, etc. and save them in x and y. */
  for (j = 0; j < d; j++) {
    for (i = 0; i < rows; i++) {
      x[d*i + j] = nodes_x[j*rows + i];
    }
  }
  for (j = 0; j < d; j++) {
    for (i = 0; i < cols; i++) {
      y[d*i + j] = nodes_y[j*cols + i];
    }
  }

  for (i = 0; i < k; i++) {
    /* Compute all entries of a pivot row. */
    piv_row = 0;
    max_col_value = 0.0;
    while (fabs(max_col_value) <= 1e-10 && piv_row < rows) {
      for (j = 0; j < cols; j++) {
        r->b[i*cols + j] = compute_entry_aca(r, d, i, piv_row, j, x + d*piv_row,
                                             y + d*j, test_function);
        if (fabs(max_col_value) < fabs(r->b[i*cols + j])) {
          max_col_value = r->b[i*cols + j];
          piv_col = j;
        }
      }
      if (fabs(max_col_value) <= 1e-10) {
        /* This row is (close to) zero. */
        piv_row += 1; 
      }
    }
    if (fabs(max_col_value) <= 1e-10 && piv_row == rows) {
      printf("aca_rkmatrix: no pivot row found, rk-matrix has rank %d\n", r->kt);
      return r;
    }
    else {
      /* Compute all entries of column piv_col and determine new piv_row. */
      for (j = 0; j < rows; j++) {
	      /*************/
        /* exercise */
        /*************/
 
        r->a[i*rows + j] = compute_entry_aca(r, d, i, j, piv_col, x + d*j, 
                                             y + d*piv_col, test_function);
        if (fabs(max_col_value) < fabs(r->a[i*rows + j])) {
          max_col_value = r->a[i*rows + j];
          piv_row = j;
        }
      }
      /* Scale column i by 1.0 / max_col_value. */
      max_col_value = 1.0 / max_col_value;
      dscal_(&rows, &max_col_value, r->a + i*rows, eins_);
      /* for (j = 0; j < rows; j++) {
        printf("r->a(%d, %d) = %f, piv_row = %d, piv_col = %d\n", j, i, r->a[i*rows + j], piv_row, piv_col);
      } */

      /* Compute all entries of row piv_row. */
      for (j = 0; j < cols; j++) { 
        r->b[i*cols + j] = compute_entry_aca(r, d, i, piv_row, j, x + d*piv_row,
                                             y + d*j, test_function);
      }
      /*
      for (j = 0; j < cols; j++) { 
        printf("r->b(%d, %d) = %f\n", i, j, r->b[i*cols + j]);
      }
      */
    }
    /* printf("piv_row = %d, piv_col = %d\n", piv_row, piv_col); */
    r->kt += 1;
  }

  /*
  printf("Full matrix block\n");
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      printf("%f ", test_function(d, &x[d*i], &y[d*j]));
    }
    printf("\n");
  }
  printf("\nComputed U, V factors\n");
  for (j = 0; j < rows; j++) {
    for (i = 0; i < r->kt; i++) {
      printf("%f ", r->a[i*rows + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < cols; j++) {
    for (i = 0; i < r->kt; i++) { 
      printf("%f ", r->b[i*cols + j]);
    }
    printf("\n");
  }
  */
  del_doubles(x);
  del_doubles(y);
  
  return r;
}

prkmatrix
aca_delta_rkmatrix(int d, int k, double delta, int rows, int cols, 
                   const double *nodes_x, const double *nodes_y, 
                   double (*test_function)(int d, const double *x, const double *y)) {
  int i, j, piv_row, piv_col;
  double max_col_value, u0, v0, ui, vi;
  double *x, *y;
  prkmatrix r;

  assert(nodes_x != NULL || d == 0 || rows == 0);
  assert(nodes_y != NULL || d == 0 || cols == 0);
  
  k = min(rows, k);
  k = min(cols, k);

  r = new_zero_rkmatrix(k, rows, cols);

  x = allocate_doubles(d * rows);
  y = allocate_doubles(d * cols);
  
  /* The rows d-dimensional nodes x_i are saved in nodes_x and the 
     cols d-dimensional nodes y_i are saved in nodes_y. 
     nodes_x and nodes_y begin with all first components, then all second
     components etc. We reorder the entries of nodes_x and nodes_y to begin
     with all d components of the first node, then all d components of the 
     second one, etc. and save them in x and y. */
  for (j = 0; j < d; j++) {
    for (i = 0; i < rows; i++) {
      x[d*i + j] = nodes_x[j*rows + i];
    }
  }
  for (j = 0; j < d; j++) {
    for (i = 0; i < cols; i++) {
      y[d*i + j] = nodes_y[j*cols + i];
    }
  }

  for (i = 0; i < k; i++) {
    /* Compute all entries of a pivot row. */
    piv_row = 0;
    max_col_value = 0.0;
    while (fabs(max_col_value) <= 1e-10 && piv_row < rows) {
      for (j = 0; j < cols; j++) {
        r->b[i*cols + j] = compute_entry_aca(r, d, i, piv_row, j, x + d*piv_row,
                                             y + d*j, test_function);
        if (fabs(max_col_value) < fabs(r->b[i*cols + j])) {
          max_col_value = r->b[i*cols + j];
          piv_col = j;
        }
      }
      if (fabs(max_col_value) <= 1e-10) {
        /* This row is (close to) zero. */
        piv_row += 1; 
      }
    }
    if (fabs(max_col_value) <= 1e-10 && piv_row == rows) {
      printf("aca_rkmatrix: no pivot row found, rk-matrix has rank %d\n", r->kt);
      return r;
    }
    else {
      /* Compute all entries of column piv_col and determine new piv_row. */
      for (j = 0; j < rows; j++) {
        r->a[i*rows + j] = compute_entry_aca(r, d, i, j, piv_col, x + d*j, 
                                             y + d*piv_col, test_function);
        if (fabs(max_col_value) < fabs(r->a[i*rows + j])) {
          max_col_value = r->a[i*rows + j];
          piv_row = j;
        }
      }
      /* Scale column i by 1.0 / max_col_value. */
      max_col_value = 1.0 / max_col_value;
      dscal_(&rows, &max_col_value, r->a + i*rows, eins_);
      /* for (j = 0; j < rows; j++) {
        printf("r->a(%d, %d) = %f, piv_row = %d, piv_col = %d\n", j, i, r->a[i*rows + j], piv_row, piv_col);
      } */

      /* Compute all entries of row piv_row. */
      for (j = 0; j < cols; j++) { 
        r->b[i*cols + j] = compute_entry_aca(r, d, i, piv_row, j, x + d*piv_row,
                                             y + d*j, test_function);
      }
      /*
      for (j = 0; j < cols; j++) { 
        printf("r->b(%d, %d) = %f\n", i, j, r->b[i*cols + j]);
      }
      */
      
      /* Check heuristic stopping criterion for adaptive rank determination. */
      if (i == 0) { 
        u0 = sqrt(ddot_(&rows, r->a, eins_, r->a, eins_));
        v0 = sqrt(ddot_(&cols, r->b, eins_, r->b, eins_));
      }
      else {
        ui = sqrt(ddot_(&rows, r->a + i*rows, eins_, r->a + i*rows, eins_));
        vi = sqrt(ddot_(&cols, r->b + i*cols, eins_, r->b + i*cols, eins_));
        if (ui * vi / (u0*v0) < delta) {
          r->kt += 1;
          return r;
        }
      }
      
    }
    /* printf("piv_row = %d, piv_col = %d\n", piv_row, piv_col); */
    r->kt += 1;
  }

  /*
  printf("Full matrix block\n");
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      printf("%f ", test_function(d, &x[d*i], &y[d*j]));
    }
    printf("\n");
  }
  printf("\nComputed U, V factors\n");
  for (j = 0; j < rows; j++) {
    for (i = 0; i < r->kt; i++) {
      printf("%f ", r->a[i*rows + j]);
    }
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < cols; j++) {
    for (i = 0; i < r->kt; i++) { 
      printf("%f ", r->b[i*cols + j]);
    }
    printf("\n");
  }
  */
  del_doubles(x);
  del_doubles(y);
  
  return r;
}
