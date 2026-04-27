#include "interpolation.h"
#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double 
lagrange_1d(int nu, int degree, const double *nodes, double x) {
  int i;
  double y = 1.0;

  assert(nu >= 0);
  assert(degree >= nu);
  assert(nodes != NULL || degree == 0);
  
  for (i = 0; i <= degree; i++) {
    if (i != nu) { 
      y *= (x - nodes[i]) / (nodes[nu] - nodes[i]);
    }
  }

  return y;
}

double 
lagrange_md(int d, const int *nu, int degree, const double *nodes, const double *x) {
  int i;
  double y = 1.0;

  assert(d >= 0);
  assert( (nu != NULL && x != NULL) || d == 0);
  
  for (i = 0; i < d; i++) {
    y *= lagrange_1d(nu[i], degree, &nodes[i * (degree+1)], x[i]);
  }

  return y;
}

double 
test_function_log(int d, const double *x, const double *y) {
  int i;
  double r = 0.0;
  
  assert(d >= 0);
  assert( (x != NULL && y != NULL) || d == 0);
  
  for (i = 0; i < d; i++) {
    r += (x[i] - y[i]) * (x[i] - y[i]);
  }
  
  r = sqrt(r);

  return r < 1e-10 ? 0.0 : -log(r);
}

double 
test_function_gaussian(int d, const double *x, const double *y) {
  int i;
  double r = 0.0;
  
  assert(d >= 0);
  assert( (x != NULL && y != NULL) || d == 0);
  
  for (i = 0; i < d; i++) {
    r += (x[i] - y[i]) * (x[i] - y[i]);
  }

  return exp(-r);
}

prkmatrix
interpolate_rkmatrix(int d, int rows, int cols, const double *nodes_x, const double *nodes_y, 
                     int degree, double (*test_function)(int d, const double *x, const double *y)) {
  prkmatrix r;
  int k, i, j;
  int *nu;
  double *nodes_cheb, *nodes_int, *xi, *yi;
  double x_min, x_max, m, b;
  
  assert(nodes_x != NULL || d == 0 || rows == 0);
  assert(nodes_y != NULL || d == 0 || cols == 0);

  /* Create new rkmatrix of rank k = (degree+1)^d. */
  k = (int) pow((double) degree + 1, (double) d); 
  r = new_rkmatrix(k, rows, cols);
  r->kt = k;

  /* Create (degree + 1) Chebyshev nodes (roots of T_{degree+1})
     on [-1, 1] (see lemma 4.5). */
  nodes_cheb = allocate_doubles(degree + 1);
  for (i = 0; i <= degree; i++) {
    nodes_cheb[i] = cos( M_PI * (2.0*i + 1.0) / (2.0 * (degree+1)) );
  }
  
  /* Create interpolation nodes nodes_int. Nodes for interpolation polynomial
     are saved columnwise for spatial dimension in nodes_int, i.e.
     first all (degree+1) first-coordinates, then all second-coordinates,
     etc. degree+1 Chebyshev nodes nodes_cheb are linearly transformed
     onto intervals containing nodes_x, see (4.3). */ 
  nodes_int = allocate_doubles((degree + 1) * d);
  for (i = 0; i < d; i++) {
    /* Compute interval [a_i, b_i] := [x_min, x_max]. */
    x_min = nodes_x[i * rows];
    x_max = x_min;
    for (j = 1; j < rows; j++) {
      x_min = min_double(x_min, nodes_x[i*rows + j]);
      x_max = max_double(x_max, nodes_x[i*rows + j]);
    } 
    /* Catch the artificial case of only a single row,
       avoid m < 5e-13. */
    if (x_max - x_min < 1e-12) { 
      x_max += 1e-12;
    }
    /* Transform Chebyshev nodes from [-1, 1] onto
       [a_i, b_i] := [x_min, x_max] (see (4.3)). */
    m = 0.5 * (x_max - x_min);
    b = 0.5 * (x_max + x_min);
    /* Compute all i-coordinates of degree+1 interpolation nodes. */
    for (j = 0; j <= degree; j++) {
      nodes_int[i*(degree+1) + j] = m * nodes_cheb[j] + b;
    }
  }

  /* Use interpolation to fill rkmatrix r. */
  nu = allocate_ints(d);
  for (i = 0; i < d; i++) {
    nu[i] = 0;
  }
  if (d == 1) {
    for (nu[0] = 0; nu[0] <= degree; nu[0]++) {
      for (i = 0; i < rows; i++) {
        r->a[nu[0]*rows + i] = lagrange_1d(nu[0], degree, nodes_int, nodes_x[i]);
      }
      for (i = 0; i < cols; i++) {
        r->b[nu[0]*cols + i] = test_function(d, nodes_int + nu[0], nodes_y + i);
      }
    }
  }
  else if (d > 1) {
    xi = allocate_doubles(d);
    yi = allocate_doubles(d);
    
    for (j = 0; j < r->kt; j++) {
      for (i = 0; i < rows; i++) {
        /* Compute i-th node in nodes_x. */
        for (k = 0; k < d; k++) { 
          xi[k] = nodes_x[i + k*rows];
        }
        r->a[j*rows + i] = lagrange_md(d, nu, degree, nodes_int, xi);
      }
      for (i = 0; i < cols; i++) {
        /* Compute i-th node in nodes_y and corresponding values in
           nodes_int (\xi_{\nu} in 4.2.2) to evaluate the test_function. */
        for (k = 0; k < d; k++) {
          yi[k] = nodes_y[i + k*cols];
          xi[k] = nodes_int[nu[k] + k * (degree+1)];
	    }
        r->b[j*cols + i] = test_function(d, xi, yi);
      }
      
      /*
      printf("%d :", j + 1);
      for (i = 0; i < d; i++) {
        printf("%d ", nu[i]);
      }
      printf("\n");
      */

      /* Change nu to compute all values of r->a and r->b. */
      for (i = d - 1; i >= 0; i--) {
        if (nu[i] < degree) {
          nu[i]++;
          break;
	    }
        else {
          nu[i] = 0;
	    }
      }
    }
    
    del_doubles(xi);
    del_doubles(yi);
  }

  del_doubles(nodes_int);
  del_doubles(nodes_cheb);
  del_ints(nu);

  return r;
}
