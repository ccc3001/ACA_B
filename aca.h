#ifndef ACA
#define ACA

#include "rkmatrix.h"

/* Compute test_function(d, x, y) - sum_{i=0}^{k-1} r->a[i*r->rows + row] * r->b[i*r->cols + col]. */
double 
compute_entry_aca(pcrkmatrix r, int d, int k, int row, int col, 
                  const double *x, const double *y, 
                  double (*test_function)(int d, const double *x, const double *y));

/* Create an rkmatrix of rank <= min{k, rows, cols} that approximates the matrix
   with entries a_ij = test_function(nodes_x[i], nodes_y[j]) using ACA.
   rows denotes the number of nodes in nodes_x and 
   cols denotes the number of nodes in nodes_y. Nodes in nodes_x and nodes_y
   are saved columnwise for spatial dimension, i.e. first all
   first-coordinates, then all second-coordinates, etc. */
prkmatrix
aca_rkmatrix(int d, int k, int rows, int cols, 
             const double *nodes_x, const double *nodes_y, 
             double (*test_function)(int d, const double *x, const double *y)); 

/* Create an rkmatrix of rank <= min{k, rows, cols} that approximates the matrix
   with entries a_ij = test_function(nodes_x[i], nodes_y[j]) using ACA.
   rows denotes the number of nodes in nodes_x and 
   cols denotes the number of nodes in nodes_y. Nodes in nodes_x and nodes_y
   are saved columnwise for spatial dimension, i.e. first all
   first-coordinates, then all second-coordinates, etc.
   This function is a modification of the function aca_rkmatrix. It is 
   extended by a heuristic stopping criterion and
   delta is the parameter for this heuristic stopping criterion. */
prkmatrix
aca_delta_rkmatrix(int d, int k, double delta, int rows, int cols, 
                   const double *nodes_x, const double *nodes_y, 
                   double (*test_function)(int d, const double *x, const double *y)); 

#endif
