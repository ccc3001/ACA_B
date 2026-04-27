#ifndef Interpolation
#define Interpolation

#include "rkmatrix.h"

/* Evaluate the nu'th Lagrange-polynomial of degree degree for 
   the nodes nodes at x. */ 
double lagrange_1d(int nu, int degree, const double *nodes, double x);

/* Evaluate the (nu_1, ..., nu_d)'th d-variate Lagrange-polynomial 
   of degree degree (per spatial dimension) for the nodes nodes at x. */ 
double lagrange_md(int d, const int *nu, int degree, const double *nodes, const double *x);

/* Evaluate d-variate test function ln(-|x-y|_2) for x != y and 0 for x == y. */
double test_function_log(int d, const double *x, const double *y);

/* Evaluate d-variate test function e^(-|x-y|_2^2). */
double test_function_gaussian(int d, const double *x, const double *y);

/* Create an rkmatrix of rank k = (degree + 1)^d that approximates the matrix
   with entries a_ij = test_function(nodes_x[i], nodes_y[j]) using an
   interpolation wrt x of the specified degree (per spatial dimension). 
   rows denotes the number of nodes in nodes_x and 
   cols denotes the number of nodes in nodes_y. Nodes in nodes_x and nodes_y
   are saved columnwise for spatial dimension, i.e. first all
   first-coordinates, then all second-coordinates, etc. */
prkmatrix
interpolate_rkmatrix(int d, int rows, int cols, const double *nodes_x, const double *nodes_y,
                     int degree, double (*test_function)(int d, const double *x, const double *y)); 

#endif
