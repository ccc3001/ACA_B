#include "rkmatrix.h"
#include "fullmatrix.h"

double
compute_entry_aca_new(pcrkmatrix r, int k, int row, int col, pcfullmatrix A);

prkmatrix
aca_rkmatrix_new(double eps, pcfullmatrix A);

prkmatrix
aca_rkmatrix_(double eps, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y, double (*test_function)(int dim, const double *x, const double *y), double **residuals_u, double **residuals_v);