#include <stdbool.h>
#include <stdlib.h>

int matrix_aca_test(int d, int L, bool runtime_benchmark, const char *test_function_name, int iter, int repeats, double eps, bool residuals_benchmark, double eps_residual, int starts_res, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y, double (*test_function)(int dim, const double *x, const double *y), bool hbaca_residual_test, double eps_hbaca, int h_baca_repeats, int seed);