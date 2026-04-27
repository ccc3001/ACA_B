#ifndef SuperMatrix
#define SuperMatrix

#include "fullmatrix.h"
#include "rkmatrix.h"

#include <stdlib.h>

typedef struct _supermatrix supermatrix;
typedef supermatrix *psupermatrix;
typedef const supermatrix *pcsupermatrix;

struct _supermatrix {
  int block_rows;           
  int block_cols;
  int rows;              
  int cols;
  prkmatrix r;        
  pfullmatrix f;      
  psupermatrix *s;    
};

/* Create a new supermatrix. */
psupermatrix 
new_supermatrix(int block_rows, int block_cols, int rows, int cols,
                pcrkmatrix r, pcfullmatrix f);

/* Create a new supermatrix (Hp-matrix structure with rank k) and
   don't set entries. */
psupermatrix 
new_hp_supermatrix(int p, int k);

/* Delete supermatrix. */
void
del_supermatrix(psupermatrix s);

/* Compute size of storage in bytes to save the supermatrix s. */
size_t
getsize_supermatrix(pcsupermatrix s);

/* Set all entries of s to v. */
void 
fill_supermatrix(psupermatrix s, double v);

/* Set all entries of s to random values in [0, 1]. */
void 
fill_random_supermatrix(psupermatrix s, int seed);

/* Fill a supermatrix hp (works only for Hp-matrix structure) to approximate
   the matrix with entries a_ij = test_function(1, nodes_x[i], nodes_x[j]) 
   using an interpolation wrt (one-dimensional) x of the specified degree to
   fill the rk offdiagonal blocks of rank k = degree + 1. */
void
fill_kernelinterpolation_hp_supermatrix(psupermatrix hp, int degree, const double *nodes_x, 
                                        double (*test_function)(int d, const double *x, const double *y));

/* Fill a supermatrix s to approximate the matrix with entries
   a_ij = test_function(d, x_i, x_j) using an interpolation wrt the first
   variable of the specified degree to fill the rk blocks of rank
   k = (degree + 1)^d. The fullmatrix blocks are filled by the function
   fill_fullmatrix. The d-dimensional nodes x_i are saved in nodes_x. 
   nodes_x begins with all first components, then all second
   components etc. offset_row, offset_col and total_size are used to call
   this function recursively. At the beginning of the recursion this
   function should be called with offset_row = offset_col = 0 and 
   total_size = max(s->rows, s->cols). */
void
fill_kernelinterpolation_supermatrix(psupermatrix s, int d, int degree, const double *nodes_x, 
                                     double (*test_function)(int d, const double *x, const double *y), 
                                     int offset_row, int offset_col, int total_size);

/* Fill a supermatrix s to approximate the matrix with entries
   a_ij = test_function(d, x_i, x_j) using ACA to fill the rk blocks of
   rank k with min(s->rows, s->cols) > k. Other rkmatrix blocks are changed
   to fullmatrix blocks. The fullmatrix blocks are filled by the function
   fill_fullmatrix. The d-dimensional nodes x_i are saved in nodes_x. 
   nodes_x begins with all first components, then all second
   components etc. offset_row, offset_col and total_size are used to call
   this function recursively. At the beginning of the recursion this
   function should be called with offset_row = offset_col = 0 and 
   total_size = max(s->rows, s->cols). */
void
fill_kernelaca_supermatrix(psupermatrix s, int d, int k, const double *nodes_x, 
                           double (*test_function)(int d, const double *x, const double *y), 
                           int offset_row, int offset_col, int total_size);

/* Fill a supermatrix s to approximate the matrix with entries
   a_ij = test_function(d, x_i, x_j) using ACA to fill the rk blocks of
   rank k with min(s->rows, s->cols) > k. Other rkmatrix blocks are changed
   to fullmatrix blocks. The fullmatrix blocks are filled by the function
   fill_fullmatrix. The d-dimensional nodes x_i are saved in nodes_x. 
   nodes_x begins with all first components, then all second
   components etc. offset_row, offset_col and total_size are used to call
   this function recursively. At the beginning of the recursion this
   function should be called with offset_row = offset_col = 0 and 
   total_size = max(s->rows, s->cols).
   This function is a modification of the function 
   fill_kernelaca_supermatrix. It calls the function aca_delta_rkmatrix
   instead of the function aca_rkmatrix and
   delta is the parameter for the heuristic stopping criterion. */
void
fill_kernelaca_delta_supermatrix(psupermatrix s, int d, int k, double delta, const double *nodes_x, 
                                 double (*test_function)(int d, const double *x, const double *y), 
                                 int offset_row, int offset_col, int total_size);

/* Set all entries of s to 0.0. */
void 
clear_supermatrix(psupermatrix s);

/* Compute w := w + S * v. */
void 
addeval_supermatrix(pcsupermatrix s, const double *v, double *w);

/* Compute w := S * v. */
void 
eval_supermatrix(pcsupermatrix s, const double *v, double *w);

/* Recompress a supermatrix by calling the function adaptive_truncate_rkmatrix
   (to truncate an rkmatrix such that the relative error in the 2-norm is smaller or equal to eps)
   for each rkmatrix-block. */
void
recompress_supermatrix(psupermatrix s, double eps);

/* Write the structure of s on the console for row_offset = col_offset = 0. */
void 
write_supermatrix(pcsupermatrix s, int row_offset, int col_offset);

/* Create a plot of the structure of s. The name of the file is given by fname. */
void
output_supermatrix(pcsupermatrix s, char *fname);

#endif
