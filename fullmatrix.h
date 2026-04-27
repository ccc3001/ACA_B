#ifndef FullMatrix
#define FullMatrix

#include <stdlib.h>

typedef struct _fullmatrix fullmatrix;
typedef fullmatrix *pfullmatrix;
typedef const fullmatrix *pcfullmatrix;

struct _fullmatrix {
  int rows;
  int cols;
  double *e; /* entries in column-major order */
};

/* Create a new fullmatrix and don't set entries. */
pfullmatrix 
new_fullmatrix(int rows, int cols);

/* Create a new fullmatrix and set all entries to random values in [0, 1]. */
pfullmatrix 
new_random_fullmatrix(int rows, int cols, int seed);

/* Create a new fullmatrix and set all entries to 0.0. */
pfullmatrix 
new_zero_fullmatrix(int rows, int cols);

/* Delete fullmatrix. */
void
del_fullmatrix(pfullmatrix f);

/* Print f. */
void
print_fullmatrix(pcfullmatrix f);

/* Compute size of storage in bytes to save a rows x cols fullmatrix. */
size_t
getsize_fullmatrix(int rows, int cols);

/* Fill matrix f such that f_ij = test_function(d, x_i, y_j). The d-dimensional
   nodes x_i, y_j are saved in nodes_x and nodes_y. nodes_x and nodes_y begin
   with all first components, then all second components etc. */
void
fill_fullmatrix(pfullmatrix f, int d, const double *nodes_x, const double *nodes_y, 
                double (*test_function)(int d, const double *x, const double *y));

/* Compute G := G + a * F. */ 
void
add_fullmatrix(pcfullmatrix f, pfullmatrix g, double a);

/* Compute w := w + F * v. */
void 
addeval_fullmatrix(pcfullmatrix f, const double *v, double *w);

/* Compute w := F * v. */
void 
eval_fullmatrix(pcfullmatrix f, const double *v, double *w);

/* Compute FG := F * G. */
void
mul_fullmatrix(pcfullmatrix f, pcfullmatrix g, pfullmatrix fg); 

#endif
