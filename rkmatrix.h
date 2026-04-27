#ifndef RkMatrix
#define RkMatrix

#include "fullmatrix.h"

#include <stdlib.h>

typedef struct _rkmatrix rkmatrix;
typedef rkmatrix *prkmatrix;
typedef const rkmatrix *pcrkmatrix;

/* rkmatrix R = A * B^T */
struct _rkmatrix {
  int k; /* maximal rank */   
  int kt; /* current rank */      
  int rows;       
  int cols;
  double *a; /* entries of A in column-major order */
  double *b; /* entries of B in column-major order */
};

/* Create a new rkmatrix, set k := k, kt := 0 and don't set entries. */
prkmatrix
new_rkmatrix(int k, int rows, int cols);

/* Create a new rkmatrix and set k := k, kt := k
   and all entries to random values in [0, 1]. */
prkmatrix 
new_random_rkmatrix(int k, int rows, int cols, int seed);

/* Create a new rkmatrix and set k := k, kt := 0 and all entries to 0.0. */
prkmatrix 
new_zero_rkmatrix(int k, int rows, int cols);

/* Delete rkmatrix. */
void
del_rkmatrix(prkmatrix r);

/* Print rkmatrix. */
void
print_rkmatrix(pcrkmatrix r);

/* Compute size of storage in bytes to save a rows x cols rkmatrix of rank k. */
size_t
getsize_rkmatrix(int k, int rows, int cols);

/* Compute w := w + R * v = w + A * B^T * v. */
void 
addeval_rkmatrix(pcrkmatrix r, const double *v, double *w);

/* Compute w := R * v = A * B^T * v. */
void 
eval_rkmatrix(pcrkmatrix r, const double *v, double *w);

/* Convert rkmatrix to fullmatrix. */
void
convertrk2_fullmatrix(pcrkmatrix r, pfullmatrix f);

/* Compute the truncation of r to rank k. If k > r->kt,
   then just copy r into rtrunc. The rest of the entries
   of rtrunc->a and rtrunc->b are not set.*/
prkmatrix
truncate_rkmatrix(pcrkmatrix r, int k);

/* Truncate r such that the relative error in the 2-norm is smaller or equal to eps. */
void
adaptive_truncate_rkmatrix(prkmatrix r, double eps);

/* Compute and return all singular values of r. */
double* 
getsigma_rkmatrix(pcrkmatrix r);

/* Compute and return the nr-th singular value of r.
   Returns 0.0, if nr > min(r->rows, r->cols). */
double
getsv_rkmatrix(pcrkmatrix r, int nr);

#endif
