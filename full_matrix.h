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
