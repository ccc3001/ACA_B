#include "fullmatrix.h"
#include "basic.h"
#include "blas.h"

#include <stdio.h>
#include <assert.h>

pfullmatrix
new_fullmatrix(int rows, int cols) {
  pfullmatrix f;
  
  assert(rows >= 0);
  assert(cols >= 0);
  
  f = (pfullmatrix) malloc(sizeof(fullmatrix));
  if (f == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
  
  f->rows = rows;
  f->cols = cols;
  f->e = allocate_doubles(rows * cols);
  
  return f;
}