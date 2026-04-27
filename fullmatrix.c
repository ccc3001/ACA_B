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

pfullmatrix 
new_random_fullmatrix(int rows, int cols, int seed) {
  pfullmatrix f;
  
  f = new_fullmatrix(rows, cols);
  fill_random_entries(rows * cols, f->e, seed);
  
  return f;
}

pfullmatrix
new_zero_fullmatrix(int rows, int cols) {
  pfullmatrix f;
  
  f = new_fullmatrix(rows, cols);
  clear_entries(rows * cols, f->e);
  
  return f;
}

void
del_fullmatrix(pfullmatrix f) {
  assert(f != NULL);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert(f->rows == 0 || f->cols == 0 || f->e != NULL);
  
  del_doubles(f->e);
  free(f);
}

void
print_fullmatrix(pcfullmatrix f) { 
  int i, j;

  assert(f != NULL);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert(f->rows == 0 || f->cols == 0 || f->e != NULL);
  
  printf("%dx%d fullmatrix\n", f->rows, f->cols);
  for (i = 0; i < f->rows; i++) {
    for (j = 0; j < f->cols; j++) {
      printf("%.4g ", f->e[j*f->rows + i]);
    }
    printf("\n");
  }
  printf("\n");
}

size_t
getsize_fullmatrix(int rows, int cols) {
  size_t size;
  
  assert(rows >= 0);
  assert(cols >= 0);
  
  size = sizeof(fullmatrix); 
  size += sizeof(double) * rows * cols;
  
  return size;
}

void
fill_fullmatrix(pfullmatrix f, int d, const double *nodes_x, const double *nodes_y, 
                double (*test_function)(int d, const double *x, const double *y)) {
  int i, j;
  double *x, *y;

  assert(f != NULL);
  assert(d >= 0);
  assert((f->rows > 0 && d > 0) || nodes_x == NULL);
  assert(f->rows == 0 || d == 0 || nodes_x != NULL);
  assert((f->cols > 0 && d > 0) || nodes_y == NULL);
  assert(f->cols == 0 || d == 0 || nodes_y != NULL);
  
  x = allocate_doubles(d * f->rows);
  y = allocate_doubles(d * f->cols);
  
  /* Reorder and copy the entries of nodes_x and nodes_y such that the 
     entries of x and y begin with all d components of the first node, 
     then all d components of the second node, etc. */
  for (j = 0; j < d; j++) {
    for (i = 0; i < f->rows; i++) {
      x[d*i + j] = nodes_x[j*f->rows + i];
    }
    for (i = 0; i < f->cols; i++) {
      y[d*i + j] = nodes_y[j*f->cols + i];
    }
  }
  
  for (i = 0; i < f->cols; i++) {
    for (j = 0; j < f->rows; j++) {
      f->e[i*f->rows + j] = test_function(d, x + d*j, y + d*i);
    }
  }

  del_doubles(x);
  del_doubles(y);
}

#ifdef USE_BLAS 
void
add_fullmatrix(pcfullmatrix f, pfullmatrix g, double a) {
  int numentries;
  
  assert(f != NULL);
  assert(g != NULL);
  assert(f->rows == g->rows);
  assert(f->cols == g->cols);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || (f->e == NULL && g->e == NULL) );
  assert( f->rows == 0 || f->cols == 0 || (f->e != NULL && g->e != NULL) );
  
  numentries = f->rows * f->cols;
  
  if (numentries > 0) {
    daxpy_(&numentries, &a, f->e, eins_, g->e, eins_);
  }
}
#else
void
add_fullmatrix(pcfullmatrix f, pfullmatrix g, double a) {
  int i;
  
  assert(f != NULL);
  assert(g != NULL);
  assert(f->rows == g->rows);
  assert(f->cols == g->cols);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || (f->e == NULL && g->e == NULL) );
  assert( f->rows == 0 || f->cols == 0 || (f->e != NULL && g->e != NULL) );
  
  for (i = 0; i < f->rows * f->cols; i++) {
    g->e[i] += a * f->e[i];
  } 
}
#endif

#ifdef USE_BLAS 
void 
addeval_fullmatrix(pcfullmatrix f, const double *v, double *w) {
  assert(f != NULL);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert(f->rows == 0 || f->cols == 0 || f->e != NULL);
  assert(f->rows == 0 || w != NULL);
  assert(f->rows > 0 || w == NULL);
  assert(f->cols == 0 || v != NULL);
  assert(f->cols > 0 || v == NULL);
  
  if(f->rows > 0 && f->cols > 0) {
    dgemv_(ntrans, &f->rows, &f->cols, deins_, f->e, &f->rows,
	       v, eins_, deins_, w, eins_);  
  }
}
#else
void 
addeval_fullmatrix(pcfullmatrix f, const double *v, double *w) {
  int i, j;
  
  assert(f != NULL);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert(f->rows == 0 || f->cols == 0 || f->e != NULL);
  assert(f->rows == 0 || w != NULL);
  assert(f->rows > 0 || w == NULL);
  assert(f->cols == 0 || v != NULL);
  assert(f->cols > 0 || v == NULL);
  
  for (i = 0; i < f->rows; i++) {
    for (j = 0; j < f->cols; j++) { 
      w[i] += f->e[i + j*f->rows] * v[j];
    }
  } 
}
#endif

void 
eval_fullmatrix(pcfullmatrix f, const double *v, double *w) {
  clear_entries(f->rows, w);
  addeval_fullmatrix(f, v, w);     
}

#ifdef USE_BLAS
void
mul_fullmatrix(pcfullmatrix f, pcfullmatrix g, pfullmatrix fg) { 
  assert(f != NULL);
  assert(g != NULL);
  assert(fg != NULL);
  assert(f->rows == fg->rows);
  assert(f->cols == g->rows);
  assert(fg->cols == g->cols);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert(g->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert( (f->rows == 0 || f->cols == 0) || f->e != NULL );
  assert( (g->rows > 0 && g->cols > 0) || g->e == NULL );
  assert( (g->rows == 0 || g->cols == 0) || g->e != NULL );
  assert( (fg->rows > 0 && fg->cols > 0) || fg->e == NULL );
  assert( (fg->rows == 0 || fg->cols == 0) || fg->e != NULL );
  
  if (f->rows > 0 && g->cols > 0) {
    if (f->cols > 0) {
      dgemm_(ntrans, ntrans, &f->rows, &g->cols, &g->rows, deins_, f->e,
             &f->rows, g->e, &g->rows, dnull_, fg->e, &fg->rows);
    }
    else {
      clear_entries(f->rows * g->cols, fg->e);
    }
  }
}
#else
void
mul_fullmatrix(pcfullmatrix f, pcfullmatrix g, pfullmatrix fg) {
  int i, j, k, ind_ij;
  
  assert(f != NULL);
  assert(g != NULL);
  assert(fg != NULL);
  assert(f->rows == fg->rows);
  assert(f->cols == g->rows);
  assert(fg->cols == g->cols);
  assert(f->rows >= 0);
  assert(f->cols >= 0);
  assert(g->cols >= 0);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert( (f->rows == 0 || f->cols == 0) || f->e != NULL );
  assert( (g->rows > 0 && g->cols > 0) || g->e == NULL );
  assert( (g->rows == 0 || g->cols == 0) || g->e != NULL );
  assert( (fg->rows > 0 && fg->cols > 0) || fg->e == NULL );
  assert( (fg->rows == 0 || fg->cols == 0) || fg->e != NULL );

  clear_entries(f->rows * g->cols, fg->e);
  
  for (i = 0; i < f->rows; i++) {
    for (j = 0; j < g->cols; j++) {
      ind_ij = i + j*f->rows;
      for (k = 0; k < f->cols; k++) {
        fg->e[ind_ij] += f->e[i + k*f->rows] * g->e[k + j*f->cols];
      }
    }
  }
}
#endif
