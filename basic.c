#include "basic.h"
#include "blas.h"

#include <stdio.h>
#include <assert.h>

double*
allocate_doubles(int n) {
  double *x;

  assert(n >= 0);
  
  if (n == 0) {
    x = NULL;
  }
  else {
    x = (double*) malloc(sizeof(double) * n);
    if (x == NULL) {
      printf("Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }
  }

  return x;
}

int*
allocate_ints(int n) {
  int *j;

  assert(n >= 0);
  
  if (n == 0) {
    j = NULL;
  }
  else {
    j = (int*) malloc(sizeof(int) * n);
    if (j == NULL) {
      printf("Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }
  }

  return j;
}

void
del_doubles(double *x) {
  if (x != NULL) {
    free(x);
  }
}

void
del_ints(int *j) {
  if (j != NULL) {
    free(j);
  }
}

void
print_entries(int n, double *x) {
  int i;
  
  assert(n >= 0);
  assert(n == 0 || x != NULL);
  
  for (i = 0; i < n; i++) {
    printf("%.4g\n", x[i]);
  }
  printf("\n");
}

#ifdef USE_BLAS
void
copy_entries(int n, const double *x, double *y) {
  assert(n >= 0);
  assert( n == 0 || (x != NULL && y != NULL) );
  
  dcopy_(&n, x, eins_, y, eins_);
}
#else
void
copy_entries(int n, const double *x, double *y) {
  int i;
  
  assert(n >= 0);
  assert( n == 0 || (x != NULL && y != NULL) );
  
  for (i = 0; i < n; i++) {
    y[i] = x[i];
  }
}
#endif

void
fill_entries(int n, double *x, double v) {
  int i;

  assert(n >= 0);
  assert(n == 0 || x != NULL);
  
  for (i = 0; i < n; i++) {
    x[i] = v;
  }
}

void
fill_random_entries(int n, double *x, int seed) {
  int i;
  
  assert(n >= 0);
  assert(n == 0 || x != NULL);
  
  /* Seed the random number generator used by the function rand. */
  srand((unsigned) seed);

  for (i = 0; i < n; i++) {
    /* rand returns an int pseudo-random number in {0, ..., RAND_MAX}. */
    x[i] = (double) rand() / (double) RAND_MAX;
  }
}

void
clear_entries(int n, double *x) {
  fill_entries(n, x, 0.0);
}

double
convertb2_mb(size_t size) {
  return (double) size / (1024.0 * 1024.0);
}

int
max(int a, int b) {
  return a > b ? a : b;
}

double
max_double(double x, double y) {
  return x > y ? x : y;
}

int
min(int a, int b) {
  return a > b ? b : a;
}

double
min_double(double x, double y) {
  return x > y ? y : x;
}
