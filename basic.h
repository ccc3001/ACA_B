#ifndef Basic
#define Basic

#include <stdlib.h>

/* Allocate storage for n doubles. */
double*
allocate_doubles(int n);

/* Allocate storage for n ints. */
int*
allocate_ints(int n);

/* Delete doubles. */
void
del_doubles(double *x);

/* Delete ints. */
void
del_ints(int *j);

/* Print the first n entries of x. */
void
print_entries(int n, double *x);

/* Copy the first n entries of x to y. */
void
copy_entries(int n, const double *x, double *y);

/* Set the first n entries of x to v. */
void
fill_entries(int n, double *x, double v);

/* Set the first n entries of x to random values in [0, 1]. */
void
fill_random_entries(int n, double *x, int seed);

/* Set the first n entries of x to 0.0. */
void
clear_entries(int n, double *x);

/* Convert size in bytes to megabytes. */
double
convertb2_mb(size_t size);

/* Compute max(a, b). */
int
max(int a, int b);

/* Compute max(x, y). */
double
max_double(double x, double y);

/* Compute min(a, b). */
int
min(int a, int b);

/* Compute min(x, y). */
double
min_double(double x, double y);

#endif
