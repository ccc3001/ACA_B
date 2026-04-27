#include "blas.h"
#include "basic.h"

#include <assert.h>

int eins_[1] = { 1 };
double deins_[1] = { 1.0 };
double dnull_[1] = { 0.0 };
double dmeins_[1] = { -1.0 };
char ttrans[1] = { 't' };
char ntrans[1] = { 'n' };

void
matrixmultrans2(int rows, int cols, int mid,
		        const double *a, const double *b, double *c) {
  assert(rows >= 0);
  assert(cols >= 0);
  assert(mid >= 0);
  assert( (rows > 0 && mid > 0) || a == NULL );
  assert(rows == 0 || mid == 0 || a != NULL);
  assert( (cols > 0 && mid > 0) || b == NULL );
  assert(cols == 0 || mid == 0 || b != NULL);
  assert( (rows > 0 && cols > 0) || c == NULL );
  assert(rows == 0 || cols == 0 || c != NULL);

  if (rows > 0 && cols > 0) {
    if (mid > 0) {
      dgemm_(ntrans, ttrans, &rows, &cols, &mid,
	         deins_, a, &rows, b, &cols, dnull_, c, &rows);
    }
    else {
      /* Return rows x cols zero matrix, if mid = 0. */
      clear_entries(rows * cols, c);
    }
  }
}

void
matrixmul(int rows, int cols, int mid,
	      const double *a, const double *b, double *c) {
  assert(rows >= 0);
  assert(cols >= 0);
  assert(mid >= 0);
  assert( (rows > 0 && mid > 0) || a == NULL );
  assert(rows == 0 || mid == 0 || a != NULL);
  assert( (mid > 0 && cols > 0) || b == NULL );
  assert(mid == 0 || cols == 0 || b != NULL);
  assert( (rows > 0 && cols > 0) || c == NULL );
  assert(rows == 0 || cols == 0 || c != NULL);

  if (rows > 0 && cols > 0) {
    if (mid > 0) {
      dgemm_(ntrans, ntrans, &rows, &cols, &mid,
	         deins_, a, &rows, b, &mid, dnull_, c, &rows);
    }
    else {
      /* Return rows x cols zero matrix, if mid = 0. */
      clear_entries(rows * cols, c);
    }
  }
}

void
qr_decomposition(const double *a, int n, int k, double *q, double *r) {
  int i, j, lwork, info;
  double *work, *tau;

  assert(n >= k);
  assert(k >= 0);
  assert( (n > 0 && k > 0) || (a == NULL && q == NULL) );
  assert( n == 0 || k == 0 || (a != NULL && q != NULL) );
  assert(k > 0 || r == NULL);
  assert(k == 0 || r != NULL);
  
  lwork = 8*k + n;
  work = allocate_doubles(lwork); 
  tau = allocate_doubles(n);
  
  /* Copy A, so A is not overwritten. */
  copy_entries(n * k, a, q);

  /* Compute the QR factorization A = Q*R, Q and R are saved in q */
  dgeqrf_(&n, &k, q, &n, tau, work, &lwork, &info);

  /* Compute the factor R of the QR factorization */
  for (i = 0; i < k; i++) {
    for (j = 0; j < k; j++) {
      /* R (k x k) is upper triangular part of q (n x k) */
      r[i + j*k] = j < i ? 0.0 : q[i + j*n];
    }
  }
  
  /* Compute the factor Q of the QR factorization */
  dorgqr_(&n, &k, &k, q, &n, tau, work, &lwork, &info);
   
  del_doubles(tau);
  del_doubles(work);
}

void
singular_value_decomposition(double *a, int n, double *u, double *s, double *v) {
  int i, j, info, lwork;
  double *work;
  double tmp;
  
  assert(n >= 0);
  assert( n > 0 || (a == NULL && u == NULL && s == NULL && v == NULL) );
  assert( n == 0 || (a != NULL && u != NULL && s != NULL && v != NULL) );
  
  lwork = 5 * n;
  work = allocate_doubles(lwork);
  
  /* Compute the singular value decomposition of A. Save V^T in v. 
     a is overwritten. */
  dgesvd_("A", "A", &n, &n, a, &n, s, u, &n, v, &n, work, &lwork, &info);

  /* Save V in v. */
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
	  tmp = v[i + j*n];
	  v[i + j*n] = v[j + i*n];
	  v[j + i*n] = tmp;
	} 
  }

  del_doubles(work);    
}
