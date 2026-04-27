#ifndef BLAS
#define BLAS

extern int eins_[1];
extern double deins_[1];
extern double dnull_[1];
extern double dmeins_[1];
extern char ttrans[1];
extern char ntrans[1];

/* ============================================================
     BLAS routines
   ============================================================ */

/* Copy x to y. */
extern void
dcopy_(const int *n,		/* Dimension */
       const double *x,		/* Source vector x */
       const int *incx,		/* Increment for x */
       double *y,		    /* Target vector y */
       const int *incy);	/* Increment for y */

/* Compute the Euclidean inner product of vectors x and y. */
extern double
ddot_(const int *n,		    /* Dimension */
      const double *x,		/* Source vector x */
      const int *incx,		/* Increment for x */
      const double *y,		/* Source Vector y */
      const int *incy);		/* Increment for y */

/* Compute x := alpha * x. */
extern void
dscal_(const int *n,		/* Dimension */
       const double *alpha,	/* Scaling factor */
       double *x,		    /* Target vector x */
       const int *incx);	/* Increment for x */

/* Compute y := y + alpha * x. */
extern void
daxpy_(const int *n,		/* Dimension */
       const double *alpha,	/* Scaling factor for x */
       const double *x,		/* Source vector x */
       const int *incx,		/* Increment for x */
       double *y,		    /* Target vector y */
       const int *incy);	/* Increment for y */

/* Compute the matrix-vector product y := beta * y + alpha * A * x. */
extern void
dgemv_(const char *trans,	/* Transpose A? */
       const int *m,		/* Number of rows of A */
       const int *n,		/* Number of columns of A */
       const double *alpha,	/* Scaling factor for A*x */
       const double *a,		/* Matrix A */
       const int *lda,		/* Column length of A (distance in memory between elements of two consecutive columns which have the same row index) */
       const double *x,		/* Vector x */
       const int *incx,		/* Increment for x */
       const double *beta,	/* Scaling factor for y */
       double *y,		    /* Result vector y */
       const int *incy);	/* Increment for y */

/* Compute the matrix-matrix product C := beta * C + alpha * A * B. */
extern void
dgemm_(const char *transa,	/* Transpose A? */
       const char *transb,	/* Transpose B? */
       const int *m,		/* Number of rows of the result C */
       const int *n,		/* Number of columns of the result C */
       const int *k,		/* Number of columns of A and rows of B */
       const double *alpha,	/* Scaling factor for A*B */
       const double *a,		/* Matrix A */
       const int *lda,		/* Column length of A (distance in memory between elements of two consecutive columns which have the same row index) */
       const double *b,		/* Matrix B */
       const int *ldb,		/* Column length of B (distance in memory between elements of two consecutive columns which have the same row index) */
       const double *beta,	/* Scaling factor for C */
       double *c,		    /* Result matrix C */
       const int *ldc);		/* Column length of C (distance in memory between elements of two consecutive columns which have the same row index) */

/* ============================================================
     LAPACK routines
   ============================================================ */

/* Compute the QR factorization A = Q*R. */
extern void
dgeqrf_(const int *m,		/* Number of rows of A */
	const int *n,		    /* Number of columns of A */
	double *a,		        /* Input: m by n matrix, Output: Entries on and above the diagonal are the factor R, entries below the diagonal define elementary reflectors */
	const int *lda,		    /* Column length of A (distance in memory between elements of two consecutive columns which have the same row index) */
	double *tau,		    /* Scalar factors of the elementary reflectors */
	double *work,		    /* Workspace, >= lwork elements */
	const int *lwork,	    /* Size of workspace, >= max(1,n) elements */
	int *info);		        /* Result code */

/* Compute the factor Q of the QR factorization. */
extern void
dorgqr_(const int *m,		/* Number of rows of Q */
	const int *n,		    /* Number of columns of Q, m >= n >= 0 */
	const int *k,		    /* Number of elementary reflectors defining Q, n >= k >= 0*/
	double *A,		        /* Input: Elementary reflectors as returned by DGEQRF, Output: Matrix Q */
	const int *lda,		    /* Column length of A (distance in memory between elements of two consecutive columns which have the same row index) */
	const double *tau,	    /* Scalar factors of the elementary reflectors */
	double *work,		    /* Workspace, >= lwork elements */
	const int *lwork,	    /* Size of workspace, >= max(1,n) */
	int *info);		        /* Result code */

/* Compute the singular value decomposition. */
extern void
dgesvd_(const char *jobu,	/* Is U to be computed? Where to store it? */
	const char *jobvt,	    /* Is V to be computed? Where to store it? */
	const int *m,		    /* Number of rows of A */
	const int *n,		    /* Number of columns of A */
	double *a,		        /* Input: Matrix A */
	const int *lda,		    /* Column length of A (distance in memory between elements of two consecutive columns which have the same row index) */
	double *s,		        /* Singular values of A, in decreasing order */
	double *u,		        /* Output: Matrix U */
	const int *ldu,		    /* Column length of U (distance in memory between elements of two consecutive columns which have the same row index) */
	double *vt,		        /* Output: Matrix V */
	const int *ldvt,	    /* Column length of V (distance in memory between elements of two consecutive columns which have the same row index) */
	double *work,		    /* Workspace, >= lwork */
	const int *lwork,	    /* Size of workspace, >= max(3*min(m,n)+max(m,n), 5*min(m,n)) */
	int *info);		        /* Result code */

/* ============================================================
     wrapper for BLAS and LAPACK routines
   ============================================================ */

/* Compute C := A * B^T. Return rows x cols zero matrix, if mid = 0.
   A is rows x mid matrix and B is cols x mid matrix. */
void
matrixmultrans2(int rows, int cols, int mid,
                const double *a, const double *b, double *c);

/* Compute C := A * B. Return rows x cols zero matrix, if mid = 0.
   A is rows x mid matrix and B is mid x cols matrix. */
void
matrixmul(int rows, int cols, int mid,
	      const double *a, const double *b, double *c); 

/* Compute compressed QR decomposition A = QR of the n times k matrix A
   with n >= k and store the result in the
   n times k matrix Q and k times k matrix R. */
void
qr_decomposition(const double *a, int n, int k, double *q, double *r);

/* Compute SVD A = U*S*V^T of the n times n matrix A
   and store the result in the
   n times n matrix U and
   array s of length n (only diagonal entries of matrix S) and
   n times n matrix V. 
   The array a will be overwritten. */
void
singular_value_decomposition(double *a, int n, double *u, double *s, double *v);

#endif
