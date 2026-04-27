#include "rkmatrix.h"
#include "basic.h"
#include "blas.h"

#include <stdio.h>
#include <assert.h>

prkmatrix
new_rkmatrix(int k, int rows, int cols) {
  prkmatrix r;
  
  assert(k >= 0);
  assert(rows >= 0);
  assert(cols >= 0);
  
  r = (prkmatrix) malloc(sizeof(rkmatrix));
  if (r == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
  
  r->k = k;
  r->kt = 0;
  r->rows = rows;
  r->cols = cols;
  r->a = allocate_doubles(rows * k);
  r->b = allocate_doubles(cols * k);
  
  return r;
}

prkmatrix
new_random_rkmatrix(int k, int rows, int cols, int seed) {
  prkmatrix r;
  
  r = new_rkmatrix(k, rows, cols);
  r->kt = k;
  fill_random_entries(rows * k, r->a, seed);
  fill_random_entries(cols * k, r->b, seed + 1);
  
  return r;
}

prkmatrix
new_zero_rkmatrix(int k, int rows, int cols) {
  prkmatrix r;
  
  r = new_rkmatrix(k, rows, cols);
  clear_entries(rows * k, r->a);
  clear_entries(cols * k, r->b);
  
  return r;
}

void
del_rkmatrix(prkmatrix r) {
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  
  del_doubles(r->a);
  del_doubles(r->b);
  free(r);
}

void
print_rkmatrix(pcrkmatrix r) { 
  int i, j;

  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  
  printf("%dx%d rkmatrix of rank %d\nA:\n", r->rows, r->cols, r->kt);
  for (i = 0; i < r->rows; i++) {
    for (j = 0; j < r->kt; j++) {
      printf("%.4g ", r->a[j*r->rows + i]);
    }
    printf("\n");
  }
  
  printf("B:\n");
  for (i = 0; i < r->cols; i++) {
    for (j = 0; j < r->kt; j++) {
      printf("%.4g ", r->b[j*r->cols + i]);
    }
    printf("\n");
  }
  printf("\n");
}

size_t
getsize_rkmatrix(int k, int rows, int cols) {
  size_t size;
  
  assert(k >= 0);
  assert(rows >= 0);
  assert(cols >= 0);
  
  size = sizeof(rkmatrix); 
  size += sizeof(double) * (rows + cols) * k;
  
  return size;
}

#ifdef USE_BLAS
void 
addeval_rkmatrix(pcrkmatrix r, const double *v, double *w) {
  double tmp;
  int i;
  
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  assert(r->rows == 0 || w != NULL);
  assert(r->rows > 0 || w == NULL);
  assert(r->cols == 0 || v != NULL);
  assert(r->cols > 0 || v == NULL);
  
  for (i = 0; i < r->kt; i++) {
    /* Compute scaling factor for the i-th column of A. */
    tmp = ddot_(&r->cols, v, eins_, r->b + i*r->cols, eins_);
    /* Scale i-th column of A and add to w. */
    daxpy_(&r->rows, &tmp, r->a + i*r->rows, eins_, w, eins_);
  }
}
#else
void 
addeval_rkmatrix(pcrkmatrix r, const double *v, double *w) {
  int i, j;
  double tmp;

  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  assert(r->rows == 0 || w != NULL);
  assert(r->rows > 0 || w == NULL);
  assert(r->cols == 0 || v != NULL);
  assert(r->cols > 0 || v == NULL);
  
  for (i = 0; i < r->kt; i++) {
    tmp = 0.0;
    for (j = 0; j < r->cols; j++) { 
      tmp += r->b[j + i*r->cols] * v[j];
    }
    for (j = 0; j < r->rows; j++) { 
      w[j] += r->a[j + i*r->rows] * tmp;
    }
  }
}
#endif

void 
eval_rkmatrix(pcrkmatrix r, const double *v, double *w) {
  clear_entries(r->rows, w);
  addeval_rkmatrix(r, v, w);
}

void
convertrk2_fullmatrix(pcrkmatrix r, pfullmatrix f) {  
  assert(f->rows == r->rows); 
  assert(f->cols == r->cols); 
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  assert(f != NULL);
  assert( (f->rows > 0 && f->cols > 0) || f->e == NULL );
  assert(f->rows == 0 || f->cols == 0 || f->e != NULL);
  
  matrixmultrans2(f->rows, f->cols, r->kt, r->a, r->b, f->e);
}

prkmatrix
truncate_rkmatrix(pcrkmatrix r, int k) {
  prkmatrix rtrunc;
  int i, n, m, kr;
  double *qa, *ra, *qb, *rb, *rarbt, *u, *s, *v;

  assert(k >= 0);
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  
  n = r->rows;
  m = r->cols;
  kr = r->kt;
  
  rtrunc = new_rkmatrix(k, n, m);
  if (kr <= k) {
    /* Just copy r into rtrunc. The rest of the entries
       of A and B are not set. */
    rtrunc->kt = kr;
    copy_entries(n * kr, r->a, rtrunc->a);
    copy_entries(m * kr, r->b, rtrunc->b);
  }
  else {
    /* Allocate storage for QR and SVD factors. */
    qa = allocate_doubles(n * kr);
    ra = allocate_doubles(kr * kr);
    qb = allocate_doubles(m * kr);
    rb = allocate_doubles(kr * kr);

    rarbt = allocate_doubles(kr * kr);

    u = allocate_doubles(kr * kr);
    v = allocate_doubles(kr * kr);
    s = allocate_doubles(kr);
      
    /* Compute QR decomposition of A and B. */
    qr_decomposition(r->a, n, kr, qa, ra);
    qr_decomposition(r->b, m, kr, qb, rb);
    
    /* Compute rarbt := R_A * R_B^T and its SVD. */
    matrixmultrans2(kr, kr, kr, ra, rb, rarbt);
    singular_value_decomposition(rarbt, kr, u, s, v);
    
    /* Compute truncation (using compressed SVD). */
    matrixmul(n, k, kr, qa, u, rtrunc->a);
    matrixmul(m, k, kr, qb, v, rtrunc->b);
    if (n < m) {
      for (i = 0; i < k; i++) {
        /* Use i-th singular value s[i] to scale i-th column of A. */
        dscal_(&n, s + i, rtrunc->a + i*n, eins_);
      }
    }
    else {
      for (i = 0; i < k; i++) {
        /* Use i-th singular value s[i] to scale i-th column of B. */
        dscal_(&m, s + i, rtrunc->b + i*m, eins_);
      }
    }
    rtrunc->kt = k;
    
    del_doubles(u);
    del_doubles(s);
    del_doubles(v);
    del_doubles(rarbt);
    del_doubles(qa);
    del_doubles(ra);
    del_doubles(qb);
    del_doubles(rb);
  }
  
  return rtrunc;
}

void
adaptive_truncate_rkmatrix(prkmatrix r, double eps) {
  int n, m, kt, i, k_new;
  double *qa, *ra, *qb, *rb, *rarbt, *u, *s, *v;
  
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  assert(eps >= 0.0);
  
  /* Only try to truncate, if the current rank is greater than 1. */
  if (r->kt > 1) {  
    n = r->rows;
    m = r->cols;
    kt = r->kt;
  
    /* Allocate storage for QR and SVD factors. */
    qa = allocate_doubles(n * kt);
    ra = allocate_doubles(kt * kt);
    qb = allocate_doubles(m * kt);
    rb = allocate_doubles(kt * kt);

    rarbt = allocate_doubles(kt * kt);

    u = allocate_doubles(kt * kt);
    v = allocate_doubles(kt * kt);
    s = allocate_doubles(kt);
      
    /* Compute QR decomposition of A and B. */
    qr_decomposition(r->a, n, kt, qa, ra);
    qr_decomposition(r->b, m, kt, qb, rb);
    
    /* Compute rarbt := R_A * R_B^T and its SVD. */
    matrixmultrans2(kt, kt, kt, ra, rb, rarbt);
    singular_value_decomposition(rarbt, kt, u, s, v);

    /* Compute the new rank k_new such that sigma_(k_new + 1) <= eps * sigma_1. */ 
    k_new = 1;  
    while (k_new < kt && s[k_new] > eps * s[0]) {
      k_new++;
    }
    /* printf("old rank: %d, new rank: %d\n", kt, k_new); */
    
    /* Only truncate, if the current rank is greater than k_new. */
    if (k_new < kt) { 
      /* Compute truncation (using compressed SVD). */
      matrixmul(n, k_new, kt, qa, u, r->a);
      matrixmul(m, k_new, kt, qb, v, r->b);
      if (n < m) {
        for (i = 0; i < k_new; i++) {
          /* Use i-th singular value s[i] to scale i-th column of A. */
          dscal_(&n, s + i, r->a + i*n, eins_);
        }
      }
      else {
        for (i = 0; i < k_new; i++) {
          /* Use i-th singular value s[i] to scale i-th column of B. */
          dscal_(&m, s + i, r->b + i*m, eins_);
        }
      }
      r->kt = k_new;
    }
    
    del_doubles(u);
    del_doubles(s);
    del_doubles(v);
    del_doubles(rarbt);
    del_doubles(qa);
    del_doubles(ra);
    del_doubles(qb);
    del_doubles(rb);
  }
}

double*
getsigma_rkmatrix(pcrkmatrix r) {
  double *atmp, *btmp, *qr_work, *tau1, *tau2, *sigma, *rarbt, *rarbtwork;
  double *u = NULL;
  double *v = NULL;
  int rows, cols, kmax, lwork, info, i, j, k;
  pfullmatrix f;
  
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  
  rows = r->rows;
  cols = r->cols;
  kmax = r->kt;

  sigma = allocate_doubles(min(rows, cols));

  /* There are 3 cases how the singular values are computed. */
  
  /* All singular values are zero. */
  if (kmax == 0) {
    clear_entries(min(rows, cols), sigma);
  }
  /* Compute singular values of F := A * B^T. */
  else if (kmax >= rows || kmax >= cols) {
    /* Convert r to fullmatrix. */
    f = new_fullmatrix(rows, cols);
    convertrk2_fullmatrix(r, f);
    
    lwork = 10 * rows * cols;
    rarbtwork = allocate_doubles(lwork);
    
    /* Compute SVD of fullmatrix (compute only sigma) and check if exit was successful. */
    dgesvd_("N", "N", &rows, &cols, f->e, &rows, sigma, u, &rows, v, &cols,
	        rarbtwork, &lwork, &info);
    if (info != 0) {
      printf("get_sv_rkmatrix: info(dgesvd)=%d\n", info);
      exit(EXIT_FAILURE);
    }
    
    del_fullmatrix(f);
    del_doubles(rarbtwork);
  }
  /* Compute QR factorizations A = Q_A * R_A of A and B = Q_B * R_B
     of B. Compute singular values of rarbt := R_A * R_B^T. */
  else {
    lwork = rows + cols + 10;
    qr_work = allocate_doubles(lwork);
    tau1 = allocate_doubles(min(rows, kmax));
    tau2 = allocate_doubles(min(cols, kmax));
    atmp = allocate_doubles(rows * kmax);
    btmp = allocate_doubles(cols * kmax);
    rarbt = allocate_doubles(kmax * kmax);
  
    /* Copy A and B, so they are not overwritten. */
    copy_entries(rows * kmax, r->a, atmp);
    copy_entries(cols * kmax, r->b, btmp);

    /* Compute QR factorization of A and check if exit was successful. */
    dgeqrf_(&rows, &kmax, atmp, &rows, tau1, qr_work, &lwork, &info);
    if (info != 0) {
      printf("get_sv_rkmatrix: info(dgeqrf,a)=%d\n", info);
      exit(EXIT_FAILURE);
    }
  
    /* Compute QR factorization of B and check if exit was successful. */
    dgeqrf_(&cols, &kmax, btmp, &cols, tau2, qr_work, &lwork, &info);
    if (info != 0) {
      printf("get_sv_rkmatrix: info(dgeqrf,b)=%d\n", info);
      exit(EXIT_FAILURE);
    }
  
    /* Compute rarbt := R_A * R_B^T. The factors R_A and R_B are
       saved in atmp and btmp on and above the diagonal. */
    clear_entries(kmax * kmax, rarbt);
    for (i = 0; i < kmax; i++) {
      for (j = 0; j < kmax; j++) {
        for (k = max(i, j); k < kmax; k++) {
          rarbt[i + kmax*j] += atmp[i + k*rows] * btmp[j + k*cols];
        }
      }
    }
  
    lwork = 10 * kmax * kmax;
    rarbtwork = allocate_doubles(lwork);
  
    /* Compute SVD of rarbt (compute only sigma) and check if exit was successful. */
    dgesvd_("N", "N", &kmax, &kmax, rarbt, &kmax, sigma, u, &kmax, v, 
            &kmax, rarbtwork, &lwork, &info);
    if (info != 0) {
      printf("get_sv_kmatrix: info(dgesvd)=%d\n", info);
      exit(EXIT_FAILURE);
    }

    del_doubles(rarbtwork);
    del_doubles(rarbt);
    del_doubles(btmp);
    del_doubles(atmp);
    del_doubles(tau2);
    del_doubles(tau1);
    del_doubles(qr_work);
  }

  return sigma;
}

double
getsv_rkmatrix(pcrkmatrix r, int nr) {
  double tmp;
  double *sigma;
  
  assert(nr >= 0);
  assert(r != NULL);
  assert(r->k >= r->kt);
  assert(r->kt >= 0);
  assert(r->rows >= 0);
  assert(r->cols >= 0);
  assert( (r->rows > 0 && r->k > 0) || r->a == NULL );
  assert(r->rows == 0 || r->k == 0 || r->a != NULL);
  assert( (r->cols > 0 && r->k > 0) || r->b == NULL );
  assert(r->cols == 0 || r->k == 0 || r->b != NULL);
  
  if (nr > r->rows && nr > r->cols) {
    tmp = 0.0;
  }
  else {
    sigma = getsigma_rkmatrix(r);
    tmp = sigma[nr];
    del_doubles(sigma);
  }
  
  return tmp;
}
