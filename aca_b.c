#include "fullmatrix.h"
#include "rkmatrix.h"
#include "aca_b.h"
#include <lapacke.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <openblas/cblas.h>
#include "basic.h"
#include "new_aca.h"

//#include "blas.h"

//floids algorithm 
int* random_unique(int n, int d) {
    if (d > n) return NULL;

    int *S = malloc(d * sizeof(int));
    if (!S) return NULL;

    int size = 0;

    for (int j = n - d; j <= n - 1; j++) {
        int t = rand() % (j + 1);

        int found = 0;
        for (int i = 0; i < size; i++) {
            if (S[i] == t) {
                found = 1;
                break;
            }
        }

        if (!found)
            S[size++] = t;
        else
            S[size++] = j;
    }

    return S;
}

void QR(pfullmatrix W, double eps,
        pfullmatrix *Q_out,
        pfullmatrix *R_out,
        int **Jbar_out,
        int *r_out)
{
    int m = W->rows;
    int n = W->cols;
    int k = (m < n) ? m : n;

    double *A = malloc(sizeof(double) * m * n);
    memcpy(A, W->e, sizeof(double) * m * n);

    double *tau = malloc(sizeof(double) * k);
    int *jpvt = calloc(n, sizeof(int));

    LAPACKE_dgeqp3(LAPACK_COL_MAJOR, m, n, A, m, jpvt, tau);

    /* convert pivots to 0-based */
    printf("jpvt:\n");
    for (int i = 0; i < n; i++){
        jpvt[i]--;
        printf("%d,",jpvt[i]);
    }
    printf("\n");
    pfullmatrix Q = new_fullmatrix(m, m);
    pfullmatrix R = new_fullmatrix(m, n);

    /* extract R */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            if (i <= j)
                R->e[j*m + i] = A[j*m + i];
            else
                R->e[j*m + i] = 0.0;
        }
    }

    /* build Q */
    LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, m, k, A, m, tau);

    memcpy(Q->e, A, sizeof(double) * m * m);

    /* rank */
    double max_diag = fabs(R->e[0*m + 0]);
    int r = 0;

    if (eps > 0) {
        for (int i = 0; i < k; i++) {
            double val = fabs(R->e[i*m + i]);
            if (val > eps * max_diag)
                r++;
        }
    } else {
        r = k;
    }

    *Q_out = Q;
    *R_out = R;
    *Jbar_out = jpvt;
    *r_out = r;

    free(A);
    free(tau);
}

pfullmatrix build_U(pfullmatrix C, int *Jbar, int r)
{
    int m = C->rows;

    pfullmatrix U = new_fullmatrix(m, r);

    for (int j = 0; j < r; j++) {
        int col = Jbar[j];

        for (int i = 0; i < m; i++) {
            U->e[j*U->rows + i] = C->e[col*C->rows + i];
        }
    }

    return U;
}

pfullmatrix build_V(pfullmatrix QtR, int r)
{
    int n = QtR->cols;

    pfullmatrix V = new_fullmatrix(r, n);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < n; j++) {
            V->e[j*V->rows + i] = QtR->e[j*QtR->rows + i];
        }
    }

    return V;
}

void LRID(pfullmatrix C,
          pfullmatrix W,
          pfullmatrix R,
          double eps,
          pfullmatrix *U_out,
          pfullmatrix *V_out,
          int *r_out,
          int **Jbar_out)
{
    pfullmatrix Q, T;
    int *Jbar;
    int r;

    QR(W, eps, &Q, &T, &Jbar, &r);

    /* U = C(:, Jbar) */
    pfullmatrix U = new_fullmatrix(C->rows, r);

    for (int j = 0; j < r; j++) {
        int col = Jbar[j];
        for (int i = 0; i < C->rows; i++)
            U->e[j*C->rows + i] = C->e[col*C->rows + i];
    }

    /* compute Q^T R */
    pfullmatrix QtR = new_fullmatrix(Q->cols, R->cols);

    cblas_dgemm(CblasColMajor,
                CblasTrans, CblasNoTrans,
                Q->cols, R->cols, Q->rows,
                1.0,
                Q->e, Q->rows,
                R->e, R->rows,
                0.0,
                QtR->e, QtR->rows);

    /* solve upper triangular safely */
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,
                   'U','N','N',
                   r, R->cols,
                   T->e, T->rows,
                   QtR->e, QtR->rows);

    pfullmatrix V = new_fullmatrix(r, R->cols);

    for (int j = 0; j < R->cols; j++)
        for (int i = 0; i < r; i++)
            V->e[j*r + i] = QtR->e[j*QtR->rows + i];

    *U_out = U;
    *V_out = V;
    *r_out = r;
    *Jbar_out = Jbar;

    del_fullmatrix(QtR);
    del_fullmatrix(Q);
    del_fullmatrix(T);
}




pfullmatrix 
cholesky(pcfullmatrix A)
{
    int n = A->rows;

    /* must be square */
    if (A->rows != A->cols) {
        printf("Cholesky error: matrix is not square\n");
        return NULL;
    }

    /* allocate result */
    pfullmatrix L = new_fullmatrix(n, n);

    /* copy A (LAPACK overwrites input) */
    memcpy(L->e, A->e, sizeof(double) * n * n);

    /* compute Cholesky: A = L * L^T */
    int info = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR,
        'L',
        n,
        L->e,
        n
    );

    if (info > 0) {
        printf("Cholesky failed: matrix is not positive definite\n");
        free(L->e);
        free(L);
        return NULL;
    }

    /* clean upper triangle */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            L->e[j*n + i] = 0.0;
        }
    }

    return L;
}


double 
frob_T1_T2(pcfullmatrix T1, pcfullmatrix T2)
{
    int n = T1->rows;

    /* temporary matrix */
    pfullmatrix M = new_fullmatrix(n, n);

    /* M = T1 * T2^T */
    cblas_dgemm(
        CblasColMajor,
        CblasNoTrans,   // T1
        CblasTrans,     // T2^T
        n, n, n,
        1.0,
        T1->e, n,
        T2->e, n,
        0.0,
        M->e, n
    );

    /* Frobenius norm */
    double sum = 0.0;
    for (int i = 0; i < n*n; i++)
        sum += M->e[i] * M->e[i];

    free(M->e);
    free(M);

    return sqrt(sum);
}

double 
LRNorm(pcfullmatrix U, pcfullmatrix V)
{
    double A_F;

    pfullmatrix G1 = NULL;
    pfullmatrix G2 = NULL;
    pfullmatrix T1 = NULL;
    pfullmatrix T2 = NULL;

    int m = U->rows;
    int k = U->cols;
    int n = V->cols;

    /* G1 = U^T U (k x k) */
    G1 = new_fullmatrix(k, k);
    cblas_dgemm(
        CblasColMajor,
        CblasTrans,
        CblasNoTrans,
        k, k, m,
        1.0,
        U->e, m,
        U->e, m,
        0.0,
        G1->e, k
    );

    /* G2 = V^T V (k x k) */
    G2 = new_fullmatrix(k, k);
    cblas_dgemm(
        CblasColMajor,
        CblasTrans,
        CblasNoTrans,
        k, k, n,
        1.0,
        V->e, n,
        V->e, n,
        0.0,
        G2->e, k
    );

    /* Cholesky */
    printf("G1 size: %d x %d\n", G1->rows, G1->cols);
    printf("G2 size: %d x %d\n", G2->rows, G2->cols);
    T1 = cholesky(G1);
    T2 = cholesky(G2);

    /* Frobenius norm */
    A_F = frob_T1_T2(T1, T2);

    /* cleanup */
    del_fullmatrix(G1);
    del_fullmatrix(G2);
    del_fullmatrix(T1);
    del_fullmatrix(T2);

    return A_F;
}
double LRnormUp(prkmatrix R,
                pfullmatrix U_bar,
                pfullmatrix V_bar,
                double nu,
                double nu_bar)
{
    int m = R->rows;
    int n = R->cols;
    int r = R->k;
    int r_bar = U_bar->cols;

    double *U = R->a;   // m × r (column-major)
    double *B = R->b;   // n × r (column-major)

    /* U^T U_bar (r × r_bar) */
    double *UtU_bar = calloc(r * r_bar, sizeof(double));

    cblas_dgemm(CblasColMajor,
                CblasTrans, CblasNoTrans,
                r, r_bar, m,
                1.0,
                U, m,
                U_bar->e, m,
                0.0,
                UtU_bar, r);

    /* V_bar * B (r_bar × r) */
    double *VbarB = calloc(r_bar * r, sizeof(double));

    cblas_dgemm(CblasColMajor,
                CblasNoTrans, CblasNoTrans,
                r_bar, r, n,
                1.0,
                V_bar->e, r_bar,
                B, n,
                0.0,
                VbarB, r_bar);

    /* Frobenius inner product */
    double sum = 0.0;

    for (int j = 0; j < r_bar; j++) {
        for (int i = 0; i < r; i++) {
            double Uij = UtU_bar[i + j*r];        // (i,j)
            double Vij = VbarB[i + j*r_bar];      // (i,j) FIXED
            sum += Uij * Vij;
        }
    }

    free(UtU_bar);
    free(VbarB);

    double s = nu*nu + nu_bar*nu_bar + 2.0 * sum;

    return sqrt(fmax(s, 0.0));
}

prkmatrix
b_aca_rkmatrix_new(double eps,int d, pcfullmatrix A){
    int k_max;
    int cols = A->cols; 
    int rows = A->rows;
    prkmatrix r;
    double u=0,v;
    pfullmatrix C,W,R,Q,T,U_k,V_k;
    int d_k,i,j,*piv_rows,*piv_cols, *J_bar;
    
    
    piv_cols=random_unique(cols,d);

    k_max=min(cols,rows);
    r= new_zero_rkmatrix(k_max,rows,cols);
    r->kt=0;
    do{
        // select I_k sceleton rows 
        printf("\nd:%d\n",d);
        printf("\npivot cols:\n");
        C= new_fullmatrix(rows,d);
        for (i=0;i<d;i++){
            printf("%d,",piv_cols[i]);
            for (j=0;j<rows;j++){
                C->e[i*rows + j]=compute_entry_aca_new(r,r->kt,j,piv_cols[i],A);
            }
        }
        printf("\n");
        pfullmatrix C_T=transpose_fullmatrix(C);
        QR(C_T,0,&Q, &T, &piv_rows, &d_k);
        
        // set R TODO: do i need to change d to d_k (side note it should not make a difference since d_k should be d )
        R=new_fullmatrix(d,cols);
        printf("\npivot rows:\n");
        for (i=0;i<d;i++){
            printf("%d",piv_rows[i]);
            for(j=0;j<cols;j++){
              R->e[j*d + i]=compute_entry_aca_new(r, r->kt, piv_rows[i], j, A);
              
            }
        }

        // set W
        W=new_fullmatrix(d,d);
        for (i=0;i<d;i++){
            for (j=0;j<d;j++){
                W->e[j*d+i]=compute_entry_aca_new(r, r->kt,piv_rows[i],piv_cols[j],A);
            }
        }

        //
        /*printf("C:");
        print_fullmatrix(C);
        printf("W");
        print_fullmatrix(W);
        printf("R");
        print_fullmatrix(R);*/
        LRID(C,W,R,eps,&U_k,&V_k,&d_k,&J_bar);
        if (d_k <= 0){
            printf("d_k is zero\n");
            break;
        }
        // update I and J (do i even need to do that? )
        //...
        printf("U_k");
        print_fullmatrix(U_k);
        printf("V_k");
        print_fullmatrix(V_k);
        printf("this is kt:%d\n",r->kt);
        v=LRNorm(U_k,V_k);

        u=LRnormUp(r,U_k,V_k,u,v);
        for(i=0;i<d_k;i++){ 
            for(j=0;j<rows;j++){
                r->a[(i + r->kt)*rows + j] = U_k->e[i*U_k->rows + j];
            }
            for(j=0;j<cols;j++){
                r->b[(i + r->kt)*cols + j] = V_k->e[i*V_k->rows + j];
                } 
        }

        r->kt+=d;//d_k

        // select J_k sceleton columns
        int r_out;
        QR(R,0,&Q,&T,&piv_cols,&r_out); //TODO: &r_out ist noch falsch 
        //r->kt+=d_k;//update k 

    }while(v>=eps*u && (r->kt+d)<k_max);

    return r;
}




