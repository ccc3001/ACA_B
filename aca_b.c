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
#include "aca.h"
#include "svd.h"
#include <omp.h>
#include <cjson/cJSON.h>
ACAResidualNode *new_residual_node(void)
{
    ACAResidualNode *n =
        malloc(sizeof(ACAResidualNode));

    memset(n, 0, sizeof(ACAResidualNode));

    return n;
}

void del_residual_node(ACAResidualNode *node)
{
    if (!node)
        return;

    for (int i = 0; i < 4; i++)
        del_residual_node(node->child[i]);

    free(node->u);
    free(node->v);
    free(node->rank_inc);

    free(node);
}
cJSON *aca_residual_node_to_json(
    ACAResidualNode *node)
{
    if (!node)
        return cJSON_CreateNull();

    cJSON *json = cJSON_CreateObject();

    cJSON_AddNumberToObject(json, "row_start", node->row_start);
    cJSON_AddNumberToObject(json, "row_size", node->row_size);
    cJSON_AddNumberToObject(json, "col_start", node->col_start);
    cJSON_AddNumberToObject(json, "col_size", node->col_size);
    cJSON_AddNumberToObject(json, "level", node->level);
    cJSON_AddNumberToObject(json, "time", node->time);

    if (node->u)
    {
        cJSON_AddItemToObject(
            json,
            "u",
            cJSON_CreateDoubleArray(
                node->u,
                node->rank_len));
    }

    if (node->v)
    {
        cJSON_AddItemToObject(
            json,
            "v",
            cJSON_CreateDoubleArray(
                node->v,
                node->rank_len));
    }

    if (node->rank_inc)
    {
        cJSON_AddItemToObject(
            json,
            "rank_inc",
            cJSON_CreateDoubleArray(
                node->rank_inc,
                node->rank_len));
    }

    cJSON *children = cJSON_CreateArray();

    for (int i = 0; i < 4; i++)
    {
        cJSON_AddItemToArray(
            children,
            aca_residual_node_to_json(
                node->child[i]));
    }

    cJSON_AddItemToObject(json, "children", children);

    return json;
}
// floids algorithm
int *random_unique(int n, int d)
{
    if (d > n)
        return NULL;

    int *S = malloc(d * sizeof(int));
    if (!S)
        return NULL;

    int size = 0;

    for (int j = n - d; j <= n - 1; j++)
    {
        int t = rand() % (j + 1);

        int found = 0;
        for (int i = 0; i < size; i++)
        {
            if (S[i] == t)
            {
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

void QR(pfullmatrix W,
        double eps, double *max_pivot,
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

    /* pivoted QR */
    LAPACKE_dgeqp3(
        LAPACK_COL_MAJOR,
        m,
        n,
        A,
        m,
        jpvt,
        tau);

    /* convert pivots to 0-based */
    for (int i = 0; i < n; i++)
        jpvt[i]--;

    /* build R */
    pfullmatrix R = new_fullmatrix(m, n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {

            if (i <= j)
                R->e[j * m + i] = A[j * m + i];
            else
                R->e[j * m + i] = 0.0;
        }
    }

    /* build Q */
    pfullmatrix Q = new_fullmatrix(m, m);

    LAPACKE_dorgqr(
        LAPACK_COL_MAJOR,
        m,
        m,
        k,
        A,
        m,
        tau);

    memcpy(Q->e, A, sizeof(double) * m * m);

    /* stable numerical rank detection */
    double pivot_abs = 1e-50;

    double first_pivot = fabs(R->e[0]);

    if (first_pivot > (*max_pivot))
    {
        *max_pivot = first_pivot;
    }

    //*max_pivot = first_pivot;
    int r = 0;

    for (int i = 0; i < k; i++)
    {

        double piv =
            fabs(R->e[i * m + i]);

        /* absolute threshold */
        if (piv < pivot_abs)
            break;

        /* relative threshold */
        // if (piv < eps * first_pivot)
        if (piv < eps * (*max_pivot))
            break;
        r++;
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

    for (int j = 0; j < r; j++)
    {
        int col = Jbar[j];

        for (int i = 0; i < m; i++)
        {
            U->e[j * U->rows + i] = C->e[col * C->rows + i];
        }
    }

    return U;
}

pfullmatrix build_V(pfullmatrix QtR, int r)
{
    int n = QtR->cols;

    pfullmatrix V = new_fullmatrix(r, n);

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V->e[j * V->rows + i] = QtR->e[j * QtR->rows + i];
        }
    }

    return V;
}
void LRID(pfullmatrix C,
          pfullmatrix W,
          pfullmatrix R,
          double eps,
          double *max_piv,
          pfullmatrix *U_out,
          pfullmatrix *V_out,
          int *r_out,
          int **Jbar_out)
{
    pfullmatrix Q;
    pfullmatrix T;

    int *Jbar;
    int r;

    printf("old max pivot changed to:%.16f\n", *max_piv);
    QR(W, eps, max_piv, &Q, &T, &Jbar, &r);

    if (r <= 0)
    {
        *r_out = 0;
        printf("LRID: QR is of rank %d \n", r);
        return;
    }

    /* ---------------------------------
       U = C(:,Jbar(1:r))
       --------------------------------- */

    pfullmatrix U =
        new_fullmatrix(C->rows, r);

    for (int j = 0; j < r; j++)
    {

        int col = Jbar[j];

        for (int i = 0; i < C->rows; i++)
        {

            U->e[j * C->rows + i] =
                C->e[col * C->rows + i];
        }
    }

    /* ---------------------------------
       build triangular T11
       --------------------------------- */

    pfullmatrix T11 =
        new_fullmatrix(r, r);

    for (int j = 0; j < r; j++)
    {

        for (int i = 0; i <= j; i++)
        {

            T11->e[j * r + i] =
                T->e[j * T->rows + i];
        }
    }

    /* ---------------------------------
       QtR = Q_r^T R
       only first r columns of Q
       --------------------------------- */

    pfullmatrix QtR =
        new_fullmatrix(r, R->cols);

    cblas_dgemm(
        CblasColMajor,
        CblasTrans,
        CblasNoTrans,
        r,
        R->cols,
        Q->rows,
        1.0,
        Q->e,
        Q->rows,
        R->e,
        R->rows,
        0.0,
        QtR->e,
        QtR->rows);

    /* ---------------------------------
       solve T11 * X = QtR
       --------------------------------- */

    int info = LAPACKE_dtrtrs(
        LAPACK_COL_MAJOR,
        'U',
        'N',
        'N',
        r,
        R->cols,
        T11->e,
        r,
        QtR->e,
        QtR->rows);

    if (info != 0)
    {

        printf("dtrtrs failed: %d\n", info);

        del_fullmatrix(U);
        del_fullmatrix(QtR);
        del_fullmatrix(Q);
        del_fullmatrix(T);
        del_fullmatrix(T11);

        *r_out = 0;
        return;
    }

    /* ---------------------------------
       V
       --------------------------------- */

    pfullmatrix V =
        new_fullmatrix(r, R->cols);

    for (int j = 0; j < R->cols; j++)
    {

        for (int i = 0; i < r; i++)
        {

            V->e[j * r + i] =
                QtR->e[j * QtR->rows + i];
        }
    }

    *U_out = U;
    *V_out = V;
    *r_out = r;
    printf("LRID: rank = %d U: %d x %d, V: %d x %d\n", r, U->rows, U->cols, V->rows, V->cols);
    *Jbar_out = Jbar;

    del_fullmatrix(QtR);
    del_fullmatrix(Q);
    del_fullmatrix(T);
    del_fullmatrix(T11);
}

pfullmatrix
cholesky(pcfullmatrix A)
{
    int n = A->rows;

    /* must be square */
    if (A->rows != A->cols)
    {
        printf("Cholesky error: matrix is not square\n");
        return NULL;
    }

    /* allocate result */
    pfullmatrix L = new_fullmatrix(n, n);

    /* copy A (LAPACK overwrites input) */
    memcpy(L->e, A->e, sizeof(double) * n * n);
    /*double reg = 1e-14;

    for (int i = 0; i < n; i++)
    {
        L->e[i * n + i] += reg;
    }*/

    /* compute Cholesky: A = L * L^T */
    int info = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR,
        'L',
        n,
        L->e,
        n);

    if (info > 0)
    {
        printf("\n=========================================\n");
        printf("Cholesky failed.\n");
        printf("Leading principal minor of order %d is not positive definite.\n", info);

        int k = info - 1;

        printf("\nPivot diagonal A(%d,%d) = %.16e\n",
               k, k, L->e[k * n + k]);

        /* symmetry check */
        double max_asym = 0.0;
        int imax = 0, jmax = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double diff = fabs(L->e[i + j * n] - L->e[j + i * n]);

                if (diff > max_asym)
                {
                    max_asym = diff;
                    imax = i;
                    jmax = j;
                }
            }
        }

        printf("\nMaximum asymmetry = %.16e at (%d,%d)\n",
               max_asym, imax, jmax);

        /* diagonal statistics */
        double mind = L->e[0];
        double maxd = L->e[0];

        for (int i = 1; i < n; i++)
        {
            double d = L->e[i * n + i];

            if (d < mind)
                mind = d;

            if (d > maxd)
                maxd = d;
        }

        printf("Diagonal range: [% .16e , % .16e]\n",
               mind, maxd);

        /* print failing row */
        printf("\nFailing row %d:\n", k);

        for (int j = 0; j < n; j++)
            printf("% .3e ", L->e[k + j * n]);

        printf("\n");

        printf("=========================================\n");

        free(L->e);
        free(L);

        return NULL;
    }

    /* clean upper triangle */
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            L->e[j * n + i] = 0.0;
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
        CblasNoTrans, // T1
        CblasTrans,   // T2^T
        n, n, n,
        1.0,
        T1->e, n,
        T2->e, n,
        0.0,
        M->e, n);

    /* Frobenius norm */
    double sum = 0.0;
    for (int i = 0; i < n * n; i++)
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
        G1->e, k);

    /* G2 = V^T V (k x k) */
    G2 = new_fullmatrix(k, k);
    cblas_dgemm(
        CblasColMajor,
        CblasNoTrans,
        CblasTrans,
        k, /* rows */
        k, /* cols */
        n, /* inner dimension */
        1.0,
        V->e,
        k, /* lda = number of rows = k */
        V->e,
        k,
        0.0,
        G2->e,
        k);

    /* Cholesky */
    // printf("G1 size: %d x %d\n", G1->rows, G1->cols);
    // printf("G2 size: %d x %d\n", G2->rows, G2->cols);
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
    int r = R->kt;
    int r_bar = U_bar->cols;

    double *U = R->a; /* m x r */
    double *B = R->b; /* n x r */

    /* U^T U_bar : r x r_bar */
    double *UtU_bar = calloc(r * r_bar, sizeof(double));

    cblas_dgemm(CblasColMajor,
                CblasTrans,
                CblasNoTrans,
                r,
                r_bar,
                m,
                1.0,
                U,
                m,
                U_bar->e,
                m,
                0.0,
                UtU_bar,
                r);

    /* V * V_bar^T : r x r_bar */
    double *VVbarT = calloc(r * r_bar, sizeof(double));

    cblas_dgemm(CblasColMajor,
                CblasTrans,
                CblasTrans,
                r,
                r_bar,
                n,
                1.0,
                B,
                n,
                V_bar->e,
                r_bar,
                0.0,
                VVbarT,
                r);

    double sum = 0.0;

    for (int j = 0; j < r_bar; j++)
    {
        for (int i = 0; i < r; i++)
        {
            sum += UtU_bar[i + j * r] *
                   VVbarT[i + j * r];
        }
    }

    free(UtU_bar);
    free(VVbarT);

    double s = nu * nu + nu_bar * nu_bar + 2.0 * sum;

    return sqrt(fmax(s, 0.0));
}

prkmatrix
b_aca_rkmatrix_new(double eps, int d, pcfullmatrix A, double **residuals_u, double **residuals_v, double **rank_increase)
{
    int rows = A->rows;
    int cols = A->cols;
    int k_max = min(rows, cols);
    *residuals_u = allocate_doubles(k_max);
    *residuals_v = allocate_doubles(k_max);
    *rank_increase = allocate_doubles(k_max);
    int i, j;
    double max_piv = 0.0;
    double u = 0.0;
    double v = 0.0;
    double piv;
    prkmatrix r = NULL;

    int *piv_cols = NULL;

    /* ---------------------------------
       initial pivot columns
       --------------------------------- */

    piv_cols = random_unique(cols, d);

    if (!piv_cols)
        return NULL;

    /* ---------------------------------
       allocate rk matrix
       --------------------------------- */

    r = new_zero_rkmatrix(k_max, rows, cols);

    if (!r)
    {
        free(piv_cols);
        return NULL;
    }

    r->kt = 0;

    /* =================================
       ACA iteration
       ================================= */
    int iter = 0;
    do
    {
        iter += 1;
        pfullmatrix C = NULL;
        pfullmatrix C_T = NULL;
        pfullmatrix R = NULL;
        pfullmatrix W = NULL;

        pfullmatrix Q = NULL;
        pfullmatrix T = NULL;

        pfullmatrix U_k = NULL;
        pfullmatrix V_k = NULL;

        int *piv_rows = NULL;
        int *J_bar = NULL;

        int d_k = 0;

        /* ---------------------------------
           build C
           --------------------------------- */

        C = new_fullmatrix(rows, d);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < rows; j++)
            {
                C->e[i * rows + j] =
                    compute_entry_aca_new(
                        r,
                        r->kt,
                        j,
                        piv_cols[i],
                        A);
            }
        }

        /* ---------------------------------
           QR(C^T)
           --------------------------------- */

        C_T = transpose_fullmatrix(C);
        piv = 0.0;
        QR(C_T, 0.0, &piv,
           &Q,
           &T,
           &piv_rows,
           &d_k);

        if (d_k <= 0)
        {
            printf("QR rank failure d_k = %d\n", d_k);

            del_fullmatrix(C);
            del_fullmatrix(C_T);
            del_fullmatrix(Q);
            del_fullmatrix(T);

            free(piv_rows);

            break;
        }

        /* ---------------------------------
           build R
           --------------------------------- */

        R = new_fullmatrix(d, cols);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < cols; j++)
            {
                R->e[j * d + i] =
                    compute_entry_aca_new(
                        r,
                        r->kt,
                        piv_rows[i],
                        j,
                        A);
            }
        }

        /* ---------------------------------
           build W
           --------------------------------- */

        W = new_fullmatrix(d, d);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < d; j++)
            {
                W->e[j * d + i] =
                    compute_entry_aca_new(
                        r,
                        r->kt,
                        piv_rows[i],
                        piv_cols[j],
                        A);
            }
        }

        /* ---------------------------------
           LRID
           --------------------------------- */

        LRID(C,
             W,
             R,
             eps,
             &max_piv,
             &U_k,
             &V_k,
             &d_k,
             &J_bar);

        if (d_k <= 0 || !U_k || !V_k)
        {
            printf("LRID failed\n");

            del_fullmatrix(C);
            del_fullmatrix(C_T);
            del_fullmatrix(R);
            del_fullmatrix(W);

            del_fullmatrix(Q);
            del_fullmatrix(T);

            free(piv_rows);
            free(J_bar);

            break;
        }

        /* ---------------------------------
           norm update
           --------------------------------- */

        v = LRNorm(U_k, V_k);

        u = LRnormUp(r,
                     U_k,
                     V_k,
                     u,
                     v);

        /* ---------------------------------
           append low-rank block
           --------------------------------- */

        for (i = 0; i < d_k; i++)
        {
            for (j = 0; j < rows; j++)
            {
                r->a[(i + r->kt) * rows + j] =
                    U_k->e[i * U_k->rows + j];
            }

            for (j = 0; j < cols; j++)
            {
                r->b[(i + r->kt) * cols + j] =
                    V_k->e[j * V_k->rows + i];
            }
        }
        (*residuals_u)[iter] = u;
        (*residuals_v)[iter] = v;
        (*rank_increase)[iter] = d_k;
        r->kt += d_k;

        /* ---------------------------------
           update pivot columns
           IMPORTANT:
           free old piv_cols first
           --------------------------------- */

        free(piv_cols);
        piv_cols = NULL;

        del_fullmatrix(Q);
        del_fullmatrix(T);

        Q = NULL;
        T = NULL;

        {
            int r_out;
            piv = 0;
            QR(R,
               0.0,
               &piv,
               &Q,
               &T,
               &piv_cols,
               &r_out);
        }

        /* ---------------------------------
           cleanup iteration
           --------------------------------- */

        del_fullmatrix(C);
        del_fullmatrix(C_T);

        del_fullmatrix(R);
        del_fullmatrix(W);

        del_fullmatrix(Q);
        del_fullmatrix(T);

        del_fullmatrix(U_k);
        del_fullmatrix(V_k);

        free(piv_rows);
        free(J_bar);

    } while (v >= eps * u &&
             (r->kt + d) < k_max);

    free(piv_cols);

    return r;
}

void h_b_aca_rkmatrix_new(double eps, int d, int L, pcfullmatrix A, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V, int *r)
{
    int rows = A->rows;
    int cols = A->cols;
    // int k_max = min(rows, cols);

    if (L == 0)
    {
        double *residuals_u, *residuals_v, *rank_increase;
        prkmatrix r_ = b_aca_rkmatrix_new(eps, d, A, &residuals_u, &residuals_v, &rank_increase);
        free(residuals_v);
        free(residuals_u);
        free(rank_increase);
        *U = new_fullmatrix(rows, r_->kt);
        pfullmatrix V_ = new_fullmatrix(cols, r_->kt);
        memcpy((*U)->e, r_->a, rows * r_->kt * sizeof(double));
        memcpy(V_->e, r_->b, r_->kt * cols * sizeof(double));
        *V = transpose_fullmatrix(V_);

        (*S) = new_zero_fullmatrix(r_->kt, r_->kt);
        for (int i = 0; i < r_->kt; i++)
            (*S)->e[i * (*S)->rows + i] = 1.0;
        (*r) = r_->kt;
        return;
    }
    if (L > 0)
    {
        int row1 = rows / 2;
        int row2 = rows - row1;
        int col1 = cols / 2;
        int col2 = cols - col1;
        pfullmatrix A11 = new_submatrix(A, 0, 0, row1, col1);
        pfullmatrix A12 = new_submatrix(A, 0, col1, row1, col2);
        pfullmatrix A21 = new_submatrix(A, row1, 0, row2, col1);
        pfullmatrix A22 = new_submatrix(A, row1, col1, row2, col2);
        pfullmatrix U11, S11, V11;
        pfullmatrix U12, S12, V12;
        pfullmatrix U21, S21, V21;
        pfullmatrix U22, S22, V22;
        int r11, r12, r21, r22;
#pragma omp parallel sections num_threads(4)
        {
#pragma omp section
            {
                h_b_aca_rkmatrix_new(eps, d, L - 1, A11, &U11, &S11, &V11, &r11);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix_new(eps, d, L - 1, A12, &U12, &S12, &V12, &r12);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix_new(eps, d, L - 1, A21, &U21, &S21, &V21, &r21);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix_new(eps, d, L - 1, A22, &U22, &S22, &V22, &r22);
            }
        }
        // U_1=[U11*S11, U12*S12]
        pfullmatrix U_11 = new_fullmatrix(U11->rows, S11->cols);
        pfullmatrix U_12 = new_fullmatrix(U12->rows, S12->cols);
        mul_fullmatrix(U11, S11, U_11);
        mul_fullmatrix(U12, S12, U_12);
        if (U_11->rows != U_12->rows)
        {
            printf("Error: U_11 and U_12 have different number of rows\n");
            return;
        }

        pfullmatrix U_1 = new_fullmatrix(U_11->rows, U_11->cols + U_12->cols);
        for (int j = 0; j < U_11->cols; j++)
        {
            for (int i = 0; i < U_11->rows; i++)
            {
                U_1->e[j * U_1->rows + i] = U_11->e[j * U_11->rows + i];
            }
        }
        for (int j = 0; j < U_12->cols; j++)
        {
            for (int i = 0; i < U_12->rows; i++)
            {
                U_1->e[(j + U_11->cols) * U_1->rows + i] = U_12->e[j * U_12->rows + i];
            }
        }

        // SVD(U_1,eps,&U1,&S1,&V1);
        pfullmatrix U1, S1, V1;
        SVD_truncated(U_1, eps, &U1, &S1, &V1);

        // diag(V11,V12)
        pfullmatrix V_1 = new_zero_fullmatrix(V11->rows + V12->rows, V11->cols + V12->cols);

        for (int j = 0; j < V11->cols; j++)
        {
            for (int i = 0; i < V11->rows; i++)
            {
                V_1->e[j * V_1->rows + i] = V11->e[j * V11->rows + i];
            }
        }
        for (int j = 0; j < V12->cols; j++)
        {
            for (int i = 0; i < V12->rows; i++)
            {
                V_1->e[(j + V11->cols) * V_1->rows + (i + V11->rows)] = V12->e[j * V12->rows + i];
            }
        }
        // pfullmatrix V1 =V1*V_1;
        pfullmatrix V1_temp = new_fullmatrix(V1->rows, V_1->cols);
        mul_fullmatrix(V1, V_1, V1_temp);
        // U_2=[U21*S21, U22*S22]
        pfullmatrix U_21 = new_fullmatrix(U21->rows, S21->cols);
        pfullmatrix U_22 = new_fullmatrix(U22->rows, S22->cols);
        mul_fullmatrix(U21, S21, U_21);
        mul_fullmatrix(U22, S22, U_22);
        if (U_21->rows != U_22->rows)
        {
            printf("Error: U_21 and U_22 have different number of rows\n");
            return;
        }
        pfullmatrix U_2 = new_fullmatrix(U_21->rows, U_21->cols + U_22->cols);
        for (int j = 0; j < U_21->cols; j++)
        {
            for (int i = 0; i < U_21->rows; i++)
            {
                U_2->e[j * U_2->rows + i] = U_21->e[j * U_21->rows + i];
            }
        }
        for (int j = 0; j < U_22->cols; j++)
        {
            for (int i = 0; i < U_22->rows; i++)
            {
                U_2->e[(j + U_21->cols) * U_2->rows + i] = U_22->e[j * U_22->rows + i];
            }
        }
        // SVD(U_2,eps,&U2,&S2,&V2);

        pfullmatrix U2, S2, V2;
        SVD_truncated(U_2, eps, &U2, &S2, &V2);

        // diag(V21,V22)
        pfullmatrix V_2 = new_zero_fullmatrix(V21->rows + V22->rows, V21->cols + V22->cols);
        for (int j = 0; j < V21->cols; j++)
        {
            for (int i = 0; i < V21->rows; i++)
            {
                V_2->e[j * V_2->rows + i] = V21->e[j * V21->rows + i];
            }
        }
        for (int j = 0; j < V22->cols; j++)
        {
            for (int i = 0; i < V22->rows; i++)
            {
                V_2->e[(j + V21->cols) * V_2->rows + (i + V21->rows)] = V22->e[j * V22->rows + i];
            }
        }
        // TODO: pfullmatrix V2 =V2*V_2;
        pfullmatrix V2_temp = new_fullmatrix(V2->rows, V_2->cols);
        mul_fullmatrix(V2, V_2, V2_temp);

        /* --------------------------------------------------------------------*/
        /*                       after for loop                                */
        /* --------------------------------------------------------------------*/

        // U_ =diag(U1,U2);
        pfullmatrix U_ = new_zero_fullmatrix(U1->rows + U2->rows, U1->cols + U2->cols);
        for (int j = 0; j < U1->cols; j++)
        {
            for (int i = 0; i < U1->rows; i++)
            {
                U_->e[j * U_->rows + i] = U1->e[j * U1->rows + i];
            }
        }
        for (int j = 0; j < U2->cols; j++)
        {
            for (int i = 0; i < U2->rows; i++)
            {
                U_->e[(j + U1->cols) * U_->rows + (i + U1->rows)] = U2->e[j * U2->rows + i];
            }
        }

        // V_=[S1*V1; S2*V2]

        pfullmatrix S1V1 = new_fullmatrix(S1->rows, V1_temp->cols);

        mul_fullmatrix(S1, V1_temp, S1V1);
        pfullmatrix S2V2 = new_fullmatrix(S2->rows, V2_temp->cols);

        mul_fullmatrix(S2, V2_temp, S2V2);

        pfullmatrix V_ = new_zero_fullmatrix(S1V1->rows + S2V2->rows, S1V1->cols);
        for (int j = 0; j < S1V1->cols; j++)
        {
            for (int i = 0; i < S1V1->rows; i++)
            {
                V_->e[j * V_->rows + i] =
                    S1V1->e[j * S1V1->rows + i];
            }
        }

        for (int j = 0; j < S2V2->cols; j++)
        {
            for (int i = 0; i < S2V2->rows; i++)
            {
                V_->e[j * V_->rows +
                      (i + S1V1->rows)] =
                    S2V2->e[j * S2V2->rows + i];
            }
        }
        // SVD(V_,eps,&U,&S,&V);
        pfullmatrix Uloc;
        pfullmatrix Sloc;
        pfullmatrix Vloc;
        SVD_truncated(V_, eps, &Uloc, &Sloc, &Vloc);
        *U = Uloc;
        *S = Sloc;
        *V = Vloc;
        (*r) = (*S)->rows;

        // *U=U*U_;

        pfullmatrix U_temp = new_fullmatrix((*U)->rows, U_->cols);

        mul_fullmatrix(U_, Uloc, U_temp);

        *U = U_temp;
        return;
    }
}

double *new_subnodes(
    const double *nodes,
    int d,
    int n,
    int start,
    int size)
{
    double *sub = allocate_doubles(d * size);

    for (int m = 0; m < d; m++)
    {
        for (int i = 0; i < size; i++)
        {
            sub[m * size + i] = nodes[m * n + start + i];
        }
    }

    return sub;
}
prkmatrix
b_aca_rkmatrix(double eps, int d, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y,
               double (*test_function)(int dim, const double *x, const double *y), double **residuals_u, double **residuals_v, double **rank_increase)
{
    int retry_count = 0;
    const int max_retries = 2;
    assert(nodes_x != NULL || d == 0 || rows == 0);
    assert(nodes_y != NULL || d == 0 || cols == 0);
    int k_max = min(rows, cols);
    *residuals_u = allocate_doubles(k_max);
    *residuals_v = allocate_doubles(k_max);
    *rank_increase = allocate_doubles(k_max);

    memset(*residuals_u, 0, k_max * sizeof(double));
    memset(*residuals_v, 0, k_max * sizeof(double));
    memset(*rank_increase, 0, k_max * sizeof(double));

    int i, j;
    double *x, *y;
    double max_piv = 0.0;
    double u = 0.0;
    double v = 0.0;
    double piv;

    prkmatrix r = NULL;

    int *piv_cols = NULL;
    int *prev_piv_cols = NULL;

    /* ---------------------------------
       initial pivot columns
       --------------------------------- */

    piv_cols = random_unique(cols, d);
    prev_piv_cols = malloc(d * sizeof(int));
    memcpy(prev_piv_cols, piv_cols, d * sizeof(int));
    /*printf("initial piv_cols: ");
    for (int i = 0; i < d; i++)
    {
        printf("%d ", piv_cols[i]);
    }
    printf("\n");*/
    if (!piv_cols)
        return NULL;

    /* ---------------------------------
       allocate rk matrix
       --------------------------------- */

    r = new_zero_rkmatrix(k_max, rows, cols);
    x = allocate_doubles(dim * rows);
    y = allocate_doubles(dim * cols);

    for (j = 0; j < dim; j++)
    {
        for (i = 0; i < rows; i++)
        {
            x[dim * i + j] = nodes_x[j * rows + i];
        }
    }
    for (j = 0; j < dim; j++)
    {
        for (i = 0; i < cols; i++)
        {
            y[dim * i + j] = nodes_y[j * cols + i];
        }
    }

    if (!r)
    {
        free(piv_cols);
        return NULL;
    }

    r->kt = 0;

    /* =================================
       ACA iteration
       ================================= */
    int iter = 0;
    do
    {

        pfullmatrix C = NULL;
        pfullmatrix C_T = NULL;
        pfullmatrix R = NULL;
        pfullmatrix W = NULL;

        pfullmatrix Q = NULL;
        pfullmatrix T = NULL;

        pfullmatrix U_k = NULL;
        pfullmatrix V_k = NULL;

        int *piv_rows = NULL;
        int *J_bar = NULL;

        int d_k = 0;

        /* ---------------------------------
           build C
           --------------------------------- */

        C = new_fullmatrix(rows, d);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < rows; j++)
            {
                C->e[i * rows + j] = compute_entry_aca(r, dim, r->kt, j, piv_cols[i], x + dim * j, y + dim * piv_cols[i], test_function);
            }
        }

        /* ---------------------------------
           QR(C^T)
           --------------------------------- */

        C_T = transpose_fullmatrix(C);
        piv = 0;
        QR(C_T, 0.0, &piv,
           &Q,
           &T,
           &piv_rows,
           &d_k);
        if (d_k <= 0)
        {
            printf("QR rank failure d_k = %d\n", d_k);

            del_fullmatrix(C);
            del_fullmatrix(C_T);
            del_fullmatrix(Q);
            del_fullmatrix(T);

            free(piv_rows);

            free(piv_cols);
            piv_cols = random_unique(cols, d);

            retry_count++;

            if (retry_count >= max_retries)
            {
                printf("Too many QR failures\n");
                break;
            }

            continue;
        }

        /* ---------------------------------
           build R
           --------------------------------- */

        R = new_fullmatrix(d, cols);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < cols; j++)
            {
                R->e[j * d + i] = compute_entry_aca(r, dim, r->kt, piv_rows[i], j, x + dim * piv_rows[i], y + dim * j, test_function);
            }
        }

        /* ---------------------------------
           build W
           --------------------------------- */

        W = new_fullmatrix(d, d);

        for (i = 0; i < d; i++)
        {
            for (j = 0; j < d; j++)
            {
                W->e[j * d + i] = compute_entry_aca(r, dim, r->kt, piv_rows[i], piv_cols[j], x + dim * piv_rows[i], y + dim * piv_cols[j], test_function);
            }
        }

        /* ---------------------------------
           LRID
           --------------------------------- */

        LRID(C,
             W,
             R,
             eps,
             &max_piv,
             &U_k,
             &V_k,
             &d_k,
             &J_bar);
        if (r->kt + d_k > k_max)
        {
            printf("LRID rank exceeds maximum rank: r->kt + d_k = %d > k_max = %d\n", r->kt + d_k, k_max);

            del_fullmatrix(C);
            del_fullmatrix(C_T);
            del_fullmatrix(R);
            del_fullmatrix(W);

            del_fullmatrix(Q);
            del_fullmatrix(T);

            free(piv_rows);
            free(J_bar);

            break;
        }

        if (d_k <= 0 || !U_k || !V_k)
        {
            printf("LRID failed\n");

            del_fullmatrix(C);
            del_fullmatrix(C_T);
            del_fullmatrix(R);
            del_fullmatrix(W);

            del_fullmatrix(Q);
            del_fullmatrix(T);

            free(piv_rows);
            free(J_bar);

            return r;
        }

        /* ---------------------------------
           norm update
           --------------------------------- */

        v = LRNorm(U_k, V_k);

        u = LRnormUp(r,
                     U_k,
                     V_k,
                     u,
                     v);

        /* ---------------------------------
           append low-rank block
           --------------------------------- */

        for (i = 0; i < d_k; i++)
        {
            for (j = 0; j < rows; j++)
            {
                r->a[(i + r->kt) * rows + j] =
                    U_k->e[i * U_k->rows + j];
            }

            for (j = 0; j < cols; j++)
            {
                r->b[(i + r->kt) * cols + j] =
                    V_k->e[j * V_k->rows + i];
            }
        }
        (*residuals_u)[iter] = u;
        (*residuals_v)[iter] = v;
        (*rank_increase)[iter] = d_k;
        iter += 1;
        r->kt += d_k;
        printf("iter=%d d=%d d_k=%d r->kt=%d u=%g v=%g\n", iter, d, d_k, r->kt, u, v);

        /* ---------------------------------
           update pivot columns
           IMPORTANT:
           free old piv_cols first
           --------------------------------- */

        free(piv_cols);
        piv_cols = NULL;

        del_fullmatrix(Q);
        del_fullmatrix(T);

        Q = NULL;
        T = NULL;

        {
            int r_out;
            piv = 0;
            QR(R,
               0.0, &piv,
               &Q,
               &T,
               &piv_cols,
               &r_out);

            /* Compare against previous pivot block */

            int same = 1;

            for (i = 0; i < d; i++)
            {
                if (piv_cols[i] != prev_piv_cols[i])
                {
                    same = 0;
                    break;
                }
            }

            if (same)
            {
                printf("Same pivot block as previous iteration. Choosing random pivots.\n");

                free(piv_cols);
                piv_cols = random_unique(cols, d);
            }

            /* Remember current pivots for next iteration */

            memcpy(prev_piv_cols,
                   piv_cols,
                   d * sizeof(int));

            /*printf("Updated pivot columns: ");
            for (i = 0; i < d; i++)
                printf("%d ", piv_cols[i]);
            printf("\n");*/
        }
        /* ---------------------------------
           cleanup iteration
           --------------------------------- */

        del_fullmatrix(C);
        del_fullmatrix(C_T);

        del_fullmatrix(R);
        del_fullmatrix(W);

        del_fullmatrix(Q);
        del_fullmatrix(T);

        del_fullmatrix(U_k);
        del_fullmatrix(V_k);

        free(piv_rows);
        free(J_bar);
    } while (v >= eps * u &&
             (r->kt + d) < k_max);
    free(prev_piv_cols);
    free(piv_cols);

    return r;
}

void h_b_aca_rkmatrix(double eps, int d, int L, ACAResidualNode *node, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V,
                      int dim, int rows, int cols, const double *nodes_x, const double *nodes_y,
                      double (*test_function)(int dim, const double *x, const double *y), int *r)
{
    node->row_size = rows;
    node->col_size = cols;
    node->level = L;

    if (L == 0)
    {
        double start = omp_get_wtime();

        double *u_res;
        double *v_res;
        double *rank_inc;
        printf(
            "LEAF L=%d rows=%d cols=%d "
            "x0=%g xlast=%g "
            "y0=%g ylast=%g\n",
            L,
            rows,
            cols,
            nodes_x[0],
            nodes_x[rows - 1],
            nodes_y[0],
            nodes_y[cols - 1]);
        prkmatrix r_ = b_aca_rkmatrix(eps, d, dim, rows, cols, nodes_x, nodes_y, test_function, &u_res, &v_res, &rank_inc);
        int hist_len = 0;
        int rank_sum = 0;

        /* determine history length from cumulative rank increments */
        while (rank_sum < r_->kt)
        {
            int inc = (int)round(rank_inc[hist_len]);
            rank_sum += inc;
            hist_len++;
        }
        if (rank_sum != r_->kt)
        {
            printf(
                "rank history mismatch: "
                "sum(rank_inc)=%d kt=%d\n",
                rank_sum,
                r_->kt);
        }

        node->rank_len = hist_len;
        node->u = allocate_doubles(hist_len);
        node->v = allocate_doubles(hist_len);
        node->rank_inc = allocate_doubles(hist_len);

        memcpy(node->u, u_res, hist_len * sizeof(double));

        memcpy(node->v, v_res, hist_len * sizeof(double));

        memcpy(node->rank_inc, rank_inc, hist_len * sizeof(double));

        free(u_res);
        free(v_res);
        free(rank_inc);
        *U = new_fullmatrix(rows, r_->kt);
        pfullmatrix V_ = new_fullmatrix(cols, r_->kt);
        memcpy((*U)->e, r_->a, rows * r_->kt * sizeof(double));
        memcpy(V_->e, r_->b, r_->kt * cols * sizeof(double));
        *V = transpose_fullmatrix(V_);

        (*S) = new_zero_fullmatrix(r_->kt, r_->kt);
        for (int i = 0; i < r_->kt; i++)
            (*S)->e[i * (*S)->rows + i] = 1.0;
        (*r) = r_->kt;
        del_fullmatrix(V_);
        del_rkmatrix(r_);
        node->time = omp_get_wtime() - start;
        return;
    }
    if (L > 0)
    {
        double setup_start = omp_get_wtime();

        int row1 = rows / 2;
        int row2 = rows - row1;

        int col1 = cols / 2;
        int col2 = cols - col1;

        double *X1 = new_subnodes(nodes_x, dim, rows, 0, row1);
        double *X2 = new_subnodes(nodes_x, dim, rows, row1, row2);

        double *Y1 = new_subnodes(nodes_y, dim, cols, 0, col1);
        double *Y2 = new_subnodes(nodes_y, dim, cols, col1, col2);
        for (int i = 0; i < 4; i++)
        {
            node->child[i] =
                new_residual_node();

            node->child[i]->parent =
                node;
        }
        node->child[0]->row_start = node->row_start;
        node->child[0]->col_start = node->col_start;

        node->child[1]->row_start = node->row_start;
        node->child[1]->col_start = node->col_start + col1;

        node->child[2]->row_start = node->row_start + row1;
        node->child[2]->col_start = node->col_start;

        node->child[3]->row_start = node->row_start + row1;
        node->child[3]->col_start = node->col_start + col1;
        node->child[0]->row_size = row1;
        node->child[0]->col_size = col1;

        node->child[1]->row_size = row1;
        node->child[1]->col_size = col2;

        node->child[2]->row_size = row2;
        node->child[2]->col_size = col1;

        node->child[3]->row_size = row2;
        node->child[3]->col_size = col2;
        for (int i = 0; i < 4; i++)
            node->child[i]->level = L - 1;
        pfullmatrix U11, S11, V11;
        pfullmatrix U12, S12, V12;
        pfullmatrix U21, S21, V21;
        pfullmatrix U22, S22, V22;
        int r11, r12, r21, r22;
        double setup_time = omp_get_wtime() - setup_start;

#pragma omp parallel sections num_threads(4)
        {
#pragma omp section
            {
                h_b_aca_rkmatrix(eps, d, L - 1, node->child[0], &U11, &S11, &V11, dim, row1, col1, X1, Y1, test_function, &r11);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix(eps, d, L - 1, node->child[1], &U12, &S12, &V12, dim, row1, col2, X1, Y2, test_function, &r12);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix(eps, d, L - 1, node->child[2], &U21, &S21, &V21, dim, row2, col1, X2, Y1, test_function, &r21);
            }
#pragma omp section
            {
                h_b_aca_rkmatrix(eps, d, L - 1, node->child[3], &U22, &S22, &V22, dim, row2, col2, X2, Y2, test_function, &r22);
            }
        }
        double merge_start = omp_get_wtime();
        free(X1);
        free(X2);
        free(Y1);
        free(Y2);
        // U_1=[U11*S11, U12*S12]

        pfullmatrix U_11 = new_fullmatrix(U11->rows, S11->cols);
        pfullmatrix U_12 = new_fullmatrix(U12->rows, S12->cols);
        mul_fullmatrix(U11, S11, U_11);
        mul_fullmatrix(U12, S12, U_12);
        if (U_11->rows != U_12->rows)
        {
            printf("Error: U_11 and U_12 have different number of rows\n");
            return;
        }

        pfullmatrix U_1 = new_fullmatrix(U_11->rows, U_11->cols + U_12->cols);
        for (int j = 0; j < U_11->cols; j++)
        {
            for (int i = 0; i < U_11->rows; i++)
            {
                U_1->e[j * U_1->rows + i] = U_11->e[j * U_11->rows + i];
            }
        }
        for (int j = 0; j < U_12->cols; j++)
        {
            for (int i = 0; i < U_12->rows; i++)
            {
                U_1->e[(j + U_11->cols) * U_1->rows + i] = U_12->e[j * U_12->rows + i];
            }
        }

        // SVD(U_1,eps,&U1,&S1,&V1);
        pfullmatrix U1, S1, V1;
        SVD_truncated(U_1, eps, &U1, &S1, &V1);
        del_fullmatrix(U_11);
        del_fullmatrix(U_12);
        del_fullmatrix(U_1);
        // diag(V11,V12)
        pfullmatrix V_1 = new_zero_fullmatrix(V11->rows + V12->rows, V11->cols + V12->cols);

        for (int j = 0; j < V11->cols; j++)
        {
            for (int i = 0; i < V11->rows; i++)
            {
                V_1->e[j * V_1->rows + i] = V11->e[j * V11->rows + i];
            }
        }
        for (int j = 0; j < V12->cols; j++)
        {
            for (int i = 0; i < V12->rows; i++)
            {
                V_1->e[(j + V11->cols) * V_1->rows + (i + V11->rows)] = V12->e[j * V12->rows + i];
            }
        }

        // pfullmatrix V1 =V1*V_1;
        pfullmatrix V1_temp = new_fullmatrix(V1->rows, V_1->cols);
        mul_fullmatrix(V1, V_1, V1_temp);
        // U_2=[U21*S21, U22*S22]
        pfullmatrix U_21 = new_fullmatrix(U21->rows, S21->cols);
        pfullmatrix U_22 = new_fullmatrix(U22->rows, S22->cols);
        mul_fullmatrix(U21, S21, U_21);
        mul_fullmatrix(U22, S22, U_22);
        if (U_21->rows != U_22->rows)
        {
            printf("Error: U_21 and U_22 have different number of rows\n");
            return;
        }
        pfullmatrix U_2 = new_fullmatrix(U_21->rows, U_21->cols + U_22->cols);
        for (int j = 0; j < U_21->cols; j++)
        {
            for (int i = 0; i < U_21->rows; i++)
            {
                U_2->e[j * U_2->rows + i] = U_21->e[j * U_21->rows + i];
            }
        }
        for (int j = 0; j < U_22->cols; j++)
        {
            for (int i = 0; i < U_22->rows; i++)
            {
                U_2->e[(j + U_21->cols) * U_2->rows + i] = U_22->e[j * U_22->rows + i];
            }
        }
        // SVD(U_2,eps,&U2,&S2,&V2);

        pfullmatrix U2, S2, V2;
        SVD_truncated(U_2, eps, &U2, &S2, &V2);
        del_fullmatrix(U_21);
        del_fullmatrix(U_22);
        del_fullmatrix(U_2);
        // diag(V21,V22)
        pfullmatrix V_2 = new_zero_fullmatrix(V21->rows + V22->rows, V21->cols + V22->cols);
        for (int j = 0; j < V21->cols; j++)
        {
            for (int i = 0; i < V21->rows; i++)
            {
                V_2->e[j * V_2->rows + i] = V21->e[j * V21->rows + i];
            }
        }
        for (int j = 0; j < V22->cols; j++)
        {
            for (int i = 0; i < V22->rows; i++)
            {
                V_2->e[(j + V21->cols) * V_2->rows + (i + V21->rows)] = V22->e[j * V22->rows + i];
            }
        }
        // TODO: pfullmatrix V2 =V2*V_2;
        pfullmatrix V2_temp = new_fullmatrix(V2->rows, V_2->cols);
        mul_fullmatrix(V2, V_2, V2_temp);
        del_fullmatrix(V_2);

        del_fullmatrix(U11);
        del_fullmatrix(S11);
        del_fullmatrix(V11);

        del_fullmatrix(U12);
        del_fullmatrix(S12);
        del_fullmatrix(V12);

        del_fullmatrix(U21);
        del_fullmatrix(S21);
        del_fullmatrix(V21);

        del_fullmatrix(U22);
        del_fullmatrix(S22);
        del_fullmatrix(V22);
        /* --------------------------------------------------------------------*/
        /*                       after for loop                                */
        /* --------------------------------------------------------------------*/

        // U_ =diag(U1,U2);
        pfullmatrix U_ = new_zero_fullmatrix(U1->rows + U2->rows, U1->cols + U2->cols);
        for (int j = 0; j < U1->cols; j++)
        {
            for (int i = 0; i < U1->rows; i++)
            {
                U_->e[j * U_->rows + i] = U1->e[j * U1->rows + i];
            }
        }
        for (int j = 0; j < U2->cols; j++)
        {
            for (int i = 0; i < U2->rows; i++)
            {
                U_->e[(j + U1->cols) * U_->rows + (i + U1->rows)] = U2->e[j * U2->rows + i];
            }
        }

        // V_=[S1*V1; S2*V2]

        pfullmatrix S1V1 = new_fullmatrix(S1->rows, V1_temp->cols);

        mul_fullmatrix(S1, V1_temp, S1V1);
        pfullmatrix S2V2 = new_fullmatrix(S2->rows, V2_temp->cols);

        mul_fullmatrix(S2, V2_temp, S2V2);

        pfullmatrix V_ = new_zero_fullmatrix(S1V1->rows + S2V2->rows, S1V1->cols);
        for (int j = 0; j < S1V1->cols; j++)
        {
            for (int i = 0; i < S1V1->rows; i++)
            {
                V_->e[j * V_->rows + i] =
                    S1V1->e[j * S1V1->rows + i];
            }
        }

        for (int j = 0; j < S2V2->cols; j++)
        {
            for (int i = 0; i < S2V2->rows; i++)
            {
                V_->e[j * V_->rows +
                      (i + S1V1->rows)] =
                    S2V2->e[j * S2V2->rows + i];
            }
        }
        // SVD(V_,eps,&U,&S,&V);
        pfullmatrix Uloc;
        pfullmatrix Sloc;
        pfullmatrix Vloc;

        SVD_truncated(V_, eps, &Uloc, &Sloc, &Vloc);
        del_fullmatrix(V_);
        *S = Sloc;
        *V = Vloc;
        (*r) = (*S)->rows;

        // *U=U*U_;

        pfullmatrix U_temp = new_fullmatrix(U_->rows, Uloc->cols);
        mul_fullmatrix(U_, Uloc, U_temp);

        *U = U_temp;
        del_fullmatrix(U_);
        del_fullmatrix(Uloc);

        del_fullmatrix(U1);
        del_fullmatrix(S1);
        del_fullmatrix(V1);

        del_fullmatrix(U2);
        del_fullmatrix(S2);
        del_fullmatrix(V2);

        del_fullmatrix(V1_temp);
        del_fullmatrix(V2_temp);
        double merge_time = omp_get_wtime() - merge_start;

        double total = setup_time + merge_time;

        printf("setup=%0.15f merge=%0.15f total=%0.15f\n",
               setup_time,
               merge_time,
               total);

        node->time = total;

        printf("stored=%0.15f\n", node->time);
        return;
    }
}