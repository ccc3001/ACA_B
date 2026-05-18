#include "fullmatrix.h"
#include "rkmatrix.h"
#include "svd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

int SVD(pfullmatrix A,
        pfullmatrix *U,
        pfullmatrix *S,
        pfullmatrix *V_T)
{
    int m = A->rows;
    int n = A->cols;
    int minmn = (m < n) ? m : n;
    int lda = m;
    int info;
    double *a_copy;

    /* allocate output matrices */

    *U = new_fullmatrix(m, m);
    *S = new_fullmatrix(m, n);
    *V_T = new_fullmatrix(n, n);

    double *sigma = malloc(minmn * sizeof(double));

    if (sigma == NULL)
    {
        printf("Allocation failed.\n");
        return -1;
    }

    /*
       LAPACK destroys input matrix,
       so create a copy
    */
    a_copy = malloc(m * n * sizeof(double));
    if (a_copy == NULL)
    {
        printf("Allocation failed.\n");
        free(sigma);
        return -1;
    }
    memcpy(a_copy, A->e, m * n * sizeof(double));

    info = LAPACKE_dgesdd(
        LAPACK_COL_MAJOR,
        'A',
        m,
        n,
        a_copy,
        lda,
        sigma,
        (*U)->e,
        m,
        (*V_T)->e,
        n);

    free(a_copy);

    if (info > 0)
    {
        printf("SVD failed.\n");

        free(sigma);

        return info;
    }

    *S = new_zero_fullmatrix(m, n);
    for (int i = 0; i < minmn; i++)
        (*S)->e[i * (*S)->rows + i] = sigma[i];

    free(sigma);

    return (int)info;
}

int SVD_truncated(
    pfullmatrix A,
    double epsilon,
    pfullmatrix *U,
    pfullmatrix *S,
    pfullmatrix *V_T)
{
    int m = A->rows;
    int n = A->cols;

    int minmn = (m < n) ? m : n;

    int lda = m;
    int info;

    double *sigma = NULL;
    double *a_copy = NULL;

    /*
        Temporary full matrices returned by LAPACK
    */

    pfullmatrix Utmp;
    pfullmatrix VTtmp;

    sigma = malloc(minmn * sizeof(double));

    if (sigma == NULL)
    {
        printf("Allocation failed.\n");
        return -1;
    }

    a_copy = malloc(m * n * sizeof(double));

    if (a_copy == NULL)
    {
        printf("Allocation failed.\n");
        free(sigma);
        return -1;
    }

    memcpy(a_copy, A->e, m * n * sizeof(double));

    /*
        Economy-size SVD
    */

    Utmp = new_fullmatrix(m, minmn);
    VTtmp = new_fullmatrix(minmn, n);

    info = LAPACKE_dgesdd(
        LAPACK_COL_MAJOR,
        'S',
        m,
        n,
        a_copy,
        lda,
        sigma,
        Utmp->e,
        m,
        VTtmp->e,
        minmn);

    free(a_copy);

    if (info > 0)
    {
        printf("SVD failed.\n");

        free(sigma);

        del_fullmatrix(Utmp);
        del_fullmatrix(VTtmp);

        return info;
    }

    /*
        Determine numerical rank
    */

    int r = 0;

    double smax = sigma[0];

    for (int i = 0; i < minmn; i++)
    {
        if (sigma[i] / smax > epsilon)
            r++;
        else
            break;
    }

    /*
        Truncate matrices
    */

    *U = new_submatrix(Utmp, 0, 0, m, r);
    *V_T = new_submatrix(VTtmp, 0, 0, r, n);

    /*
        Build diagonal Sigma matrix
    */

    *S = new_zero_fullmatrix(r, r);

    for (int i = 0; i < r; i++)
    {
        (*S)->e[i * r + i] = sigma[i];
    }

    /*
        Cleanup temporary matrices
    */

    del_fullmatrix(Utmp);
    del_fullmatrix(VTtmp);

    free(sigma);

    return r;
}