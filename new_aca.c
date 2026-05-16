#include "rkmatrix.h"
#include "fullmatrix.h"
#include "basic.h"
#include "blas.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
double
compute_entry_aca_new(pcrkmatrix r, int k, int row, int col, pcfullmatrix A)
{
    int i;
    double value;

    assert(r != NULL);
    assert(r->k >= k);
    assert((r->rows > 0 && r->k > 0) || r->a == NULL);
    assert(r->rows == 0 || r->k == 0 || r->a != NULL);
    assert((r->cols > 0 && r->k > 0) || r->b == NULL);
    assert(r->cols == 0 || r->k == 0 || r->b != NULL);
    assert(row >= 0);
    assert(col >= 0);

    value = A->e[col * A->rows + row];
    for (i = 0; i < k; i++)
    {
        value -= r->a[i * r->rows + row] * r->b[i * r->cols + col];
    }

    return value;
}

prkmatrix
aca_rkmatrix_new(double eps, pcfullmatrix A)
{
    int rows = A->rows;
    int cols = A->cols;
    int d = 0;
    d = min(cols, rows);
    prkmatrix r;
    r = new_zero_rkmatrix(d, rows, cols);
    double u_k_2, v_k_2;
    double v;
    double u = 0;
    int piv_col = rand() % cols;
    int piv_row, i;
    int k = 0;
    r->kt = 0;
    double max_col_value, max_row_value;
    do
    {
        // printf("iteration:%d\n",k);
        // compute all entries for a pivot col.
        max_col_value = 0.0;
        u_k_2 = 0;
        // printf("\npiv_col:%d\n",piv_col);
        for (i = 0; i < rows; i++)
        {

            r->a[k * rows + i] = compute_entry_aca_new(r, k, i, piv_col, A);
            if (fabs(max_col_value) < fabs(r->a[k * rows + i]))
            {
                max_col_value = r->a[k * rows + i];
                piv_row = i;
            }
            u_k_2 += r->a[k * rows + i] * r->a[k * rows + i];
        }
        // printf("\npiv_row:%d\n",piv_row);
        // printf("\nmax_col_value:%.4g\n",max_col_value);
        if (fabs(max_col_value) < 1e-30)
        {
            printf("max_col value is below 1e-14\n");
            break;
        }
        max_col_value = 1 / max_col_value;
        dscal_(&rows, &max_col_value, r->a + k * rows, eins_);
        u_k_2 = u_k_2 * max_col_value * max_col_value;
        // compute all entries for a pivot row
        max_row_value = 0.0;
        v_k_2 = 0;
        for (i = 0; i < cols; i++)
        {
            r->b[k * cols + i] = compute_entry_aca_new(r, k, piv_row, i, A);
            if (fabs(max_row_value) < fabs(r->b[k * cols + i]))
            {
                max_row_value = r->b[k * cols + i];
                piv_col = i;
            }
            v_k_2 += r->b[k * cols + i] * r->b[k * cols + i];
        }
        double cross_term = 0.0;

        for (int j = 0; j < k; j++)
        {
            double dot_u = 0.0;
            double dot_v = 0.0;

            for (int i = 0; i < rows; i++)
            {
                dot_u += r->a[k * rows + i] * r->a[j * rows + i];
            }

            for (int i = 0; i < cols; i++)
            {
                dot_v += r->b[j * cols + i] * r->b[k * cols + i];
            }

            cross_term += 2 * dot_u * dot_v;
        }

        v = v_k_2 * u_k_2;
        u = u + v + cross_term;

        k += 1;
        r->kt += 1;
        // printf("rkt:%d\n",r->kt);
        // printf("%.4g<%.4g\n",sqrt(u)*eps,sqrt(v));
        // printf("%d<%d\n",k,d);
    } while (sqrt(u) * eps <= sqrt(v) && k < d);

    return r;
}
