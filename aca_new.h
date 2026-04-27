#include "rkmatrix.h"
#include "fullmatrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

prkmatrix aca(int rows, int cols, pfullmatrix A, int seed, double eta)
{

    int k = MIN(rows, cols);
    prkmatrix A_aprox = new_rkmatrix(k, rows, cols);
    int piv_row, piv_col, i;
    double *x, *y, pivot_value;

    piv_row = rand() % rows;
    for (int i = 0; i < k; i++)
    {
        int max_row_value = 0;
        compute_row(A, piv_row, cols, x);
        for (int j = 0; j < cols; j++)
        {
        }
    }

    return A_aprox;
}
