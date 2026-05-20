#include "fullmatrix.h"
#include "rkmatrix.h"

int *random_unique(int n, int d);

void QR(pfullmatrix W, double eps, pfullmatrix *Q_out, pfullmatrix *T_out,
        int **Jbar_out, int *r_out);

pfullmatrix
build_U(pfullmatrix C, int *Jbar, int r);

pfullmatrix
build_V(pfullmatrix QtR, int r);

void LRID(pfullmatrix C, pfullmatrix W, pfullmatrix R, double eps,
          pfullmatrix *U_out, pfullmatrix *V_out, int *r_out, int **Jbar_out);

pfullmatrix
cholesky(pcfullmatrix A);

double
frob_T1_T2(pcfullmatrix T1, pcfullmatrix T2);

double
LRnormUp(prkmatrix R, pfullmatrix U_bar, pfullmatrix V_bar,
         double nu, double nu_bar);

prkmatrix
b_aca_rkmatrix_new(double eps, int d, pcfullmatrix A, double **residuals);

void h_b_aca_rkmatrix_new(double eps, int d, int L, pcfullmatrix A, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V, int *r);
