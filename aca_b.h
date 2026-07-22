#include "fullmatrix.h"
#include "rkmatrix.h"
#include <cjson/cJSON.h>
typedef struct ACAResidualNode
{
    /* geometry */
    int row_start;
    int row_size;

    int col_start;
    int col_size;

    int level;

    /* hierarchy */
    struct ACAResidualNode *parent;
    struct ACAResidualNode *child[4];

    /* ACA history */
    double *u;
    double *v;
    double *rank_inc;
    int rank_len;
    double time;

} ACAResidualNode;
cJSON *aca_residual_node_to_json(ACAResidualNode *node);
ACAResidualNode *new_residual_node(void);
void del_residual_node(ACAResidualNode *node);

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
b_aca_rkmatrix_new(double eps, int d, pcfullmatrix A, double **residuals_u, double **residuals_v, double **rank_increase);

void h_b_aca_rkmatrix_new(double eps, int d, int L, pcfullmatrix A, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V, int *r);

prkmatrix
b_aca_rkmatrix(double eps, int d, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y,
               double (*test_function)(int dim, const double *x, const double *y), double **residuals_u, double **residuals_v, double **rank_increase);

double *new_subnodes(
    const double *nodes,
    int d,
    int n,
    int start,
    int size);
void h_b_aca_rkmatrix(
    double eps, int d, int L, ACAResidualNode *node, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V, int dim, int rows, int cols, const double *nodes_x, const double *nodes_y, double (*test_function)(int dim, const double *x, const double *y), int *r);
