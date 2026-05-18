#include "fullmatrix.h"

int SVD(pfullmatrix A, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V);
int SVD_truncated(pfullmatrix A, double epsilon, pfullmatrix *U, pfullmatrix *S, pfullmatrix *V_T);