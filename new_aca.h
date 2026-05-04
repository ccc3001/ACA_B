#include "rkmatrix.h"
#include "fullmatrix.h"

double 
compute_entry_aca_new(pcrkmatrix r, int k, int row, int col, pcfullmatrix A);


prkmatrix 
aca_rkmatrix_new( double eps, pcfullmatrix A);