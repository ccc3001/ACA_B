#include "supermatrix.h"
#include "aca.h"
#include "basic.h"
#include "interpolation.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

psupermatrix 
new_supermatrix(int block_rows, int block_cols, int rows, int cols,
                pcrkmatrix r, pcfullmatrix f) {
  psupermatrix s;
  
  assert(block_rows >= 0);
  assert(block_cols >= 0);
  assert(rows >= 0);
  assert(cols >= 0);
  assert(r == NULL || f == NULL);
  assert( (r == NULL && f == NULL) || (block_rows <= 1 && block_cols <= 1) );
  
  s = (psupermatrix) malloc(sizeof(supermatrix));
  if (s == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
  
  s->block_rows = block_rows;
  s->block_cols = block_cols;
  s->rows = rows;
  s->cols = cols;
  s->r = (prkmatrix) r;
  s->f = (pfullmatrix) f;
  
  if (r == NULL && f == NULL) { 
    s->s = (psupermatrix*) malloc(sizeof(psupermatrix) * block_rows * block_cols);
    if (s->s == NULL) {
      printf("Memory allocation failed.\n");
      exit(EXIT_FAILURE);
    }
  }
  else {
    s->s = NULL;
  }
  
  return s;
}

psupermatrix 
new_hp_supermatrix(int p, int k) {
  psupermatrix hp;
  pfullmatrix f;
  prkmatrix r1, r2;
  int n_half;

  assert(p >= 0);
  assert(k >= 0);
  
  if (p == 0) {
    /* Diagonal block is 1x1 fullmatrix. */
    f = new_fullmatrix(1, 1);
    hp = new_supermatrix(1, 1, 1, 1, NULL, f);
  }
  else {
    n_half = pow(2, p - 1);
    hp = new_supermatrix(2, 2, 2 * n_half, 2 * n_half, NULL, NULL);
    hp->s[0] = new_hp_supermatrix(p - 1, k);
    r1 = new_rkmatrix(k, n_half, n_half);
    hp->s[1] = new_supermatrix(1, 1, n_half, n_half, r1, NULL);
    r2 = new_rkmatrix(k, n_half, n_half);
    hp->s[2] = new_supermatrix(1, 1, n_half, n_half, r2, NULL);
    hp->s[3] = new_hp_supermatrix(p - 1, k);
  }

  return hp;
}

void
del_supermatrix(psupermatrix s) {
  int i;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );

  if (s->s != NULL) {
    for (i = 0; i < s->block_rows * s->block_cols; i++) {
      del_supermatrix(s->s[i]);
    }
    free(s->s);
  }
  else if (s->r != NULL) {
    del_rkmatrix(s->r);
  }
  else {
    del_fullmatrix(s->f);
  }
  
  free(s);
}

size_t
getsize_supermatrix(pcsupermatrix s) {
  size_t size;
  int i;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  
  size = sizeof(supermatrix); 
  
  if (s->s != NULL) {
    for (i = 0; i < s->block_rows * s->block_cols; i++) {
      size += getsize_supermatrix(s->s[i]);
    }
  }
  else if (s->r != NULL) {
    size += getsize_rkmatrix(s->r->kt, s->r->rows, s->r->cols);
  }
  else {
    size += getsize_fullmatrix(s->f->rows, s->f->cols);
  }
  
  return size;
}

void 
fill_supermatrix(psupermatrix s, double v) {
  int i;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  
  if (s->s != NULL) {
    for (i = 0; i < s->block_rows * s->block_cols; i++) {
      fill_supermatrix(s->s[i], v);
    }
  }
  else if (s->r != NULL) {
    fill_entries(s->r->rows * s->r->k, s->r->a, v);
    fill_entries(s->r->cols * s->r->k, s->r->b, v);
    s->r->kt = s->r->k;
  }
  else {
    fill_entries(s->f->rows * s->f->cols, s->f->e, v);
  }
}

void 
fill_random_supermatrix(psupermatrix s, int seed) {
  int i;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  
  if (s->s != NULL) {
    for (i = 0; i < s->block_rows * s->block_cols; i++) {
      fill_random_supermatrix(s->s[i], seed + 2);
    }
  }
  else if (s->r != NULL) {
    fill_random_entries(s->r->rows * s->r->k, s->r->a, seed);
    fill_random_entries(s->r->cols * s->r->k, s->r->b, seed + 1);
    s->r->kt = s->r->k;
  }
  else {
    fill_random_entries(s->f->rows * s->f->cols, s->f->e, seed);
  }
}

void
fill_kernelinterpolation_hp_supermatrix(psupermatrix hp, int degree, const double *nodes_x,
                                        double (*test_function)(int d, const double *x, const double *y)) {
  int n_half;
  
  assert(hp->r == NULL);
  
  if (hp->s != NULL) {
    n_half = hp->rows / 2;
    /* Fill diagonal blocks recursively. */
    fill_kernelinterpolation_hp_supermatrix(hp->s[0], degree, nodes_x, test_function);
    fill_kernelinterpolation_hp_supermatrix(hp->s[3], degree, nodes_x + n_half, test_function);
    /* Fill offdiagonal Rk blocks using interpolation. */
    del_rkmatrix(hp->s[1]->r);
    hp->s[1]->r = interpolate_rkmatrix(1, n_half, n_half, nodes_x + n_half,
                                       nodes_x, degree, test_function);
    del_rkmatrix(hp->s[2]->r);
    hp->s[2]->r = interpolate_rkmatrix(1, n_half, n_half, nodes_x, nodes_x + n_half, 
                                       degree, test_function);
  }
  else {
    /* Evaluate test_function(nodes_x, nodes_x) to fill full 1x1 diagonal block. */
    hp->f->e[0] = test_function(1, nodes_x, nodes_x);
  }
}

void
fill_kernelinterpolation_supermatrix(psupermatrix s, int d, int degree, const double *nodes_x, 
                                     double (*test_function)(int d, const double *x, const double *y), 
                                     int offset_row, int offset_col, int total_size) {
  int i, j, offset_row_new, offset_col_new;
  double *nodes_row, *nodes_col;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(d >= 0);
  assert((s->rows > 0 && s->cols > 0 && d > 0) || nodes_x == NULL);
  assert(s->rows == 0 || s->cols == 0 || d == 0 || nodes_x != NULL);
  assert(offset_row >= 0);
  assert(offset_col >= 0);
  assert(total_size >= max(offset_row, offset_col));
  assert(total_size >= max(s->rows, s->cols));
  
  if (s->s != NULL) {
    /* Fill supermatrix blocks recursively. */
    offset_col_new = offset_col;
    for (i = 0; i < s->block_cols; i++) {
      offset_row_new = offset_row;
      for (j = 0; j < s->block_rows; j++) {
        fill_kernelinterpolation_supermatrix(s->s[i*s->block_rows + j], d, 
                                             degree, nodes_x, test_function, 
                                             offset_row_new, offset_col_new,
                                             total_size);
        offset_row_new += s->s[i*s->block_rows + j]->rows;
      }
      offset_col_new += s->s[i * s->block_rows]->cols;
    }
  }
  else {
    nodes_row = allocate_doubles(d * s->rows);
    nodes_col = allocate_doubles(d * s->cols);
    
    /* Copy only the needed entries of nodes_x. */
    for (i = 0; i < s->rows; i++) {
      for (j = 0; j < d; j++) {
        nodes_row[j*s->rows + i] = nodes_x[j*total_size + offset_row + i];
      }
    }
    for (i = 0; i < s->cols; i++) {
      for (j = 0; j < d; j++) {
        nodes_col[j*s->cols + i] = nodes_x[j*total_size + offset_col + i];
      }
    }

    if (s->r != NULL) {
      /* Fill rkmatrix block using interpolation. */
      del_rkmatrix(s->r);
      s->r = interpolate_rkmatrix(d, s->rows, s->cols, nodes_row,
                                  nodes_col, degree, test_function);
    }
    else {
      /* Fill fullmatrix block through evaluation of the test_function. */
      fill_fullmatrix(s->f, d, nodes_row, nodes_col, test_function);
    }

    del_doubles(nodes_row);
    del_doubles(nodes_col);
  }
}

void
fill_kernelaca_supermatrix(psupermatrix s, int d, int k, const double *nodes_x, 
                           double (*test_function)(int d, const double *x, const double *y), 
                           int offset_row, int offset_col, int total_size) {
  int i, j, offset_row_new, offset_col_new;
  double *nodes_row, *nodes_col;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(d >= 0);
  assert((s->rows > 0 && s->cols > 0 && d > 0) || nodes_x == NULL);
  assert(s->rows == 0 || s->cols == 0 || d == 0 || nodes_x != NULL);
  assert(offset_row >= 0);
  assert(offset_col >= 0);
  assert(total_size >= max(offset_row, offset_col));
  assert(total_size >= max(s->rows, s->cols));
  
  if (s->s != NULL) {
    /* Fill supermatrix blocks recursively. */  
    offset_col_new = offset_col;
    for (i = 0; i < s->block_cols; i++) {
      offset_row_new = offset_row;
      for (j = 0; j < s->block_rows; j++) {
        fill_kernelaca_supermatrix(s->s[i*s->block_rows + j], d, k, nodes_x,
                                   test_function, offset_row_new, offset_col_new,
                                   total_size);
        offset_row_new += s->s[i*s->block_rows + j]->rows;
      }
      offset_col_new += s->s[i * s->block_rows]->cols;
    }
  }
  else {
    nodes_row = allocate_doubles(d * s->rows);
    nodes_col = allocate_doubles(d * s->cols);
    
    /* Copy only the needed entries of nodes_x. */
    for (i = 0; i < s->rows; i++) {
      for (j = 0; j < d; j++) {
        nodes_row[j*s->rows + i] = nodes_x[j*total_size + offset_row + i];
      }
    }
    for (i = 0; i < s->cols; i++) {
      for (j = 0; j < d; j++) {
        nodes_col[j*s->cols + i] = nodes_x[j*total_size + offset_col + i];
      }
    }

    if (s->r != NULL) {
      del_rkmatrix(s->r);
      if (min(s->rows, s->cols) > k) {
        /* Fill rkmatrix block using ACA. */
        s->r = aca_rkmatrix(d, k, s->rows, s->cols, nodes_row, nodes_col, test_function);
      }
      else {
        /* The rank k is too large. Change rkmatrix block to fullmatrix block. */
        s->r = NULL;
        s->f = new_fullmatrix(s->rows, s->cols);
      }
    }
    if (s->f != NULL) {
      /* Fill fullmatrix block through evaluation of the test_function. */
      fill_fullmatrix(s->f, d, nodes_row, nodes_col, test_function);
    }

    del_doubles(nodes_row);
    del_doubles(nodes_col);
  }
}

void
fill_kernelaca_delta_supermatrix(psupermatrix s, int d, int k, double delta, const double *nodes_x, 
                                 double (*test_function)(int d, const double *x, const double *y), 
                                 int offset_row, int offset_col, int total_size) {
  int i, j, offset_row_new, offset_col_new;
  double *nodes_row, *nodes_col;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(d >= 0);
  assert((s->rows > 0 && s->cols > 0 && d > 0) || nodes_x == NULL);
  assert(s->rows == 0 || s->cols == 0 || d == 0 || nodes_x != NULL);
  assert(offset_row >= 0);
  assert(offset_col >= 0);
  assert(total_size >= max(offset_row, offset_col));
  assert(total_size >= max(s->rows, s->cols));
  
  if (s->s != NULL) {
    /* Fill supermatrix blocks recursively. */  
    offset_col_new = offset_col;
    for (i = 0; i < s->block_cols; i++) {
      offset_row_new = offset_row;
      for (j = 0; j < s->block_rows; j++) {
        fill_kernelaca_delta_supermatrix(s->s[i*s->block_rows + j], d, k, delta,
                                         nodes_x, test_function, offset_row_new,
                                         offset_col_new, total_size);
        offset_row_new += s->s[i*s->block_rows + j]->rows;
      }
      offset_col_new += s->s[i * s->block_rows]->cols;
    }
  }
  else {
    nodes_row = allocate_doubles(d * s->rows);
    nodes_col = allocate_doubles(d * s->cols);
    
    /* Copy only the needed entries of nodes_x. */
    for (i = 0; i < s->rows; i++) {
      for (j = 0; j < d; j++) {
        nodes_row[j*s->rows + i] = nodes_x[j*total_size + offset_row + i];
      }
    }
    for (i = 0; i < s->cols; i++) {
      for (j = 0; j < d; j++) {
        nodes_col[j*s->cols + i] = nodes_x[j*total_size + offset_col + i];
      }
    }

    if (s->r != NULL) {
      del_rkmatrix(s->r);
      if (min(s->rows, s->cols) > k) {
        /* Fill rkmatrix block using ACA. */
        s->r = aca_delta_rkmatrix(d, k, delta, s->rows, s->cols, nodes_row, nodes_col, test_function);
      }
      else {
        /* The rank k is too large. Change rkmatrix block to fullmatrix block. */
        s->r = NULL;
        s->f = new_fullmatrix(s->rows, s->cols);
      }
    }
    if (s->f != NULL) {
      /* Fill fullmatrix block through evaluation of the test_function. */
      fill_fullmatrix(s->f, d, nodes_row, nodes_col, test_function);
    }

    del_doubles(nodes_row);
    del_doubles(nodes_col);
  }
}

void 
clear_supermatrix(psupermatrix s) {
  int i;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  
  if (s->s != NULL) {
    for (i = 0; i < s->block_rows * s->block_cols; i++) {
      clear_supermatrix(s->s[i]);
    }
  }
  else if (s->r != NULL) {
    clear_entries(s->r->rows * s->r->k, s->r->a);
    clear_entries(s->r->cols * s->r->k, s->r->b);
    s->r->kt = 0;
  }
  else {
    clear_entries(s->f->rows * s->f->cols, s->f->e);
  }
}

void 
addeval_supermatrix(pcsupermatrix s, const double *v, double *w) {
  int i, j, ij, offsetv, offsetw;
  
  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  
  if (s->s != NULL) {
    offsetw = 0;
    for (i = 0; i < s->block_rows; i++) {
      ij = i;
      offsetv = 0;
      for (j = 0; j < s->block_cols; j++) {
        addeval_supermatrix(s->s[ij], v + offsetv, w + offsetw);
        offsetv += s->s[ij]->cols;
        ij += s->block_rows;
      }
      offsetw += s->s[i]->rows;     
    }
  }
  else if (s->r != NULL) {
    addeval_rkmatrix(s->r, v, w);
  }
  else {
    addeval_fullmatrix(s->f, v, w);
  }
}

void 
eval_supermatrix(pcsupermatrix s, const double *v, double *w) {
  clear_entries(s->rows, w);
  addeval_supermatrix(s, v, w);
}

void
recompress_supermatrix(psupermatrix s, double eps) {
  int i, j;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(eps >= 0.0);
  
  if (s->s != NULL) {
    for (j = 0; j < s->block_cols; j++) {
      for (i = 0; i < s->block_rows; i++) { 
        recompress_supermatrix(s->s[j*s->block_rows + i], eps);
      }
    }
  }
  else if (s->r != NULL) { 
    adaptive_truncate_rkmatrix(s->r, eps);  
  }
}

void 
write_supermatrix(pcsupermatrix s, int row_offset, int col_offset) {
  int i, j, new_row_offset;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(row_offset >= 0);
  assert(col_offset >= 0);
  
  if (s->s != NULL) {
    printf("[%d, ..., %d] x [%d, ..., %d] supermatrix, %d block rows, %d block cols\n",
           row_offset, row_offset + s->rows - 1, col_offset, col_offset + s->cols - 1,
           s->block_rows, s->block_cols);
    for (i = 0; i < s->block_cols; i++) {
      new_row_offset = row_offset;
      for (j = 0; j < s->block_rows; j++) {
        write_supermatrix(s->s[i*s->block_cols + j], new_row_offset, col_offset);
        new_row_offset += s->s[i*s->block_cols + j]->rows;
      }
      col_offset += s->s[i * s->block_cols]->cols;
    }
  }
  else {
    if (s->r != NULL) {
      printf("[%d, ..., %d] x [%d, ..., %d] rkmatrix\n", row_offset,
             row_offset + s->rows - 1, col_offset, col_offset + s->cols - 1);
    }
    else {
      printf("[%d, ..., %d] x [%d, ..., %d] fullmatrix\n", row_offset,
             row_offset + s->rows - 1, col_offset, col_offset + s->cols - 1);
    }
  }
}

static void
outputleaf_rkmatrix(pcrkmatrix r, FILE *fp, int row_offset, int col_offset,
                    int row_max, int col_max) {
  int rows, cols;

  assert(r != NULL);
  assert(fp != NULL);
  assert(row_offset >= 0);
  assert(col_offset >= 0);
  assert(row_max > row_offset);
  assert(col_max > col_offset);
  
  rows = r->rows;
  cols = r->cols;

  fprintf(fp, "%d %d m\n", col_offset, row_max - row_offset - rows);
  fprintf(fp, "%d %d l\n", col_offset, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - row_offset - rows);
  fprintf(fp, "%d %d l\n", col_offset, row_max - row_offset - rows);
  fprintf(fp, "fl\n");
  fprintf(fp, "bl\n");  
  if ((r->kt < 10 && (1.0 * r->cols) / col_max > 0.05) || (r->kt < 100 && (1.0 * r->cols) / col_max > 0.1)) {  
    fprintf(fp, "%.1f %.1f m (%d) show\n", col_offset + cols * 0.2, row_max - 0.9 * rows - row_offset, r->kt); 
  }
  /* border */
  fprintf(fp, "%d %d m\n", col_offset, row_max - rows - row_offset);
  fprintf(fp, "%d %d l\n", col_offset, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - rows - row_offset);
  fprintf(fp, "%d %d l\n", col_offset, row_max - rows - row_offset);
  fprintf(fp, "cs\n");
}

static void
outputleaf_fullmatrix(pcfullmatrix f, FILE *fp, int row_offset, 
                      int col_offset, int row_max, int col_max) {
  int rows, cols;

  assert(f != NULL);
  assert(fp != NULL);
  assert(row_offset >= 0);
  assert(col_offset >= 0);
  assert(row_max > row_offset);
  assert(col_max > col_offset);

  rows = f->rows;
  cols = f->cols;

  fprintf(fp, "%.3f %.3f m\n", (double) col_offset, (double) row_max - rows - row_offset);
  fprintf(fp, "%.3f %.3f l\n", (double) col_offset, (double) row_max - rows - row_offset + rows);
  fprintf(fp, "%.3f %.3f l\n", (double) col_offset + cols, (double) row_max - row_offset);
  fprintf(fp, "%.3f %.3f l\n", (double) col_offset + cols, (double) row_max - rows - row_offset);
  fprintf(fp, "%.3f %.3f l\n", (double) col_offset, (double) row_max - rows - row_offset);
  fprintf(fp, "fi\n");
  fprintf(fp, "cs\n");

  /* border */
  fprintf(fp, "bl\n");     
  fprintf(fp, "%d %d m\n", col_offset, row_max - rows - row_offset);
  fprintf(fp, "%d %d l\n", col_offset, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - row_offset);
  fprintf(fp, "%d %d l\n", col_offset + cols, row_max - rows - row_offset);
  fprintf(fp, "%d %d l\n", col_offset, row_max - rows - row_offset);
  fprintf(fp, "cs\n");  
}

static void
outputleaf_supermatrix(pcsupermatrix s, FILE *fp, int row_offset, 
                       int col_offset, int row_max, int col_max) {
  int i, j, new_col_offset;

  assert(s != NULL);
  assert(s->r == NULL || s->f == NULL);
  assert(s->r == NULL || s->s == NULL);
  assert(s->f == NULL || s->s == NULL);
  assert( (s->r == NULL && s->f == NULL) || (s->block_rows <= 1 && s->block_cols <= 1) );
  assert(fp != NULL);
  assert(row_offset >= 0);
  assert(col_offset >= 0);
  assert(row_max > row_offset);
  assert(col_max > col_offset);
  
  if (s->s != NULL) {
    for (i = 0; i < s->block_rows; i++) {
      new_col_offset = col_offset;
      for (j = 0; j < s->block_cols; j++) {
        outputleaf_supermatrix(s->s[i + j*s->block_rows], fp, row_offset,
                               new_col_offset, row_max, col_max);
        new_col_offset += s->s[j * s->block_rows]->cols; 
      }
      row_offset += s->s[i]->rows;
    }
  }
  else if (s->r != NULL) {
    outputleaf_rkmatrix(s->r, fp, row_offset, col_offset, row_max, col_max);    
  }
  else {
    outputleaf_fullmatrix(s->f, fp, row_offset, col_offset, row_max, col_max);
  }
}

void
output_supermatrix(pcsupermatrix s, char *fname) {
  double scale;
  FILE *fp;

  assert(s != NULL);
  assert(fname != NULL);
  
  fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("Opening file failed.\n");
    exit(EXIT_FAILURE);
  }
  
  scale = 500.0 / (double) max(s->rows, s->cols);
  fprintf(fp, "%%!PS-Adobe-2.0-2.0 EPSF-2.0\n");
  fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", 0.0, 0.0,
          (double)(s->cols) * scale + 20, (double)(s->rows) * scale + 20);
  fprintf(fp, "%f dup translate\n", 10.0);
  fprintf(fp, "%f dup scale\n", scale);
  if (scale < 5) {
    fprintf(fp, "%f setlinewidth\n", 0.02 * scale);
  }
  else {
    fprintf(fp, "%f setlinewidth\n", 0.1);
  }
  fprintf(fp, "/Helvetica findfont %f scalefont setfont\n", 12.0 / scale);
  fprintf(fp, "0 setgray\n");
  fprintf(fp, "/bl{0 setgray} def\n");
  fprintf(fp, "/fb{0.95 setgray fill} def\n");
  fprintf(fp, "/fr{1.0 0.6 0.6  setrgbcolor fill} def\n");
  fprintf(fp, "/fg{0.5 setgray fill} def\n");
  fprintf(fp, "/fl{0.6 1.0 0.6  setrgbcolor fill} def\n");
  fprintf(fp, "/fi{0.8 0.0 0.0  setrgbcolor fill} def\n");
  fprintf(fp, "/fu{0.3 0.3 1.0  setrgbcolor fill} def\n");
  fprintf(fp, "/cs{closepath stroke} def\n");
  fprintf(fp, "/m{moveto} def\n");
  fprintf(fp, "/l{lineto} def\n");
  fprintf(fp, "/hf{/Helvetica findfont} def\n");
  fprintf(fp, "/sf{scalefont setfont} def\n");

  outputleaf_supermatrix(s, fp, 0, 0, s->rows, s->cols);
  fprintf(fp, "showpage\n");
  fclose(fp);
}
