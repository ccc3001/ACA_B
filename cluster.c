#include "cluster.h"
#include "basic.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

pcluster
new_cluster(int start, int size, int d, int sons) {
  pcluster tau;

  assert(start >= 0);
  assert(size >= sons);
  assert(d >= 0);
  assert(sons >= 0);
  
  tau = (pcluster) malloc(sizeof(cluster));
  if (tau == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  tau->start = start;
  tau->size = size;
  tau->d = d;
  tau->sons = sons;

  tau->bmin = allocate_doubles(d);
  tau->bmax = allocate_doubles(d);
  tau->son = (sons > 0 ? (pcluster*) malloc(sizeof(pcluster) * sons) : NULL);
  if (sons > 0 && tau->son == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  return tau;
}

pclustertree
new_clustertree(int nidx) {
  pclustertree ct;
  int i;

  assert(nidx >= 0);

  ct = (pclustertree) malloc(sizeof(clustertree));
  if (ct == NULL) {
    printf("Memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  ct->nidx = nidx;
  ct->clusters = 0;
  ct->depth = 0;
  ct->permute = allocate_ints(nidx);
  for (i = 0; i < nidx; i++) {
    ct->permute[i] = i;
  }
  
  return ct;
}

void
del_cluster(pcluster tau) {
  int i;
  
  assert(tau != NULL);
  assert(tau->start >= 0);
  assert(tau->size >= tau->sons);
  assert(tau->d >= 0);
  assert(tau->sons >= 0);
  assert(tau->sons > 0 || tau->son == NULL);
  assert(tau->sons == 0 || tau->son != NULL);
  assert((tau->size > 0 && tau->d > 0) || (tau->bmin == NULL && tau->bmax == NULL));
  assert(tau->size == 0 || tau->d == 0 || (tau->bmin != NULL && tau->bmax != NULL));

  if(tau->sons > 0) {
    for (i = 0; i < tau->sons; i++) {
      del_cluster(tau->son[i]);
    }
    free(tau->son);
  }

  del_doubles(tau->bmax);
  del_doubles(tau->bmin);
  free(tau);
}

void
del_clustertree(pclustertree ct) {
  assert(ct != NULL);
  assert(ct->nidx >= 0);
  assert(ct->nidx > 0 || ct->permute == NULL);
  assert(ct->nidx == 0 || ct->permute != NULL);


  del_cluster(ct->root);
  del_ints(ct->permute);
  free(ct);
}

void
set_boundingbox_cluster(pcluster tau, const double *nodes_x) {
  int i, j;
  int n = tau->size;
  int d = tau->d;

  assert(tau != NULL);
  assert(tau->size >= 0);
  assert(tau->d >= 0);
  assert((tau->size > 0 && tau->d > 0) || (tau->bmin == NULL && tau->bmax == NULL && nodes_x == NULL));
  assert(tau->size == 0 || tau->d == 0 || (tau->bmin != NULL && tau->bmax != NULL && nodes_x != NULL));
  
  for (j = 0; j < d; j++) {
     tau->bmin[j] = nodes_x[j * n];
     tau->bmax[j] = nodes_x[j * n];
  }
  for (i = 1; i < n; i++) {
    for (j = 0; j < d; j++) {
      tau->bmin[j] = min_double(tau->bmin[j], nodes_x[j*n + i]);
      tau->bmax[j] = max_double(tau->bmax[j], nodes_x[j*n + i]);
    }
  }
}

double
diameter_cluster(pccluster tau) {
  double sum;
  int i;

  assert(tau != NULL);
  assert(tau->bmax[0] >= tau->bmin[0]);
  sum = pow(tau->bmax[0] - tau->bmin[0], 2.0);

  for (i = 1; i < tau->d; i++) {
    assert(tau->bmax[i] >= tau->bmin[i]);
    sum += pow(tau->bmax[i] - tau->bmin[i], 2.0);
  }

  return sqrt(sum);
}

double
distance_cluster(pccluster tau, pccluster sigma) {
  double sum;
  int j;

  assert(tau != NULL);
  assert(sigma != NULL);
  assert(tau->d == sigma->d);

  sum = 0.0;
  for (j = 0; j < tau->d; j++) {
    assert(tau->bmax[j] >= tau->bmin[j]);
    assert(sigma->bmax[j] >= sigma->bmin[j]);

    if(tau->bmax[j] < sigma->bmin[j]) {
      sum += pow(sigma->bmin[j] - tau->bmax[j], 2.0);
    }
    else if(sigma->bmax[j] < tau->bmin[j]) {
      sum += pow(tau->bmin[j] - sigma->bmax[j], 2.0);
    }
  }

  return sqrt(sum);
}

pcluster
split_geometric_cluster(int start, int size, int d, const double *nodes_x, int nmin, int *permute) {
  pcluster tau;
  int i, j, k_max, counter_t0, counter_t1;
  int *nodes_permute;
  double k_diff, k_mid;
  double *nodes_x0, *nodes_x1;

  assert(size > nmin);
  assert(nmin >= 1);
  assert(size > 0 || permute == NULL);
  assert(size == 0 || permute != NULL);
  
  /* Create new cluster with 2 sons and compute its bounding box. */
  tau = new_cluster(start, size, d, 2);
  set_boundingbox_cluster(tau, nodes_x);
  
  /* Find spatial direction k_max to split cluster. */
  k_max = 0; 
  k_diff = tau->bmax[0] - tau->bmin[0];
  for (i = 1; i < d; i++) {
    if (tau->bmax[i] - tau->bmin[i] > k_diff) {
      k_max = i;
      k_diff = tau->bmax[i] - tau->bmin[i];
    }
  }
  
  /* Compute midpoint of the bounding box in the direction k_max
     and use it to split nodes in two sets. */
  k_mid = 0.5 * (tau->bmax[k_max] + tau->bmin[k_max]);
  counter_t0 = 0; 
  counter_t1 = size - 1;
  nodes_permute = allocate_ints(size);
  for (i = 0; i < size; i++) {
    if (nodes_x[k_max*size + i] < k_mid) { 
      nodes_permute[counter_t0++] = i;
    }
    else {
      nodes_permute[counter_t1--] = i;
    }
  }
  assert(counter_t1 + 1 == counter_t0);
  
  /* Save nodes corresponding to the same set in the same array. */
  nodes_x0 = allocate_doubles(d * counter_t0);
  nodes_x1 = allocate_doubles(d * (size - counter_t0));
  for (j = 0; j < d; j++) { 
    for (i = 0; i < counter_t0; i++) {
      nodes_x0[j*counter_t0 + i] = nodes_x[j*size + nodes_permute[i]];
    }
  }
  for (j = 0; j < d; j++) {
    for (i = counter_t0; i < size; i++) { 
      nodes_x1[j*(size - counter_t0) + (i - counter_t0)] = nodes_x[j*size + nodes_permute[i]];
    }
  }
  
  /* Update global permutation vector permute. */
  for (i = 0; i < size; i++) { 
    nodes_permute[i] = permute[start + nodes_permute[i]];
  }
  for (i = 0; i < size; i++) { 
    permute[start + i] = nodes_permute[i];
  }
  del_ints(nodes_permute);

  /* Check if the first set is small enough. If not, split it further. */
  if (counter_t0 > nmin) {
    tau->son[0] = split_geometric_cluster(start, counter_t0, d, nodes_x0, nmin, permute); 
  }
  else {
    tau->son[0] = new_cluster(start, counter_t0, d, 0);
    set_boundingbox_cluster(tau->son[0], nodes_x0);
  }
  
  /* Check if the second set is small enough. If not, split it further. */
  if (size - counter_t0 > nmin) {
    tau->son[1] = split_geometric_cluster(start + counter_t0, size - counter_t0,
                                          d, nodes_x1, nmin, permute); 
  }
  else {
    tau->son[1] = new_cluster(start + counter_t0, size - counter_t0, d, 0);
    set_boundingbox_cluster(tau->son[1], nodes_x1);
  }

  del_doubles(nodes_x0);
  del_doubles(nodes_x1);

  return tau;
}

pclustertree
create_geometric_clustertree(int n, int d, const double *nodes_x, int nmin) {
  pclustertree ct;
  
  ct = new_clustertree(n);
  if (n <= nmin) {
    ct->root = new_cluster(0, n, d, 0);
    set_boundingbox_cluster(ct->root, nodes_x);
  }
  else {
    ct->root = split_geometric_cluster(0, n, d, nodes_x, nmin, ct->permute);
  }
  
  return ct;
}

psupermatrix
build_supermatrix_from_cluster(pccluster rc, pccluster cc, double eta) {
  prkmatrix r;
  pfullmatrix f;
  psupermatrix s;
  double diam_rc, diam_cc, dist;
  int admissible, i, j;

  diam_rc = diameter_cluster(rc);
  diam_cc = diameter_cluster(cc);
  dist = distance_cluster(rc, cc);

  admissible = min_double(diam_rc, diam_cc) <= eta*dist;
  if (rc == cc) {
    /* Avoid that a 1xn or nx1 diagonal block is represented by a rkmatrix. */
    admissible = 0;
  }

  if (admissible) {
    /* An admissible leaf is represented by a rkmatrix. */
    r = new_rkmatrix(0, rc->size, cc->size);
    s = new_supermatrix(1, 1, rc->size, cc->size, r, NULL);
  }
  else {
    if (rc->sons == 0 || cc->sons == 0) {
      /* A not admissible leaf is represented by a fullmatrix. */
      f = new_fullmatrix(rc->size, cc->size);
      s = new_supermatrix(1, 1, rc->size, cc->size, NULL, f);
    }
    else {
      /* It isn't a leaf. The function is called recursively for the sons. */
      s = new_supermatrix(rc->sons, cc->sons, rc->size, cc->size, NULL, NULL);
      for (i = 0; i < rc->sons; i++) {
        for (j = 0; j < cc->sons; j++) {
	      s->s[i + j*rc->sons] = build_supermatrix_from_cluster(rc->son[i], cc->son[j], eta);
        }
      }
    }
  }
      
  return s;
}

void
reorder_nodes(int n, int d, double *nodes_x, const int *permute) {
  int i, j;
  double *nodes_new;

  assert(n >= 0);
  assert(d >= 0);
  assert((n > 0 && d > 0) || nodes_x == NULL);
  assert(n == 0 || d == 0 || nodes_x != NULL);
  assert(n > 0 || permute == NULL);
  assert(n == 0 || permute != NULL);
  
  nodes_new = allocate_doubles(n * d);
  for (j = 0; j < d; j++) {
    for (i = 0; i < n; i++) {
      nodes_new[j*n + i] = nodes_x[j*n + permute[i]];
    }
  }
  for (i = 0; i < d * n; i++) {
    nodes_x[i] = nodes_new[i];
  }
  
  del_doubles(nodes_new);
}

void
output_cluster(pccluster tau, int level) {
  int i;

  assert(tau != NULL);
  assert(tau->start >= 0);
  assert(tau->size >= tau->sons);
  assert(tau->sons >= 0);
  assert(tau->sons > 0 || tau->son == NULL);
  assert(tau->sons == 0 || tau->son != NULL);
  assert(level >= 0);
  
  if (tau->sons == 0) {
    printf("Leaf: ");
  }
  
  printf("level %d cluster: {%d, ..., %d}\n", level, tau->start, tau->start + tau->size - 1);

  for (i = 0; i < tau->sons; i++) {
    output_cluster(tau->son[i], level + 1);
  }
}

void
output_clustertree(pcclustertree ct) {
  assert(ct != NULL);
  assert(ct->root != NULL);
  
  output_cluster(ct->root, 0); 
}
