#ifndef Cluster
#define Cluster

#include "supermatrix.h"

typedef struct _cluster cluster;
typedef cluster *pcluster;
typedef const cluster *pccluster;

typedef struct _clustertree clustertree;
typedef clustertree *pclustertree;
typedef const clustertree *pcclustertree;

/* A cluster is a node of a cluster tree. */
struct _cluster {
  int start; /* First element of the associated index set. */
  int size; /* Size of the associated index set. */
  double *bmin; /* Lower bounds of the bounding box. */
  double *bmax; /* Upper bounds of the bounding box. */
  int d; /* Spatial dimension of the bounding box (number of entries in bmin and bmax). */
  int sons; /* Number of son clusters. */
  pcluster *son; /* Pointer to son clusters. */
};

struct _clustertree {
  int nidx; /* Number of indices. */
  int clusters; /* Number of clusters. */
  int depth; /* Depth of the clustertree. */
  int *permute; /* Permutation vector from clustering. */
  pcluster root; /* Pointer to root cluster. */
};

/* Create a new cluster and don't set entries of bmin and bmax. 
   If sons = 0, then son is set to Null. If sons > 0, then
   storage is allocated for son, but the entries aren't set. */
pcluster
new_cluster(int start, int size, int d, int sons);

/* Create a new clustertree. Storage for nidx entries in
   permute is allocated. The entries are set wrt 
   the equation permute[i] = i. The pointer
   to the root cluster isn't set. */
pclustertree
new_clustertree(int nidx);

/* Delete cluster. */
void
del_cluster(pcluster tau);

/* Delete clustertree. */
void
del_clustertree(pclustertree ct);

/* Use nodes nodes_x to compute the bounding box of the cluster tau.
   Nodes in nodes_x are saved columnwise for spatial dimension, i.e.
   first all first-coordinates, then all second-coordinates, etc. */
void
set_boundingbox_cluster(pcluster tau, const double *nodes_x);

/* Compute the diameter (Euclidean norm) of the bounding box
   corresponding to the cluster tau. */
double
diameter_cluster(pccluster tau);

/* Compute the distance (Euclidean norm) between the bounding boxes
   corresponding to the clusters tau and sigma. */
double
distance_cluster(pccluster tau, pccluster sigma);

/* Create a cluster (with 2 sons) based on geometric bisection
   of the size nodes (given by nodes_x) in d-dimensional space.
   The index of the first element of the associated index set is start.
   The initial permutation vector is permute.
   If son[0]->size > nmin or son[1]->size > nmin, then the 
   function is called again on the corresponding son cluster. */
pcluster
split_geometric_cluster(int start, int size, int d, const double *nodes_x, int nmin, int *permute);

/* Create a clustertree (with maximal leaf size nmin) based on geometric
   bisection of the n nodes (given by nodes_x) in d-dimensional space. */
pclustertree
create_geometric_clustertree(int n, int d, const double *nodes_x, int nmin);

/* Create a supermatrix structure corresponding to the block rc x cc.
   Blocks satisfying min(diam(rc), diam(cc)) <= eta dist(rc, cc) and
   rc != cc are represented by rkmatrix objects. The remainining
   blocks are represented by fullmatrix objects.
   rc: Cluster corresponding to the row of the supermatrix.
   cc: Cluster corresponding to the column of the supermatrix.
   eta: Admissibility parameter. */
psupermatrix
build_supermatrix_from_cluster(pccluster rc, pccluster cc, double eta);

/* Reorder the n d-dimensional nodes in nodes_x according to the 
   permutation given by permute. nodes_x begins with all first
   components, then all second components etc. */
void
reorder_nodes(int n, int d, double *nodes_x, const int *permute);

/* Output the cluster tau and prints that it is level level. */
void
output_cluster(pccluster tau, int level);

/* Output clustertree. */
void
output_clustertree(pcclustertree ct);

#endif
