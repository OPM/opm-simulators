#ifndef PARTITION_H_INCLUDED
#define PARTITION_H_INCLUDED

int
partition_unif_idx(int ndims, int nc,
                   const int *fine_d,
                   const int *coarse_d,
                   const int *idx,
                   int *p);

int
partition_compress(int n, int *p);

#endif  /* PARTITION_H_INLCUDED */
