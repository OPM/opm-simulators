/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

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


int
partition_allocate_inverse(int nc, int max_blk,
                           int **pi, int **inverse);

void
partition_deallocate_inverse(int *pi, int *inverse);

void
partition_invert(int nc, const int *p,
                 int *pi, int *inverse);

void
partition_localidx(int nblk, const int *pi, const int *inverse,
                   int *localidx);


int
partition_split_disconnected(int nc, int nneigh, const int *neigh,
                             int *p);
#endif  /* PARTITION_H_INLCUDED */
