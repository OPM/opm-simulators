/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARTITION_H_INCLUDED
#define PARTITION_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif  /* PARTITION_H_INLCUDED */
