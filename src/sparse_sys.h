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

#ifndef OPM_SPARSE_SYS_HEADER_INCLUDED
#define OPM_SPARSE_SYS_HEADER_INCLUDED

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


struct CSRMatrix
{
    size_t      m;
    size_t      nnz;

    int        *ia;
    int        *ja;

    double     *sa;
};



struct CSRMatrix *
csrmatrix_new_count_nnz(size_t m);

struct CSRMatrix *
csrmatrix_new_known_nnz(size_t m, size_t nnz);

size_t
csrmatrix_new_elms_pushback(struct CSRMatrix *A);

size_t
csrmatrix_elm_index(int i, int j, const struct CSRMatrix *A);

void
csrmatrix_sortrows(struct CSRMatrix *A);

void
csrmatrix_delete(struct CSRMatrix *A);

void
csrmatrix_zero(struct CSRMatrix *A);

/* ---------------------------------------------------------------------- */
/* v = zeros([n, 1]) */
/* ---------------------------------------------------------------------- */
void
vector_zero(size_t n, double *v);

void
csrmatrix_write(const struct CSRMatrix *A, const char *fn);

void
vector_write(size_t n, const double *v, const char *fn);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_SPARSE_SYS_HEADER_INCLUDED */
