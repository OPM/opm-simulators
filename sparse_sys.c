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

#include <assert.h>
#include <stdlib.h>

#include "sparse_sys.h"


/* ---------------------------------------------------------------------- */
struct CSRMatrix *
csrmatrix_new_count_nnz(size_t m)
/* ---------------------------------------------------------------------- */
{
    struct CSRMatrix *new;

    assert (m > 0);

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->ia = malloc((m + 1) * sizeof *new->ia);

        if (new->ia != NULL) {
            new->m   = m;
            new->n   = 0;
            new->nnz = 0;

            new->ja  = NULL;
            new->sa  = NULL;

            /* MAT_SIZE_T might not be 'int' so no memset() here */
            for (m = 0; m <= new->m; m++) { new->ia[m] = 0; }
        } else {
            csrmatrix_delete(new);
            new = NULL;
        }
    }

    return new;
}


/* Allocate CSR matrix, known nnz.  Allocation only.  Caller must
 * build sparsity structure before using in global assembly.
 *
 * Returns fully allocated structure if successful and NULL otherwise. */
/* ---------------------------------------------------------------------- */
struct CSRMatrix *
csrmatrix_new_known_nnz(size_t m, size_t nnz)
/* ---------------------------------------------------------------------- */
{
    struct CSRMatrix *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->ia = malloc((m + 1) * sizeof *new->ia);
        new->ja = malloc(nnz     * sizeof *new->ja);
        new->sa = malloc(nnz     * sizeof *new->sa);

        if ((new->ia == NULL) || (new->ja == NULL) || (new->sa == NULL)) {
            csrmatrix_delete(new);
            new = NULL;
        } else {
            new->m   = m;
            new->n   = 0;
            new->nnz = nnz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
size_t
csrmatrix_new_elms_pushback(struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    assert (A->ia[0] == 0);     /* Elems for row 'i' in bin i+1 ... */

    for (i = 1; i <= A->m; i++) {
        A->ia[0] += A->ia[i];
        A->ia[i]  = A->ia[0] - A->ia[i];
    }

    A->nnz = A->ia[0];
    assert (A->nnz > 0);        /* Else not a real system. */

    A->ia[0] = 0;

    A->ja = malloc(A->nnz * sizeof *A->ja);
    A->sa = malloc(A->nnz * sizeof *A->sa);

    if ((A->ja == NULL) || (A->sa == NULL)) {
        free(A->sa);   A->sa = NULL;
        free(A->ja);   A->ja = NULL;

        A->nnz = 0;
    }

    return A->nnz;
}


/* ---------------------------------------------------------------------- */
static int
cmp_row_elems(const void *a0, const void *b0)
/* ---------------------------------------------------------------------- */
{
    MAT_SIZE_T a, b;

    a = *(const MAT_SIZE_T * const) a0;
    b = *(const MAT_SIZE_T * const) b0;

    return (int)a - (int)b;
}


/* ---------------------------------------------------------------------- */
void
csrmatrix_sortrows(struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    /* O(A->nnz * log(average nnz per row)) \approx O(A->nnz) */
    for (i = 0; i < A->m; i++) {
        qsort(A->ja        + A->ia[i] ,
              A->ia[i + 1] - A->ia[i] ,
              sizeof A->ja  [A->ia[i]],
              cmp_row_elems);
    }
}


/* ---------------------------------------------------------------------- */
size_t
csrmatrix_elm_index(size_t i, MAT_SIZE_T j, const struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    MAT_SIZE_T *p;

    p = bsearch(&j, A->ja + A->ia[i], A->ia[i + 1] - A->ia[i],
                sizeof A->ja[A->ia[i]], cmp_row_elems);

    assert (p != NULL);

    return p - A->ja;
}


/* ---------------------------------------------------------------------- */
void
csrmatrix_delete(struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    if (A != NULL) {
        free(A->sa);
        free(A->ja);
        free(A->ia);
    }

    free(A);
}


/* ---------------------------------------------------------------------- */
void
csrmatrix_zero(struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < A->nnz; i++) { A->sa[i] = 0.0; }
}


/* ---------------------------------------------------------------------- */
/* v = zeros([n, 1]) */
/* ---------------------------------------------------------------------- */
void
vector_zero(size_t n, double *v)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] = 0.0; }
}
