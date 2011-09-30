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
#include <stdio.h>
#include <stdlib.h>

#include "sparse_sys.h"


/* ---------------------------------------------------------------------- */
struct CSRMatrix *
csrmatrix_new_count_nnz(size_t m)
/* ---------------------------------------------------------------------- */
{
    size_t            i;
    struct CSRMatrix *new;

    assert (m > 0);

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->ia = malloc((m + 1) * sizeof *new->ia);

        if (new->ia != NULL) {
            for (i = 0; i < m + 1; i++) { new->ia[i] = 0; }

            new->m   = m;
            new->nnz = 0;

            new->ja  = NULL;
            new->sa  = NULL;
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
    return *(const int * const)a0 - *(const int * const)b0;
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
csrmatrix_elm_index(int i, int j, const struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    int *p;

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


/* ---------------------------------------------------------------------- */
void
csrmatrix_write(const struct CSRMatrix *A, const char *fn)
/* ---------------------------------------------------------------------- */
{
    FILE *fp;

    fp = fopen(fn, "wt");

    if (fp != NULL) {
        csrmatrix_write_stream(A, fp);
    }

    fclose(fp);
}


/* ---------------------------------------------------------------------- */
void
csrmatrix_write_stream(const struct CSRMatrix *A, FILE *fp)
/* ---------------------------------------------------------------------- */
{
    size_t  i, j;
    for (i = j = 0; i < A->m; i++) {
        for (; j < (size_t) (A->ia[i + 1]); j++) {
            fprintf(fp, "%lu %lu %26.18e\n",
                    (unsigned long) (i + 1),
                    (unsigned long) (A->ja[j] + 1),
                    A->sa[j]);
        }
    }
}


/* ---------------------------------------------------------------------- */
void
vector_write(size_t n, const double *v, const char *fn)
/* ---------------------------------------------------------------------- */
{
    FILE *fp;

    fp = fopen(fn, "wt");

    if (fp != NULL) {
        vector_write_stream(n, v, fp);
    }
    
    fclose(fp);
}


/* ---------------------------------------------------------------------- */
void
vector_write_stream(size_t n, const double *v, FILE *fp)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < n; i++) {
        fprintf(fp, "%26.18e\n", v[i]);
    }
}
