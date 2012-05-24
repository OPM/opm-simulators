/*===========================================================================
//
// File: call_umfpack.c
//
// Created: 2011-10-05 19:55:34+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <suitesparse/umfpack.h>

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/call_umfpack.h>

struct CSCMatrix {
    UF_long  n;
    UF_long  nnz;

    UF_long *p;
    UF_long *i;
    double  *x;
};


/* ---------------------------------------------------------------------- */
static void
csc_deallocate(struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    if (csc != NULL) {
        free(csc->x);
        free(csc->i);
        free(csc->p);
    }

    free(csc);
}


/* ---------------------------------------------------------------------- */
static struct CSCMatrix *
csc_allocate(UF_long n, UF_long nnz)
/* ---------------------------------------------------------------------- */
{
    struct CSCMatrix *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->p = malloc((n + 1) * sizeof *new->p);
        new->i = malloc(nnz     * sizeof *new->i);
        new->x = malloc(nnz     * sizeof *new->x);

        if ((new->p == NULL) || (new->i == NULL) || (new->x == NULL)) {
            csc_deallocate(new);
            new = NULL;
        } else {
            new->n   = n;
            new->nnz = nnz;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
csr_to_csc(const int        *ia,
           const int        *ja,
           const double     *sa,
           struct CSCMatrix *csc)
/* ---------------------------------------------------------------------- */
{
    UF_long i, nz;

    /* Clear garbage, prepare for counting */
    for (i = 0; i <= csc->n; i++) { csc->p[i] = 0; }

    /* Count column connections */
    for (nz = 0; nz < csc->nnz; nz++) {
        csc->p[ ja[nz] + 1 ] += 1;
    }

    /* Define column start pointers */
    for (i = 1; i <= csc->n; i++) {
        csc->p[0] += csc->p[i];
        csc->p[i]  = csc->p[0] - csc->p[i];
    }

    assert (csc->p[0] == csc->nnz);

    /* Fill matrix whilst defining column end pointers */
    for (i = nz = 0; i < csc->n; i++) {
        for (; nz < ia[i + 1]; nz++) {
            csc->i[ csc->p[ ja[nz] + 1 ] ] = i;      /* Insertion sort */
            csc->x[ csc->p[ ja[nz] + 1 ] ] = sa[nz]; /* Insert mat elem */

            csc->p        [ ja[nz] + 1 ]  += 1;      /* Advance col ptr */
        }
    }

    assert (csc->p[csc->n] == csc->nnz);

    csc->p[0] = 0;
}


/* ---------------------------------------------------------------------- */
static void
solve_umfpack(struct CSCMatrix *csc, const double *b, double *x)
/* ---------------------------------------------------------------------- */
{
    void *Symbolic, *Numeric;
    double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

    umfpack_dl_defaults(Control);

    umfpack_dl_symbolic(csc->n, csc->n, csc->p, csc->i, csc->x,
                        &Symbolic, Control, Info);
    umfpack_dl_numeric (csc->p, csc->i, csc->x,
                        Symbolic, &Numeric, Control, Info);

    umfpack_dl_free_symbolic(&Symbolic);

    umfpack_dl_solve(UMFPACK_A, csc->p, csc->i, csc->x, x, b,
                     Numeric, Control, Info);

    umfpack_dl_free_numeric(&Numeric);
}


/*---------------------------------------------------------------------------*/
void
call_UMFPACK(struct CSRMatrix *A, const double *b, double *x)
/*---------------------------------------------------------------------------*/
{
    struct CSCMatrix *csc;

    csc = csc_allocate(A->m, A->ia[A->m]);

    if (csc != NULL) {
        csr_to_csc(A->ia, A->ja, A->sa, csc);

        solve_umfpack(csc, b, x);
    }

    csc_deallocate(csc);
}

