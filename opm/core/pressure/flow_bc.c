/*
  Copyright 2010, 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/pressure/flow_bc.h>


/* ---------------------------------------------------------------------- */
/* Compute an appropriate array dimension to minimise total number of
 * (re-)allocations. */
/* ---------------------------------------------------------------------- */
static size_t
alloc_size(size_t n, size_t c)
/* ---------------------------------------------------------------------- */
{
    if (c < n) {
        c *= 2;                 /* log_2(n) allocations */

        if (c < n) {
            c = n;              /* Typically for the first few allocs */
        }
    }

    return c;
}


/* ---------------------------------------------------------------------- */
/* Put structure in a well-defined, initial state */
/* ---------------------------------------------------------------------- */
static void
initialise_structure(struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    fbc->nbc  = 0;
    fbc->cpty = 0;

    fbc->type  = NULL;
    fbc->value = NULL;
    fbc->face  = NULL;
}


/* ---------------------------------------------------------------------- */
static int
expand_tables(size_t                         nbc,
              struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    int     ok;
    size_t  alloc_sz;
    void   *p1, *p2, *p3;

    ok = nbc <= fbc->cpty;

    if (! ok) {
        alloc_sz = alloc_size(nbc, fbc->cpty);

        p1 = realloc(fbc->type , alloc_sz * sizeof *fbc->type );
        p2 = realloc(fbc->value, alloc_sz * sizeof *fbc->value);
        p3 = realloc(fbc->face , alloc_sz * sizeof *fbc->face );

        ok = (p1 != NULL) && (p2 != NULL) && (p3 != NULL);

        if (ok) {
            fbc->type  = p1;
            fbc->value = p2;
            fbc->face  = p3;

            fbc->cpty  = alloc_sz;
        } else {
            free(p3);  free(p2);  free(p1);
        }
    }

    return ok;
}


/* ======================================================================
 * Public interface below separator
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
/* Allocate a 'FlowBoundaryConditions' structure, initially capable of
 * managing 'nbc' individual boundary conditions.  */
/* ---------------------------------------------------------------------- */
struct FlowBoundaryConditions *
flow_conditions_construct(size_t nbc)
/* ---------------------------------------------------------------------- */
{
    int                            ok;
    struct FlowBoundaryConditions *fbc;

    fbc = malloc(1 * sizeof *fbc);

    if (fbc != NULL) {
        initialise_structure(fbc);

        ok = expand_tables(nbc, fbc);

        if (! ok) {
            flow_conditions_destroy(fbc);

            fbc = NULL;
        }
    }

    return fbc;
}


/* ---------------------------------------------------------------------- */
/* Release memory resources managed by 'fbc', including the containing
 * 'struct' pointer, 'fbc'. */
/* ---------------------------------------------------------------------- */
void
flow_conditions_destroy(struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    if (fbc != NULL) {
        free(fbc->face );
        free(fbc->value);
        free(fbc->type );
    }

    free(fbc);
}


/* ---------------------------------------------------------------------- */
/* Append a new boundary condition to existing set.
 *
 * Return one (1) if successful, and zero (0) otherwise. */
/* ---------------------------------------------------------------------- */
int
flow_conditions_append(enum FlowBCType                type ,
                       int                            face ,
                       double                         value,
                       struct FlowBoundaryConditions *fbc  )
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = expand_tables(fbc->nbc + 1, fbc);

    if (ok) {
        fbc->type [ fbc->nbc ] = type ;
        fbc->value[ fbc->nbc ] = value;
        fbc->face [ fbc->nbc ] = face ;

        fbc->nbc += 1;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
/* Clear existing set of boundary conditions */
/* ---------------------------------------------------------------------- */
void
flow_conditions_clear(struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    assert (fbc != NULL);

    fbc->nbc = 0;
}
