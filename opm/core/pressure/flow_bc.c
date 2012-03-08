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
#include <string.h>

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
static int
initialise_structure(struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    fbc->nbc       = 0;
    fbc->cond_cpty = 0;
    fbc->face_cpty = 0;

    fbc->cond_pos  = malloc(1 * sizeof *fbc->cond_pos);
    fbc->type      = NULL;
    fbc->value     = NULL;
    fbc->face      = NULL;

    return fbc->cond_pos != NULL;
}


/* ---------------------------------------------------------------------- */
static int
expand_tables(size_t                         nbc,
              size_t                         nf ,
              struct FlowBoundaryConditions *fbc)
/* ---------------------------------------------------------------------- */
{
    int     ok_cond, ok_face;
    size_t  alloc_sz;
    void   *p1, *p2, *p3, *p4;

    ok_cond = nbc <= fbc->cond_cpty;
    ok_face = nf  <= fbc->face_cpty;

    if (! ok_cond) {
        alloc_sz = alloc_size(nbc, fbc->cond_cpty);

        p1 = realloc(fbc->type    , (alloc_sz + 0) * sizeof *fbc->type    );
        p2 = realloc(fbc->value   , (alloc_sz + 0) * sizeof *fbc->value   );
        p3 = realloc(fbc->cond_pos, (alloc_sz + 1) * sizeof *fbc->cond_pos);

        ok_cond = (p1 != NULL) && (p2 != NULL) && (p3 != NULL);

        if (p1 != NULL) { fbc->type     = p1; }
        if (p2 != NULL) { fbc->value    = p2; }
        if (p3 != NULL) { fbc->cond_pos = p3; }

        if (ok_cond) {
            fbc->cond_cpty = alloc_sz;
        }
    }

    if (! ok_face) {
        alloc_sz = alloc_size(nf, fbc->face_cpty);

        p4 = realloc(fbc->face, alloc_sz * sizeof *fbc->face);

        ok_face = p4 != NULL;

        if (ok_face) {
            fbc->face      = p4;
            fbc->face_cpty = alloc_sz;
        }
    }

    return ok_cond && ok_face;
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
        ok = initialise_structure(fbc);

        ok = ok && expand_tables(nbc, nbc, fbc);

        if (! ok) {
            flow_conditions_destroy(fbc);

            fbc = NULL;
        } else {
            fbc->cond_pos[0] = 0;
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
        free(fbc->face    );
        free(fbc->cond_pos);
        free(fbc->value   );
        free(fbc->type    );
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
    return flow_conditions_append_multi(type, 1, &face, value, fbc);
}


/* ---------------------------------------------------------------------- */
/* Append a new boundary condition that affects multiple interfaces.
 *
 * Return one (1) if successful, and zero (0) otherwise. */
/* ---------------------------------------------------------------------- */
int
flow_conditions_append_multi(enum FlowBCType                type  ,
                             size_t                         nfaces,
                             const int                     *faces ,
                             double                         value ,
                             struct FlowBoundaryConditions *fbc   )
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t nbc;

    nbc = fbc->nbc;

    ok  = expand_tables(nbc + 1, fbc->cond_pos[ nbc ] + nfaces, fbc);

    if (ok) {
        memcpy(fbc->face + fbc->cond_pos[ nbc ],
               faces, nfaces * sizeof *faces);

        fbc->type [ nbc ] = type;
        fbc->value[ nbc ] = value;

        fbc->cond_pos[ nbc + 1 ] = fbc->cond_pos[ nbc ] + nfaces;

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
