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

#include <stdlib.h>

#include "flow_bc.h"


/* Create structure to hold flow boundary conditions for 'nf' faces.
 *
 * Return fully allocated structure or NULL in case of allocation
 * failure. */
/* ---------------------------------------------------------------------- */
flowbc_t *
allocate_flowbc(size_t nf)
/* ---------------------------------------------------------------------- */
{
    size_t    i;
    flowbc_t *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->type  = malloc(nf * sizeof *new->type);
        new->bcval = malloc(nf * sizeof *new->bcval);

        if ((new->type == NULL) || (new->bcval == NULL)) {
            deallocate_flowbc(new);
            new = NULL;
        } else {
            for (i = 0; i < nf; i++) {
                new->type [i] = UNSET;
                new->bcval[i] = 0.0;
            }
        }
    }

    return new;
}


/* Release memory resources for dynamically allocated flow boundary
 * condition structure. */
/* ---------------------------------------------------------------------- */
void
deallocate_flowbc(flowbc_t *fbc)
/* ---------------------------------------------------------------------- */
{
    if (fbc != NULL) {
        free(fbc->bcval);
        free(fbc->type );
    }

    free(fbc);
}
