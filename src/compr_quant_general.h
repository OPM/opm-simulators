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

#ifndef OPM_COMPR_QUANT_HEADER_INCLUDED
#define OPM_COMPR_QUANT_HEADER_INCLUDED

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct compr_quantities_gen {
    int     nphases;      /* Number of phases/components */

    double *Ac;           /* RB^{-1} per cell */
    double *dAc;          /* d/dp (RB^{-1}) per cell */
    double *Af;           /* RB^{-1} per face */
    double *phasemobf;    /* Phase mobility per face */
    double *voldiscr;     /* Volume discrepancy per cell */
};

struct compr_quantities_gen *
compr_quantities_gen_allocate(size_t nc, size_t nf, int np);

void
compr_quantities_gen_deallocate(struct compr_quantities_gen *cq);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_QUANT_HEADER_INCLUDED */
