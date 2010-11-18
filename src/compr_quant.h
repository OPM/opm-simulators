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

#include "grid.h"

#ifdef __cplusplus
extern "C" {
#endif

struct compr_quantities {
    int     nphases;

    double *totcompr;
    double *voldiscr;
    double *Ac;                 /* RB^{-1} per cell */
    double *Af;                 /* RB^{-1} per face */
    double *phasemobf;          /* Phase mobility per face */
};

void
compr_flux_term(grid_t       *G,
                const double *fflux,
                const double *zeta,
                double       *Biv);

void
compr_accum_term(size_t        nc,
                 double        dt,
                 const double *porevol,
                 const double *totcompr,
                 double       *P);

void
compr_src_add_press_accum(size_t        nc,
                          const double *p0,
                          const double *P,
                          double       *src);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_QUANT_HEADER_INCLUDED */
