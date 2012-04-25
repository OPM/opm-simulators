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

#ifndef OPM_IFS_TPFA_HEADER_INCLUDED
#define OPM_IFS_TPFA_HEADER_INCLUDED

#include <opm/core/grid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ifs_tpfa_impl;
struct CSRMatrix;
struct FlowBoundaryConditions;
struct Wells;

struct ifs_tpfa_data {
    struct CSRMatrix     *A;
    double               *b;
    double               *x;

    struct ifs_tpfa_impl *pimpl;
};

struct ifs_tpfa_solution {
    double *cell_press;
    double *face_flux ;

    double *well_press;         /* BHP */
    double *well_flux ;         /* Perforation (total) fluxes */
};

struct ifs_tpfa_forces {
    const double                        *src;
    const struct FlowBoundaryConditions *bc ;

    const struct Wells *W     ;
    const double       *totmob;
    const double       *wdp   ;
};


struct ifs_tpfa_data *
ifs_tpfa_construct(struct UnstructuredGrid *G,
                   struct Wells            *W);

void
ifs_tpfa_assemble(struct UnstructuredGrid      *G     ,
                  const struct ifs_tpfa_forces *F     ,
                  const double                 *trans ,
                  const double                 *gpress,
                  struct ifs_tpfa_data         *h     );

void
ifs_tpfa_assemble_comprock(struct UnstructuredGrid      *G        ,
                           const struct ifs_tpfa_forces *F        ,
                           const double                 *trans    ,
                           const double                 *gpress   ,
                           const double                 *porevol  ,
                           const double                 *rock_comp,
                           const double                  dt       ,
                           const double                 *pressure ,
                           struct ifs_tpfa_data         *h        );
void
ifs_tpfa_press_flux(struct UnstructuredGrid      *G    ,
                    const struct ifs_tpfa_forces *F    ,
                    const double                 *trans,
                    struct ifs_tpfa_data         *h    ,
                    struct ifs_tpfa_solution     *soln );

void
ifs_tpfa_destroy(struct ifs_tpfa_data *h);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_IFS_TPFA_HEADER_INCLUDED */
