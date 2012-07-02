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

/**
 * \file
 * Interfaces and data structures to assemble a system of simultaneous linear
 * equations discretising a flow problem that is either incompressible or
 * features rock compressibility using the two-point flux approximation method.
 *
 * Includes support for reconstructing the Darcy flux field as well as well
 * connection fluxes.
 */

#include <opm/core/grid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ifs_tpfa_impl;
struct CSRMatrix;
struct FlowBoundaryConditions;
struct Wells;

/**
 * Main data structure presenting a view of an assembled system of simultaneous
 * linear equations which may be solved using external software.
 */
struct ifs_tpfa_data {
    struct CSRMatrix     *A;      /**< Coefficient matrix */
    double               *b;      /**< Right-hand side */
    double               *x;      /**< Solution */

    struct ifs_tpfa_impl *pimpl;  /**< Internal management structure */
};

/**
 * Solution variables.
 */
struct ifs_tpfa_solution {
    double *cell_press;  /**< Cell pressures */
    double *face_flux ;  /**< Interface fluxes */

    double *well_press;  /**< Bottom-hole pressures for each well */
    double *well_flux ;  /**< Well connection total fluxes */
};

/**
 * Driving forces pertaining to a particular model setup.
 */
struct ifs_tpfa_forces {
    const double                        *src;  /**< Explicit source terms */
    const struct FlowBoundaryConditions *bc ;  /**< Boundary conditions */

    const struct Wells *W     ;  /**< Well topology */
    const double       *totmob;  /**< Total mobility in each cell */
    const double       *wdp   ;  /**< Gravity adjustment at each perforation */
};


/**
 * Allocate TPFA management structure capable of assembling a system of
 * simultaneous linear equations corresponding to a particular grid and well
 * configuration.
 *
 * @param[in] G Grid.
 * @param[in] W Well topology.
 * @return Fully formed TPFA management structure if successful, @c NULL in case
 * of allocation failure.
 */
struct ifs_tpfa_data *
ifs_tpfa_construct(struct UnstructuredGrid *G,
                   struct Wells            *W);


/**
 *
 * @param[in]     G
 * @param[in]     F
 * @param[in]     trans
 * @param[in]     gpress
 * @param[in,out] h
 * @return
 */
int
ifs_tpfa_assemble(struct UnstructuredGrid      *G     ,
                  const struct ifs_tpfa_forces *F     ,
                  const double                 *trans ,
                  const double                 *gpress,
                  struct ifs_tpfa_data         *h     );

int
ifs_tpfa_assemble_comprock(struct UnstructuredGrid      *G        ,
                           const struct ifs_tpfa_forces *F        ,
                           const double                 *trans    ,
                           const double                 *gpress   ,
                           const double                 *porevol  ,
                           const double                 *rock_comp,
                           const double                  dt       ,
                           const double                 *pressure ,
                           struct ifs_tpfa_data         *h        );
int
ifs_tpfa_assemble_comprock_increment(struct UnstructuredGrid      *G        ,
				     const struct ifs_tpfa_forces *F        ,
				     const double                 *trans    ,
				     const double                 *gpress   ,
				     const double                 *porevol  ,
				     const double                 *rock_comp,
				     const double                  dt       ,
				     const double                 *prev_pressure ,
				     const double                 *initial_porevolume,
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
