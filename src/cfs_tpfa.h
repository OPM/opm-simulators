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

#ifndef OPM_CFS_TPFA_HEADER_INCLUDED
#define OPM_CFS_TPFA_HEADER_INCLUDED

#include "grid.h"
#include "flow_bc.h"
#include "well.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cfs_tpfa_impl;
struct CSRMatrix;
struct compr_quantities;

struct cfs_tpfa_data {
    struct CSRMatrix     *A;
    double               *b;
    double               *x;

    struct cfs_tpfa_impl *pimpl;
};


struct cfs_tpfa_data *
cfs_tpfa_construct(grid_t *G, well_t *W, int nphases);

void
cfs_tpfa_assemble(grid_t                  *G,
                  double                   dt,
                  well_t                  *W,
                  flowbc_t                *bc,
                  const double            *src,
                  struct compr_quantities *cq,
                  const double            *trans,
                  const double            *gravcap_f,
                  well_control_t          *wctrl,
                  struct completion_data  *wdata,
                  const double            *cpress0,
                  const double            *porevol,
                  struct cfs_tpfa_data    *h);

void
cfs_tpfa_press_flux(grid_t                 *G,
                    flowbc_t               *bc,
                    well_t                 *W,
                    int                     np,
                    const double           *trans,
                    const double           *pmobf,
                    const double           *gravcap_f,
                    struct completion_data *wdata,
                    struct cfs_tpfa_data   *h,
                    double                 *cpress,
                    double                 *fflux,
                    double                 *wpress,
                    double                 *wflux);

void
cfs_tpfa_fpress(grid_t       *G,
                flowbc_t     *bc,
                int           np,
                const double *htrans,
                const double *pmobf,
                const double *gravcap_f,
                const double *cpress,
                const double *fflux,
                double       *fpress);

void
cfs_tpfa_retrieve_masstrans(grid_t               *G,
                            int                   np,
                            struct cfs_tpfa_data *h,
                            double               *masstrans_f);

void
cfs_tpfa_retrieve_gravtrans(grid_t               *G,
                            int                   np,
                            struct cfs_tpfa_data *h,
                            double               *gravtrans_f);

void
cfs_tpfa_expl_mass_transport(grid_t               *G,
                             well_t               *W,
                             int                   np,
                             double                dt,
                             const double         *porevol,
                             struct cfs_tpfa_data *h,
                             double               *surf_vol);

void
cfs_tpfa_destroy(struct cfs_tpfa_data *h);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_CFS_TPFA_HEADER_INCLUDED */
