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

#ifdef __cplusplus
extern "C" {
#endif

struct cfs_tpfa_res_impl;
struct CSRMatrix;
struct compr_quantities_gen;
struct compr_src;
struct WellCompletions;
struct WellControls;
struct completion_data;

struct cfs_tpfa_res_wells {
    struct WellCompletions *conn;
    struct WellControls    *ctrl;
    struct completion_data *data;
};

struct cfs_tpfa_res_forces {
    struct cfs_tpfa_res_wells *W  ;
    struct compr_src          *src;
};

struct cfs_tpfa_res_data {
    struct CSRMatrix         *J;
    double                   *F;

    struct cfs_tpfa_res_impl *pimpl;
};


struct cfs_tpfa_res_data *
cfs_tpfa_res_construct(grid_t                 *G      ,
                       struct WellCompletions *wconn  ,
                       int                     nphases);

void
cfs_tpfa_res_destroy(struct cfs_tpfa_res_data *h);

void
cfs_tpfa_res_assemble(grid_t                      *G,
                      double                       dt,
                      struct cfs_tpfa_res_forces  *forces,
                      const double                *zc,
                      struct compr_quantities_gen *cq,
                      const double                *trans,
                      const double                *gravcap_f,
                      const double                *cpress,
                      const double                *wpress,
                      const double                *porevol,
                      struct cfs_tpfa_res_data    *h);

void
cfs_tpfa_res_flux(grid_t       *G,
                  int           np,
                  const double *trans,
                  const double *pmobf,
                  const double *gravcap_f,
                  const double *cpress,
                  double       *fflux);

void
cfs_tpfa_res_fpress(grid_t                   *G,
                    flowbc_t                 *bc,
                    int                       np,
                    const double             *htrans,
                    const double             *pmobf,
                    const double             *gravcap_f,
                    struct cfs_tpfa_res_data *h,
                    const double             *cpress,
                    const double             *fflux,
                    double                   *fpress);

#if 0
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

double
cfs_tpfa_impes_maxtime(grid_t                  *G,
                       struct compr_quantities *cq,
                       const double            *trans,
                       const double            *porevol,
                       struct cfs_tpfa_data    *h,
                       const double            *dpmobf,
                       const double            *surf_dens,
                       const double            *gravity);

void
cfs_tpfa_expl_mass_transport(grid_t                 *G,
                             well_t                 *W,
                             struct completion_data *wdata,
                             int                     np,
                             double                  dt,
                             const double           *porevol,
                             struct cfs_tpfa_data   *h,
                             double                 *surf_vol);

#endif

#ifdef __cplusplus
}
#endif

#endif  /* OPM_CFS_TPFA_HEADER_INCLUDED */
