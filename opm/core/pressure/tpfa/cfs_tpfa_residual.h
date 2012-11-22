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

#include <opm/core/grid.h>
#include <opm/core/newwells.h>

#include <opm/core/pressure/tpfa/compr_source.h>

/**
 * \file
 * Public interface to assembler for (compressible) discrete pressure system
 * based on two-point flux approximation method.  The assembler implements a
 * residual formulation that for a single cell \f$i\f$ reads
 * \f[
 * \mathit{pv}_i\cdot (1 - \sum_\alpha A_i^{-1}(p^{n+1},z^n)z^n) +
 * \Delta t\sum_\alpha A_i^{-1}(p^{n+1},z^n) \Big(\sum_j A_{ij}v_{ij}^{n+1} -
 * q_i^{n+1}\Big) = 0
 * \f]
 * in which \f$\mathit{pv}_i\f$ is the (constant or pressure-dependent)
 * pore-volume of cell \f$i\f$.  Moreover, \f$\Delta t\f$ is the time step size,
 * \f$n\f$ denotes the time level, \f$A\f$ is the pressure and mass-dependent
 * fluid matrix that converts phase volumes at reservoir conditions into
 * component volumes at surface conditions and \f$v_{ij}\f$ is the vector of
 * outward (with respect to cell \f$i\f$) phase fluxes across the \f$ij\f$
 * cell-interface.
 *
 * This module's usage model is intended to be
 * -# Construct assembler
 * -# for (each time step)
 *    -# while (pressure not converged)
 *       -# Assemble appropriate (Jacobian) system of linear equations
 *       -# Solve Jacobian system to derive pressure increments
 *       -# Include increments into current state
 *       -# Check convergence
 *    -# Derive fluxes and, optionally, interface pressures
 *    -# Solve transport by some means
 * -# Destroy assembler
 */

#ifdef __cplusplus
extern "C" {
#endif

struct cfs_tpfa_res_impl;
struct CSRMatrix;
struct compr_quantities_gen;

/**
 * Type encapsulating well topology and completion data (e.g., phase mobilities
 * per connection (perforation)).
 */
struct cfs_tpfa_res_wells {
    /**
     * All wells pertaining to a particular linear system assembly.
     * Must include current controls/targets and complete well topology.
     */
    struct Wells          *W   ;

    /**
     * Completion data describing the fluid state at the current time level.
     */
    struct CompletionData *data;
};


/**
 * Type encapsulating all driving forces affecting the discrete pressure system.
 */
struct cfs_tpfa_res_forces {
    struct cfs_tpfa_res_wells *wells;  /**< Wells */
    struct compr_src          *src  ;  /**< Explicit source terms */
};


/**
 * Result structure that presents the fully assembled system of linear
 * equations, linearised around the current pressure point.
 */
struct cfs_tpfa_res_data {
    struct CSRMatrix         *J; /**< Jacobian matrix */
    double                   *F; /**< Residual vector (right-hand side) */

    struct cfs_tpfa_res_impl *pimpl; /**< Internal management structure */
};


/**
 * Construct assembler for system of linear equations.
 *
 * @param[in] G        Grid
 * @param[in] wells    Well description.  @c NULL in case of no wells.
 *                     For backwards compatibility, the constructor also
 *                     interprets <CODE>(wells != NULL) &&
 *                     (wells->W == NULL)</CODE> as "no wells present", but new
 *                     code should use <CODE>wells == NULL</CODE> to signify
 *                     "no wells".
 * @param[in] nphases  Number of active fluid phases in this simulation run.
 *                     Needed to correctly size various internal work arrays.
 * @return Fully formed assembler structure suitable for forming systems of
 * linear equations using, e.g., function cfs_tpfa_res_assemble().  @c NULL in
 * case of allocation failure.  Must be destroyed using function
 * cfs_tpfa_res_destroy().
 */
struct cfs_tpfa_res_data *
cfs_tpfa_res_construct(struct UnstructuredGrid   *G      ,
                       struct cfs_tpfa_res_wells *wells  ,
                       int                        nphases);


/**
 * Destroy assembler for system of linear equations.
 *
 * Disposes of all resources acquired in a previous call to construction
 * function cfs_tpfa_res_construct().  Note that the statement
 * <CODE>
 *   cfs_tpfa_res_destroy(NULL)
 * </CODE>
 * is supported and benign (i.e., behaves like <CODE>free(NULL)</CODE>).
 *
 * @param[in,out] h On input - assembler obtained through a previous call to
 * construction function cfs_tpfa_res_construct().  On output - invalid pointer.
 */
void
cfs_tpfa_res_destroy(struct cfs_tpfa_res_data *h);


/**
 * Assemble system of linear equations by linearising the residual around the
 * current pressure point.  Assume incompressible rock (i.e., that the
 * pore-volume is independent of pressure).
 *
 * The fully assembled system is presented in <CODE>h->J</CODE> and
 * <CODE>h->F</CODE> and must be solved separately using external software.
 *
 * @param[in]     G         Grid.
 * @param[in]     dt        Time step size \f$\Delta t\f$.
 * @param[in]     forces    Driving forces.
 * @param[in]     zc        Component volumes, per pore-volume, at surface
 *                          conditions for all components in all cells stored
 *                          consecutively per cell.  Array of size
 *                          <CODE>G->number_of_cells * cq->nphases</CODE>.
 * @param[in]     cq        Compressible quantities describing the current fluid
 *                          state.  Fields @c Ac, @c dAc, @c Af, and
 *                          @c phasemobf must be valid.
 * @param[in]     trans     Background transmissibilities as defined by function
 *                          tpfa_trans_compute().
 * @param[in]     gravcap_f Discrete gravity and capillary forces.
 * @param[in]     cpress    Cell pressures.  One scalar value per grid cell.
 * @param[in]     wpress    Well (bottom-hole) pressures.  One scalar value per
 *                          well.  @c NULL in case of no wells.
 * @param[in]     porevol   Pore-volumes.  One (positive) scalar value for each
 *                          grid cell.
 * @param[in,out] h         On input-a valid (non-@c NULL) assembler obtained
 *                          from a previous call to constructor function
 *                          cfs_tpfa_res_construct().  On output-valid assembler
 *                          that includes the new system of linear equations in
 *                          its @c J and @c F fields.
 *
 * @return 1 if the assembled matrix was adjusted to remove a singularity.  This
 * happens if all fluids are incompressible and there are no pressure conditions
 * on wells or boundaries.  Otherwise return 0.
 */
int
cfs_tpfa_res_assemble(struct UnstructuredGrid     *G,
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


/**
 * Assemble system of linear equations by linearising the residual around the
 * current pressure point.  Assume compressible rock (i.e., that the pore-volume
 * depends on pressure).
 *
 * The fully assembled system is presented in <CODE>h->J</CODE> and
 * <CODE>h->F</CODE> and must be solved separately using external software.
 *
 * @param[in]     G         Grid.
 * @param[in]     dt        Time step size \f$\Delta t\f$.
 * @param[in]     forces    Driving forces.
 * @param[in]     zc        Component volumes, per pore-volume, at surface
 *                          conditions for all components in all cells stored
 *                          consecutively per cell.  Array of size
 *                          <CODE>G->number_of_cells * cq->nphases</CODE>.
 * @param[in]     cq        Compressible quantities describing the current fluid
 *                          state.  Fields @c Ac, @c dAc, @c Af, and
 *                          @c phasemobf must be valid.
 * @param[in]     trans     Background transmissibilities as defined by function
 *                          tpfa_trans_compute().
 * @param[in]     gravcap_f Discrete gravity and capillary forces.
 * @param[in]     cpress    Cell pressures.  One scalar value per grid cell.
 * @param[in]     wpress    Well (bottom-hole) pressures.  One scalar value per
 *                          well.  @c NULL in case of no wells.
 * @param[in]     porevol   Pore-volumes.  One (positive) scalar value for each
 *                          grid cell.
 * @param[in]     porevol0  Pore-volumes at start of time step (i.e., at time
 *                          level \f$n\f$).  One (positive) scalar value for
 *                          each grid cell.
 * @param[in]     rock_comp Rock compressibility.  One non-negative scalar for
 *                          each grid cell.
 * @param[in,out] h         On input-a valid (non-@c NULL) assembler obtained
 *                          from a previous call to constructor function
 *                          cfs_tpfa_res_construct().  On output-valid assembler
 *                          that includes the new system of linear equations in
 *                          its @c J and @c F fields.
 *
 * @return 1 if the assembled matrix was adjusted to remove a singularity.  This
 * happens if all fluids are incompressible, the rock is incompressible, and
 * there are no pressure conditions on wells or boundaries.  Otherwise return 0.
 */
int
cfs_tpfa_res_comprock_assemble(
                      struct UnstructuredGrid     *G,
                      double                       dt,
                      struct cfs_tpfa_res_forces  *forces,
                      const double                *zc,
                      struct compr_quantities_gen *cq,
                      const double                *trans,
                      const double                *gravcap_f,
                      const double                *cpress,
                      const double                *wpress,
                      const double                *porevol,
                      const double                *porevol0,
                      const double                *rock_comp,
                      struct cfs_tpfa_res_data    *h);


/**
 * Derive interface (total) Darcy fluxes from (converged) pressure solution.
 *
 * @param[in]  G         Grid
 * @param[in]  forces    Driving forces.  Must correspond to those used when
 *                       forming the system of linear equations, e.g., in the
 *                       call to function cfs_tpfa_res_assemble().
 * @param[in]  np        Number of fluid phases (and components).
 * @param[in]  trans     Background transmissibilities as defined by function
 *                       tpfa_trans_compute().  Must correspond to equally
 *                       named parameter of the assembly functions.
 * @param[in]  pmobc     Phase mobilities stored consecutively per cell with
 *                       phase index cycling the most rapidly.  Array of size
 *                       <CODE>G->number_of_cells * np</CODE>.
 * @param[in]  pmobf     Phase mobilities stored consecutively per interface
 *                       with phase index cycling the most rapidly.  Array of
 *                       size <CODE>G->number_of_faces * np</CODE>.
 * @param[in]  gravcap_f Discrete gravity and capillary forces.
 * @param[in]  cpress    Cell pressure values.  One (positive) scalar for each
 *                       grid cell.
 * @param[in]  wpress    Well (bottom-hole) pressure values.  One (positive)
 *                       scalar value for each well.  @c NULL in case of no
 *                       wells.
 * @param[out] fflux     Total Darcy interface fluxes.  One scalar value for
 *                       each interface (inter-cell connection).  Array of size
 *                       <CODE>G->number_of_faces</CODE>.
 * @param[out] wflux     Total Darcy well connection fluxes.  One scalar value
 *                       for each well connection (perforation).  Array of size
 *                       <CODE>forces->wells->W->well_connpos
 *                       [forces->wells->W->number_of_wells]</CODE>.
 */
void
cfs_tpfa_res_flux(struct UnstructuredGrid    *G        ,
                  struct cfs_tpfa_res_forces *forces   ,
                  int                         np       ,
                  const double               *trans    ,
                  const double               *pmobc    ,
                  const double               *pmobf    ,
                  const double               *gravcap_f,
                  const double               *cpress   ,
                  const double               *wpress   ,
                  double                     *fflux    ,
                  double                     *wflux    );


/**
 * Derive interface pressures from (converged) pressure solution.
 *
 * @param[in]  G         Grid
 * @param[in]  np        Number of fluid phases (and components).
 * @param[in]  htrans    Background one-sided ("half") transmissibilities as
 *                       defined by function tpfa_htrans_compute().
 * @param[in]  pmobf     Phase mobilities stored consecutively per interface
 *                       with phase index cycling the most rapidly.  Array of
 *                       size <CODE>G->number_of_faces * np</CODE>.
 * @param[in]  gravcap_f Discrete gravity and capillary forces.
 * @param[in]  h         System assembler.  Must correspond to the assembler
 *                       state used to form the final system of linear equations
 *                       from which the converged pressure solution was derived.
 * @param[in]  cpress    Cell pressure values.  One (positive) scalar for each
 *                       grid cell.
 * @param[in]  fflux     Total Darcy interface fluxes.  One scalar value for
 *                       each interface (inter-cell connection).  Array of size
 *                       <CODE>G->number_of_faces</CODE>.  Typically computed
 *                       using function cfs_tpfa_res_flux().
 * @param[out] fpress    Interface pressure values.  One (positive) scalar for
 *                       each interface.  Array of size
 *                       <CODE>G->number_of_faces</CODE>.
 */
void
cfs_tpfa_res_fpress(struct UnstructuredGrid  *G,
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
cfs_tpfa_retrieve_masstrans(struct UnstructuredGrid *G,
                            int                      np,
                            struct cfs_tpfa_data    *h,
                            double                  *masstrans_f);

void
cfs_tpfa_retrieve_gravtrans(struct UnstructuredGrid *G,
                            int                      np,
                            struct cfs_tpfa_data    *h,
                            double                  *gravtrans_f);

double
cfs_tpfa_impes_maxtime(struct UnstructuredGrid *G,
                       struct compr_quantities *cq,
                       const double            *trans,
                       const double            *porevol,
                       struct cfs_tpfa_data    *h,
                       const double            *dpmobf,
                       const double            *surf_dens,
                       const double            *gravity);

void
cfs_tpfa_expl_mass_transport(struct UnstructuredGrid *G,
                             well_t                  *W,
                             struct completion_data  *wdata,
                             int                      np,
                             double                   dt,
                             const double            *porevol,
                             struct cfs_tpfa_data    *h,
                             double                  *surf_vol);

#endif

#ifdef __cplusplus
}
#endif

#endif  /* OPM_CFS_TPFA_HEADER_INCLUDED */
