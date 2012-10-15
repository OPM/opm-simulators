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

/**
 * \file
 * Module defining derived fluid quantities needed to discretise compressible
 * and miscible pressure (flow) problems.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Aggregate structure that represents an atomic view of the current fluid
 * state.  These quantities are used directly in the cfs_tpfa_residual module to
 * discretise a particular, linearised flow problem.
 */
struct compr_quantities_gen {
    /**
     * Number of fluid phases.  The pressure solvers also assume that the number
     * of fluid components (at surface conditions) equals the number of fluid
     * phases.
     */
    int nphases;

    /**
     * Pressure and mass-dependent fluid matrix that converts phase volumes at
     * reservoir conditions into component volumes at surface conditions.  Obeys
     * the defining formula
     * \f[
     * A = RB^{-1}
     * \f]
     * in which \f$R\f$ denotes the miscibility ratios (i.e., the dissolved
     * gas-oil ratio, \f$R_s\f$ and the evaporated oil-gas ratio, \f$R_v\f$)
     * collected as a matrix and \f$B\f$ is the diagonal matrix of
     * formation-volume expansion factors.  The function is sampled in each grid
     * cell.  Array of size <CODE>nphases * nphases * nc</CODE>.
     */
    double *Ac;

    /**
     * Derivative of \f$A\f$ with respect to pressure,
     * \f[
     * \frac{\partial A}{\partial p} = \frac{\partial R}{\partial p}B^{-1} +
     * R\frac{\partial B^{-1}}{\partial p} = (\frac{\partial R}{\partial p} -
     * A\frac{\partial B}{\partial p})B^{-1}
     * \f]
     * sampled in each grid cell.  Array of size
     * <CODE>nphases * nphases * nc</CODE>.
     */
    double *dAc;

    /**
     * Fluid matrix sampled at each interface.  Possibly as a result of an
     * upwind calculation.  Array of size <CODE>nphases * nphases * nf</CODE>.
     */
    double *Af;

    /**
     * Phase mobility per interface.  Possibly defined through an upwind
     * calculation.  Array of size <CODE>nphases * nf</CODE>.
     */
    double *phasemobf;

    /**
     * Deceptively named "volume-discrepancy" term.  Array of size @c nc.
     * Unused by the cfs_tpfa_residual module.
     */
    double *voldiscr;
};


/**
 * Create management structure capable of storing derived fluid quantities
 * pertaining to a reservoir model of a particular size.
 *
 * The resources acquired using function compr_quantities_gen_allocate() must
 * be released using the destructor function compr_quantities_gen_deallocate().
 *
 * @param[in] nc Number of grid cells
 * @param[in] nf Number of grid interfaces
 * @param[in] np Number of fluid phases (and components)
 *
 * @return Fully formed management structure.  Returns @c NULL in case of
 * allocation failure.
 */
struct compr_quantities_gen *
compr_quantities_gen_allocate(size_t nc, size_t nf, int np);


/**
 * Release resources acquired in a previous call to constructor function
 * compr_quantities_gen_allocate().
 *
 * Note that
 * <CODE>
 * compr_quantities_gen_deallocate(NULL)
 * </CODE>
 * is supported and benign (i.e., this statement behaves as
 * <CODE>free(NULL)</CODE>).
 *
 * @param[in,out] cq On input - compressible quantity management structure
 * obtained through a previous call to construction function
 * compr_quantities_gen_allocate().  On output - invalid pointer.
 */
void
compr_quantities_gen_deallocate(struct compr_quantities_gen *cq);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_QUANT_HEADER_INCLUDED */
