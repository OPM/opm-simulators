/*===========================================================================
//
// File: compr_source.h
//
// Created: 2011-10-19 19:14:30+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_COMPR_SOURCE_H_HEADER
#define OPM_COMPR_SOURCE_H_HEADER

/**
 * \file
 * Data structures and support routines needed to represent explicit,
 * compressible source terms.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Collection of explicit, compressible source terms.
 */
struct compr_src {
    /**
     * Number of source terms.
     */
    int nsrc;

    /**
     * Source term capacity.  Client code should treat this member as read-only.
     * The field is used in internal memory management.
     */
    int cpty;

    /**
     * Number of fluid phases.
     */
    int nphases;

    /**
     * Cells influenced by explicit source terms.  Array of size @c cpty, the
     * @c nsrc first elements (only) of which are valid.
     */
    int *cell;

    /**
     * Total Darcy rate of inflow (measured at reservoir conditions) of each
     * individual source term.  Sign convention: Positive rate into reservoir
     * (i.e., injection) and negative rate out of reservoir (production).
     * Array of size @c cpty, the @c nsrc first elements (only) of which are
     * valid.
     */
    double *flux;

    /**
     * Injection composition for all explicit source terms.  Not referenced for
     * production sources (i.e., those terms for which <CODE>->flux[]</CODE> is
     * negative).  Array of size <CODE>nphases * cpty</CODE>, the
     * <CODE>nphases * nsrc</CODE> of which (only) are valid.
     */
    double *saturation;
};


/**
 * Create a management structure that is, initially, capable of storing a
 * specified number of sources defined by a particular number a fluid phases.
 *
 * @param[in] np    Number of fluid phases.  Must be positive to actually
 *                  allocate any sources.
 * @param[in] nsrc  Initial management capacity.  If positive, attempt to
 *                  allocate that number of source terms.  Otherwise, the
 *                  initial capacity is treated as (and set to) zero.
 * @return Fully formed management structure if <CODE>np > 0</CODE> and
 * allocation success.  @c NULL otherwise.  The resources must be released using
 * destructor function compr_src_deallocate().
 */
struct compr_src *
compr_src_allocate(int np, int nsrc);


/**
 * Release memory resources acquired in a previous call to constructor function
 * compr_src_allocate() and, possibly, source term insertion function
 * append_compr_source_term().
 *
 * @param[in,out] src On input - source term management structure obtained
 * through a previous call to construction function compr_src_allocate().
 * On output - invalid pointer.
 */
void
compr_src_deallocate(struct compr_src *src);


/**
 * Insert a new explicit source term into an existing collection.
 *
 * @param[in]     c   Cell influenced by this source term.
 * @param[in]     np  Number of fluid phases.  Used for consistency checking
 *                    only.
 * @param[in]     v   Source term total reservoir volume Darcy flux.  Positive
 *                    if the source term is to be interpreted as an injection
 *                    source and negative otherwise.
 * @param[in]     sat Injection composition for this source term.  Array of size
 *                    @c np.  Copied to internal storage, but the actual numbers
 *                    are not inspected unless <CODE>v > 0.0</CODE>.
 * @param[in,out] src On input - source term management structure obtained
 * through a previous call to construction function compr_src_allocate() and,
 * possibly, another call to function append_compr_source_term().  On output -
 * source term collection that includes the new source term if successful and
 * unchanged if not.
 *
 * @return One (1, true) if successful (i.e., if the new source term could be
 * included in the existing collection) and zero (0, false) if not.
 */
int
append_compr_source_term(int               c  ,
                         int               np ,
                         double            v  ,
                         const double     *sat,
                         struct compr_src *src);


/**
 * Empty source term collection while maintaining existing capacity.
 *
 * @param[in,out] src On input - an existing source term collection with a
 * given number of sources and source capacity.  On output - an empty source
 * term collection (i.e., <CODE>src->nsrc == 0</CODE>) with an unchanged
 * capacity.
 */
void
clear_compr_source_term(struct compr_src *src);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_SOURCE_H_HEADER */
