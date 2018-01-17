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

#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_SOURCE_H_HEADER */
