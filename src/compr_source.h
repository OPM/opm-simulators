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

#ifdef __cplusplus
extern "C" {
#endif

struct compr_src {
    int     nsrc;
    int     cpty;

    int     nphases;

    int    *cell;
    double *flux;
    double *saturation;
};


struct compr_src *
compr_src_allocate(int np, int nsrc);

void
compr_src_deallocate(struct compr_src *src);

int
append_compr_source_term(int               c  ,
                         int               np ,
                         double            v  ,
                         const double     *sat,
                         struct compr_src *src);

void
clear_compr_source_term(struct compr_src *src);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_COMPR_SOURCE_H_HEADER */
