// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This test is identical to the simulation of the lens problem that uses the
 *        element centered finite volume discretization in conjunction with automatic
 *        differentiation (lens_immiscible_ecfv_ad).
 *
 * The only difference is that it uses multiple compile units in order to ensure that
 * eWoms code can be used within libraries that use the same type tag within multiple
 * compile units. This file represents the first compile unit and just defines an startup
 * function for the simulator.
 */
#include "config.h"

#include "lens_immiscible_ecfv_ad.hh"

#include <opm/models/utils/start.hh>

// fake forward declaration to prevent esoteric compiler warning
int mainCU1(int argc, char **argv);

int mainCU1(int argc, char **argv)
{
    using ProblemTypeTag = Opm::Properties::TTag::LensProblemEcfvAd;
    return Opm::start<ProblemTypeTag>(argc, argv);
}
