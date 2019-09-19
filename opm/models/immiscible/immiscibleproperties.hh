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
 * \ingroup ImmiscibleModel
 *
 * \brief Defines the properties required for the immiscible
 *        multi-phase model.
 */
#ifndef EWOMS_IMMISCIBLE_PROPERTIES_HH
#define EWOMS_IMMISCIBLE_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/io/vtkenergymodule.hh>

BEGIN_PROPERTIES

//!The fluid systems including the information about the phases
NEW_PROP_TAG(FluidSystem);
//! Specify whether energy should be considered as a conservation quantity or not
NEW_PROP_TAG(EnableEnergy);

// these properties only make sense for the ImmiscibleTwoPhase type tag

//! The wetting phase for two-phase models
NEW_PROP_TAG(WettingPhase);
//! The non-wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase);

// these properties only make sense for the ImmiscibleSinglePhase type tag

//! The fluid used by the model
NEW_PROP_TAG(Fluid);


END_PROPERTIES

#endif
