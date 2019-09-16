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
 * \ingroup RichardsModel
 *
 * \brief Contains the property declarations for the Richards model.
 */
#ifndef EWOMS_RICHARDS_PROPERTIES_HH
#define EWOMS_RICHARDS_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>

// \{
BEGIN_PROPERTIES

//! The fluid system used for the problem
NEW_PROP_TAG(FluidSystem);

//! The fluid used as the wetting phase (by default, we set the fluid
//! system to the immiscible one, which requires this property.)
NEW_PROP_TAG(WettingFluid);

//! The fluid used as the non-wetting phase (by default, we set the
//! fluid system to the immiscible one, which requires this property.)
NEW_PROP_TAG(NonWettingFluid);

//! Index of the fluid which represents the wetting phase
NEW_PROP_TAG(LiquidPhaseIndex);

//! Index of the fluid which represents the non-wetting phase
NEW_PROP_TAG(GasPhaseIndex);

//! Index of the component which constitutes the liquid
NEW_PROP_TAG(LiquidComponentIndex);

//! Index of the component which constitutes the gas
NEW_PROP_TAG(GasComponentIndex);

// \}

END_PROPERTIES

#endif
