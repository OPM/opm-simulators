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
namespace Opm::Properties {

//! The fluid used as the wetting phase (by default, we set the fluid
//! system to the immiscible one, which requires this property.)
template<class TypeTag, class MyTypeTag>
struct WettingFluid { using type = UndefinedProperty; };

//! The fluid used as the non-wetting phase (by default, we set the
//! fluid system to the immiscible one, which requires this property.)
template<class TypeTag, class MyTypeTag>
struct NonWettingFluid { using type = UndefinedProperty; };

//! Index of the fluid which represents the wetting phase
template<class TypeTag, class MyTypeTag>
struct LiquidPhaseIndex { using type = UndefinedProperty; };

//! Index of the fluid which represents the non-wetting phase
template<class TypeTag, class MyTypeTag>
struct GasPhaseIndex { using type = UndefinedProperty; };

//! Index of the component which constitutes the liquid
template<class TypeTag, class MyTypeTag>
struct LiquidComponentIndex { using type = UndefinedProperty; };

//! Index of the component which constitutes the gas
template<class TypeTag, class MyTypeTag>
struct GasComponentIndex { using type = UndefinedProperty; };

// \}

} // namespace Opm::Properties

#endif
