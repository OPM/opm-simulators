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
 * \ingroup BlackOilModel
 *
 * \brief Declares the properties required by the black oil model.
 */
#ifndef EWOMS_BLACK_OIL_PROPERTIES_HH
#define EWOMS_BLACK_OIL_PROPERTIES_HH

#include <opm/models/common/multiphasebaseproperties.hh>

namespace Opm::Properties {

//! Specifies if the simulation should write output files that are
//! compatible with those produced by the commercial Eclipse simulator
template<class TypeTag, class MyTypeTag>
struct EnableEclipseOutput { using type = UndefinedProperty; };
//! Enable the ECL-blackoil extension for solvents. ("Second gas")
template<class TypeTag, class MyTypeTag>
struct EnableSolvent { using type = UndefinedProperty; };
//! Enable the ECL-blackoil extension for extended BO. ("Second gas" - alternative approach)
template<class TypeTag, class MyTypeTag>
struct EnableExtbo { using type = UndefinedProperty; };
//! Enable the ECL-blackoil extension for polymer.
template<class TypeTag, class MyTypeTag>
struct EnablePolymer { using type = UndefinedProperty; };
//! Enable the tracking polymer molecular weight tracking and related functionalities
template<class TypeTag, class MyTypeTag>
struct EnablePolymerMW { using type = UndefinedProperty; };
//! Enable surface volume scaling
template<class TypeTag, class MyTypeTag>
struct BlackoilConserveSurfaceVolume { using type = UndefinedProperty; };
//! Enable the ECL-blackoil extension for foam
template<class TypeTag, class MyTypeTag>
struct EnableFoam { using type = UndefinedProperty; };
//! Enable the ECL-blackoil extension for salt
template<class TypeTag, class MyTypeTag>
struct EnableBrine { using type = UndefinedProperty; };


//! Allow the spatial and temporal domains to exhibit non-constant temperature
//! in the black-oil model
template<class TypeTag, class MyTypeTag>
struct EnableTemperature { using type = UndefinedProperty; };

//! The relative weight of the residual of the energy equation compared to the mass
//! residuals
//!
//! this is basically a hack to work around the limitation that the convergence criterion
//! of unmodified dune-istl linear solvers cannot weight the individual equations. if the
//! energy equation is not scaled, its absolute value is normally several orders of
//! magnitude larger than that of the mass balance equations
template<class TypeTag, class MyTypeTag>
struct BlackOilEnergyScalingFactor { using type = UndefinedProperty; };


} // namespace Opm::Properties

#endif
