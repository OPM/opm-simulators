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
 * \brief Defines some fundamental parameters for all models.
 */
#ifndef EWOMS_BASIC_PARAMETERS_HH
#define EWOMS_BASIC_PARAMETERS_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

//! grid resolution
struct CellsX { static constexpr unsigned value = 1; };
struct CellsY { static constexpr unsigned value = 1; };
struct CellsZ { static constexpr unsigned value = 1; };

//! domain size
template<class Scalar>
struct DomainSizeX { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct DomainSizeY { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct DomainSizeZ { static constexpr Scalar value = 1.0; };

//! The default value for the simulation's end time
template<class TypeTag, class MyTypeTag>
struct EndTime { using type = Properties::UndefinedProperty; };

//! Name of the grid file
struct GridFile { static constexpr auto value = ""; };

//! Property which tells the Vanguard how often the grid should be refined
//! after creation.
struct GridGlobalRefinements { static constexpr unsigned value = 0; };

//! The default value for the simulation's initial time step size
template<class TypeTag, class MyTypeTag>
struct InitialTimeStepSize { using type = Properties::UndefinedProperty; };

//! Set a value for the ParameterFile property
struct ParameterFile { static constexpr auto value = ""; };

//! By default, do not force any time steps
struct PredeterminedTimeStepsFile { static constexpr auto value = ""; };

/*!
 * \brief Print all parameters on startup?
 *
 * 0 means 'no', 1 means 'yes', 2 means 'print only to logfiles'. The
 * default is 2.
 */
template<class TypeTag, class MyTypeTag>
struct PrintParameters { using type = Properties::UndefinedProperty; };

//! The default value for the simulation's restart time
template<class TypeTag, class MyTypeTag>
struct RestartTime { using type = Properties::UndefinedProperty; };

} // namespace Opm:Parameters

#endif
