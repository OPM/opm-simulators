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
 * \copydoc Opm::FlowProblem
 */
#ifndef OPM_FLOW_PROBLEM_PARAMETERS_HPP
#define OPM_FLOW_PROBLEM_PARAMETERS_HPP

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

// Enable the additional checks even if compiled in debug mode (i.e., with the NDEBUG
// macro undefined). Next to a slightly better performance, this also eliminates some
// print statements in debug mode.
template<class TypeTag, class MyTypeTag>
struct EnableDebuggingChecks { using type = Properties::UndefinedProperty; };

// Enable partial compensation of systematic mass losses via the source term of the next time
// step
template<class TypeTag, class MyTypeTag>
struct EnableDriftCompensation { using type = Properties::UndefinedProperty; };

// implicit or explicit pressure in rock compaction
template<class TypeTag, class MyTypeTag>
struct ExplicitRockCompaction { using type = Properties::UndefinedProperty; };

// Parameterize equilibration accuracy
template<class TypeTag, class MyTypeTag>
struct NumPressurePointsEquil { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct OutputMode { using type = Properties::UndefinedProperty; };

// The number of time steps skipped between writing two consequtive restart files
template<class TypeTag, class MyTypeTag>
struct RestartWritingInterval { using type = Properties::UndefinedProperty; };

} // namespace Opm::Parameters

#endif // OPM_FLOW_PROBLEM_PARAMETERS_HPP
