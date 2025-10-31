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

#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>

namespace Opm::Parameters {

// Enable partial compensation of systematic mass losses via
// the source term of the next time step
struct EnableDriftCompensation { static constexpr bool value = false; };

// implicit or explicit pressure in rock compaction
struct ExplicitRockCompaction { static constexpr bool value = false; };

// Whether or not to check saturation function consistency requirements.
struct CheckSatfuncConsistency { static constexpr bool value = true; };

// Maximum number of reported failures for each saturation function
// consistency check.
struct NumSatfuncConsistencySamplePoints { static constexpr int value = 5; };

// Parameterize equilibration accuracy
struct NumPressurePointsEquil
{ static constexpr int value = ParserKeywords::EQLDIMS::DEPTH_NODES_P::defaultValue; };

struct OutputMode { static constexpr auto value = "all"; };

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
struct RestartWritingInterval { static constexpr int value = 0xffffff; }; // disable

// Path to the config file containing all Hybrid Newton parameters
struct HyNeConfigFile { static constexpr auto value = "hybridNewtonConfig.json"; };
// Wheter or not to use Hybrid Newton nonlinear preconditioning
struct UseHyNe { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

template<class Scalar>
void registerFlowProblemParameters();

}

#endif // OPM_FLOW_PROBLEM_PARAMETERS_HPP
