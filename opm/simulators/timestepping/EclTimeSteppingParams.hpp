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
#ifndef OPM_ECL_TIMESTEPPING_PARAMS_HPP
#define OPM_ECL_TIMESTEPPING_PARAMS_HPP

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

struct EnableTuning { static constexpr bool value = false; };
template<class Scalar>
struct SolverGrowthFactor { static constexpr Scalar value = 2.0; };

template<class Scalar>
struct SolverMaxGrowth { static constexpr Scalar value = 3.0; };

template<class Scalar>
struct SolverMinTimeStep { static constexpr Scalar value = 1e-12; };

template<class Scalar>
struct SolverMaxTimeStepInDays { static constexpr Scalar value = 365.0; };

template<class Scalar>
struct SolverRestartFactor { static constexpr Scalar value = 0.33; };

template<class Scalar>
struct TimeStepAfterEventInDays { static constexpr Scalar value = -1.0; };

} // namespace Opm::Properties

namespace Opm {

template<class Scalar>
void registerEclTimeSteppingParameters()
{
    Parameters::Register<Parameters::EnableTuning>
        ("Honor some aspects of the TUNING keyword.");
    Parameters::Register<Parameters::SolverGrowthFactor<Scalar>>
        ("The factor time steps are elongated after a successful substep");
    Parameters::Register<Parameters::SolverMaxGrowth<Scalar>>
        ("The maximum factor time steps are elongated after a report step");
    Parameters::Register<Parameters::SolverMaxTimeStepInDays<Scalar>>
        ("The maximum size of a time step in days");
    Parameters::Register<Parameters::SolverMinTimeStep<Scalar>>
        ("The minimum size of a time step in days for field and "
         "metric and hours for lab. If a step cannot converge without "
         "getting cut below this step size the simulator will stop");
    Parameters::Register<Parameters::SolverRestartFactor<Scalar>>
        ("The factor time steps are elongated after restarts");
    Parameters::Register<Parameters::TimeStepAfterEventInDays<Scalar>>
        ("Time step size of the first time step after an event "
         "occurs during the simulation in days");
}

} // namespace Opm

#endif // OPM_ECL_TIME_STEPPING_PARAMS_HPP
