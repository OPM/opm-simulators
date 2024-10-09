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

#include <config.h>

#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#include <opm/models/utils/parametersystem.hpp>

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

template void registerEclTimeSteppingParameters<double>();

#if FLOW_INSTANTIATE_FLOAT
template void registerEclTimeSteppingParameters<float>();
#endif

} // namespace Opm
