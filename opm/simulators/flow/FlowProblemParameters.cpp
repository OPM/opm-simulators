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
#include <opm/simulators/flow/FlowProblemParameters.hpp>

#include <opm/models/common/multiphasebaseparameters.hh>
#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/EclWriter.hpp>

#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisParameters.hpp>
#endif

namespace Opm {

template<class Scalar>
void registerFlowProblemParameters()
{
    Parameters::Register<Parameters::EnableWriteAllSolutions>
       ("Write all solutions to disk instead of only the ones for the "
        "report steps");
#if HAVE_DAMARIS
    Parameters::Register<Parameters::EnableDamarisOutput>
        ("Write a specific variable using Damaris in a separate core");
#endif
    Parameters::Register<Parameters::EclOutputDoublePrecision>
        ("Tell the output writer to use double precision. Useful for 'perfect' restarts");
    Parameters::Register<Parameters::RestartWritingInterval>
        ("The frequencies of which time steps are serialized to disk");
    Parameters::Register<Parameters::EnableDriftCompensation>
        ("Enable partial compensation of systematic mass losses via "
         "the source term of the next time step");
    Parameters::Register<Parameters::OutputMode>
        ("Specify which messages are going to be printed. "
         "Valid values are: none, log, all (default)");
    Parameters::Register<Parameters::NumPressurePointsEquil>
        ("Number of pressure points (in each direction) in tables used for equilibration");
    Parameters::Hide<Parameters::NumPressurePointsEquil>(); // Users will typically not need to modify this parameter..
    Parameters::Register<Parameters::ExplicitRockCompaction>
        ("Use pressure from end of the last time step when evaluating rock compaction");
    Parameters::Hide<Parameters::ExplicitRockCompaction>(); // Users will typically not need to modify this parameter..

    Parameters::Register<Parameters::CheckSatfuncConsistency>
        ("Whether or not to check saturation function consistency requirements");

    Parameters::Register<Parameters::NumSatfuncConsistencySamplePoints>
        ("Maximum number of reported failures for each individual saturation function consistency check");
    
    Parameters::Register<Parameters::HyNeConfigFile>
        ("Use config files for Hybrid Newton");

    Parameters::Register<Parameters::UseHyNe>
        ("Wheter or not to use Hybrid Newton");

    // By default, stop it after the universe will probably have stopped
    // to exist. (the ECL problem will finish the simulation explicitly
    // after it simulated the last episode specified in the deck.)
    Parameters::SetDefault<Parameters::EndTime<Scalar>>(1e100);
    // The chosen value means that the size of the first time step is the
    // one of the initial episode (if the length of the initial episode is
    // not millions of trillions of years, that is...)
    Parameters::SetDefault<Parameters::InitialTimeStepSize<Scalar>>(3600*24);
    // Disable the VTK output by default for this problem ...
    Parameters::SetDefault<Parameters::EnableVtkOutput>(false);
    // the cache for intensive quantities can be used for ECL problems and also yields a
    // decent speedup...
    Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);
    // the cache for the storage term can also be used and also yields a decent speedup
    Parameters::SetDefault<Parameters::EnableStorageCache>(true);
    // the default for the allowed volumetric error for oil per second
    Parameters::SetDefault<Parameters::NewtonTolerance<Scalar>>(1e-2);
    Parameters::SetDefault<Parameters::EnableGravity>(true);
}

template void registerFlowProblemParameters<double>();

#if FLOW_INSTANTIATE_FLOAT
template void registerFlowProblemParameters<float>();
#endif

}
