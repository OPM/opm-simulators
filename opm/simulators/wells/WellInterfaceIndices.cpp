/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <cassert>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
WellInterfaceIndices<FluidSystem,Indices,Scalar>::
WellInterfaceIndices(const Well& well,
                     const ParallelWellInfo& parallel_well_info,
                     const int time_step,
                     const typename WellInterfaceFluidSystem<FluidSystem>::RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components,
                     const int num_phases,
                     const int index_of_well,
                     const std::vector<PerforationData>& perf_data)
    : WellInterfaceFluidSystem<FluidSystem>(well,
                                            parallel_well_info,
                                            time_step,
                                            rate_converter,
                                            pvtRegionIdx,
                                            num_components,
                                            num_phases,
                                            index_of_well,
                                            perf_data)
{
}

template<class FluidSystem, class Indices, class Scalar>
int
WellInterfaceIndices<FluidSystem,Indices,Scalar>::
flowPhaseToEbosCompIdx(const int phaseIdx) const
{
    const auto& pu = this->phaseUsage();
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
        return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
        return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
        return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

    // for other phases return the index
    return phaseIdx;
}

template<class FluidSystem, class Indices, class Scalar>
int
WellInterfaceIndices<FluidSystem,Indices,Scalar>::
ebosCompIdxToFlowCompIdx(const unsigned compIdx) const
{
    const auto& pu = this->phaseUsage();
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == compIdx)
        return pu.phase_pos[Water];
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == compIdx)
        return pu.phase_pos[Oil];
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == compIdx)
        return pu.phase_pos[Gas];

    // for other phases return the index
    return compIdx;
}

template<class FluidSystem, class Indices, class Scalar>
double
WellInterfaceIndices<FluidSystem,Indices,Scalar>::
scalingFactor(const int phaseIdx) const
{
    const auto& pu = this->phaseUsage();
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
        return 1.0;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
        return 1.0;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
        return 0.01;
    if (Indices::enableSolvent && phaseIdx == Indices::contiSolventEqIdx )
        return 0.01;

    // we should not come this far
    assert(false);
    return 1.0;
}

#define INSTANCE( ...) \
template class WellInterfaceIndices<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                                    __VA_ARGS__, \
                                    double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

} // namespace Opm
