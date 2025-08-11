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

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <cassert>

namespace Opm
{

template<class FluidSystem, class Indices>
WellInterfaceIndices<FluidSystem,Indices>::
WellInterfaceIndices(const Well& well,
                     const ParallelWellInfo<Scalar>& parallel_well_info,
                     const int time_step,
                     const ModelParameters& param,
                     const typename WellInterfaceFluidSystem<FluidSystem, Indices>::RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components,
                     const int num_phases,
                     const int index_of_well,
                     const std::vector<PerforationData<Scalar>>& perf_data)
    : WellInterfaceFluidSystem<FluidSystem, Indices>(well,
                                            parallel_well_info,
                                            time_step,
                                            param,
                                            rate_converter,
                                            pvtRegionIdx,
                                            num_components,
                                            num_phases,
                                            index_of_well,
                                            perf_data)
{
}

template<class FluidSystem, class Indices>
int
WellInterfaceIndices<FluidSystem,Indices>::
flowPhaseToModelCompIdx(const int phaseIdx) const
{
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx) == phaseIdx) {
            return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    }
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx) == phaseIdx) {
        return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    }
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx) == phaseIdx) {
        return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    }

    // for other phases return the index
    return phaseIdx;
}

template<class FluidSystem, class Indices>
int
WellInterfaceIndices<FluidSystem,Indices>::
modelCompIdxToFlowCompIdx(const int compIdx) const
{
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == compIdx)
        return FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == compIdx)
        return FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == compIdx)
        return FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);

    // for other phases return the index
    return compIdx;
}


template<typename FluidSystem, class Indices>
int
WellInterfaceIndices<FluidSystem, Indices>::
flowPhaseToModelPhaseIdx(const int phaseIdx) const
{
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx) == phaseIdx)
        return FluidSystem::waterPhaseIdx;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx) == phaseIdx)
        return FluidSystem::oilPhaseIdx;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx) == phaseIdx)
        return FluidSystem::gasPhaseIdx;

    // for other phases return the index
    return phaseIdx;
}


template<class FluidSystem, class Indices>
typename WellInterfaceIndices<FluidSystem,Indices>::Scalar
WellInterfaceIndices<FluidSystem,Indices>::
scalingFactor(const int phaseIdx) const
{
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx) == phaseIdx)
        return 1.0;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx) == phaseIdx)
        return 1.0;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx) == phaseIdx)
        return 0.01;
    if (Indices::enableSolvent && phaseIdx == Indices::contiSolventEqIdx )
        return 0.01;

    // we should not come this far
    assert(false);
    return 1.0;
}

#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(WellInterfaceIndices, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(WellInterfaceIndices, float)
#endif

} // namespace Opm
