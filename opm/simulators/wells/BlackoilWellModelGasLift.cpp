/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/BlackoilWellModelGasLift.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

namespace Opm {

template<typename FluidSystem, typename Indices>
void BlackoilWellModelGasLiftGeneric<FluidSystem, Indices>::
gliftDebug([[maybe_unused]] const std::string& msg,
           [[maybe_unused]] DeferredLogger& deferred_logger) const
{
    if constexpr (glift_debug) {
        if (terminal_output_) {
            const std::string message =
                fmt::format("  GLIFT (DEBUG) : BlackoilWellModel : {}", msg);
            deferred_logger.debug(message);
        }
    }
}

template<typename FluidSystem, typename Indices>
void BlackoilWellModelGasLiftGeneric<FluidSystem, Indices>::
gliftDebugShowALQ(const std::vector<WellInterfaceGeneric<FluidSystem, Indices>*>& well_container,
                  const WellState<FluidSystem, Indices>& wellState,
                  DeferredLogger& deferred_logger)
{
    for (const auto& well : well_container) {
        if (well->isProducer()) {
            const auto alq = wellState.well(well->name()).alq_state.get();
            const std::string msg = fmt::format("ALQ_REPORT : {} : {}",
                                                well->name(), alq);
            gliftDebug(msg, deferred_logger);
        }
    }
}

// If a group has any production rate constraints, and/or a limit
// on its total rate of lift gas supply,  allocate lift gas
// preferentially to the wells that gain the most benefit from
// it. Lift gas increments are allocated in turn to the well that
// currently has the largest weighted incremental gradient. The
// procedure takes account of any limits on the group production
// rate or lift gas supply applied to any level of group.
template<typename FluidSystem, typename Indices>
void BlackoilWellModelGasLiftGeneric<FluidSystem, Indices>::
gasLiftOptimizationStage2(const Parallel::Communication& comm,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          WellState<FluidSystem, Indices>& wellState,
                          GroupState<Scalar>& groupState,
                          GLiftProdWells& prod_wells,
                          GLiftOptWells& glift_wells,
                          GasLiftGroupInfo<FluidSystem, Indices>& group_info,
                          GLiftWellStateMap& glift_well_state_map,
                          const int episodeIndex,
                          DeferredLogger& deferred_logger)

{
    OPM_TIMEFUNCTION();
    GasLiftStage2 glift {episodeIndex,
                         comm,
                         schedule,
                         summaryState,
                         deferred_logger,
                         wellState,
                         groupState,
                         prod_wells,
                         glift_wells,
                         group_info,
                         glift_well_state_map,
                         this->glift_debug
    };
    glift.runOptimize();
}

    template<class Scalar>
    using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultIndexTraits>;

#define INSTANTIATE(T,...) \
    template class BlackoilWellModelGasLiftGeneric<FS<T>, __VA_ARGS__>;

#define INSTANTIATE_TYPE(T)                                                  \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

    INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
    INSTANTIATE_TYPE(float)
#endif

}
