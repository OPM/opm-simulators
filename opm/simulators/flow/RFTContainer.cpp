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
#include <opm/simulators/flow/RFTContainer.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/output/data/Wells.hpp>

#include <algorithm>
#include <tuple>

namespace {

template<class Scalar>
void gatherAndUpdateMap(std::map<std::size_t, Scalar>& local_map,
                        const Opm::Parallel::Communication& comm)
{
    std::vector<std::pair<int, Scalar>> pairs(local_map.begin(), local_map.end());
    std::vector<std::pair<int, Scalar>> all_pairs;
    std::vector<int> offsets;

    std::tie(all_pairs, offsets) = Opm::allGatherv(pairs, comm);

    // Update maps on all ranks
    for (auto i = static_cast<std::size_t>(offsets[0]); i < all_pairs.size(); ++i) {
        const auto& key_value = all_pairs[i];
        if (auto candidate = local_map.find(key_value.first); candidate != local_map.end()) {
            const Scalar prev_value = candidate->second;
            candidate->second = std::max(prev_value, key_value.second);
        }
        else {
            local_map[key_value.first] = key_value.second;
        }
    }
}

}

namespace Opm {

template<class FluidSystem>
void RFTContainer<FluidSystem>::
addToWells(data::Wells& wellDatas,
           const std::size_t reportStepNum,
           const Parallel::Communication& comm)
{
    const auto& rft_config = schedule_[reportStepNum].rft_config();
    if (!rft_config.active()) {
        return;
    }

    if (comm.size() > 1) {
        gatherAndUpdateMap(oilConnectionPressures_, comm);
        gatherAndUpdateMap(waterConnectionSaturations_, comm);
        gatherAndUpdateMap(gasConnectionSaturations_, comm);
    }

    for (const auto& wname : this->schedule_[reportStepNum].well_order()) {
        // Don't bother with wells not on this process.
        if (!wellIsOwnedByCurrentRank_(wname)) {
            continue;
        }

        const auto& [wellDataPos, inserted] = wellDatas.try_emplace(wname);

        if (inserted) {
            // New well.  Typically a shut/inactive one.  Allocate
            // data::Connection result objects for this well.
            const auto& conns = this->schedule_[reportStepNum]
                .wells(wname).getConnections();

            auto& xcon = wellDataPos->second.connections;

            xcon.reserve(conns.size());
            for (const auto& conn : conns) {
                xcon.emplace_back().index = conn.global_index();
            }
        }

        auto cond_assign = [](double& dest, unsigned idx, const auto& map)
        {
            auto it = map.find(idx);
            if (it != map.end()) {
                dest = it->second;
            }
        };

        std::for_each(wellDataPos->second.connections.begin(),
                      wellDataPos->second.connections.end(),
                      [&cond_assign, this](auto& connectionData)
                      {
                          const auto index = connectionData.index;
                          cond_assign(connectionData.cell_pressure, index,
                                      oilConnectionPressures_);
                          cond_assign(connectionData.cell_saturation_water, index,
                                      waterConnectionSaturations_);
                          cond_assign(connectionData.cell_saturation_gas, index,
                                      gasConnectionSaturations_);
                       });
    }

    oilConnectionPressures_.clear();
    waterConnectionSaturations_.clear();
    gasConnectionSaturations_.clear();
}

template<class FluidSystem>
void RFTContainer<FluidSystem>::
allocate(const std::size_t reportStepNum)
{
    // Well RFT data
    const auto& rft_config = schedule_[reportStepNum].rft_config();
    if (!rft_config.active()) {
        return;
    }

    for (const auto& wname : schedule_[reportStepNum].well_order()) {
        // don't bother with wells not on this process
        if (!wellOnCurrentRank_(wname)) {
            continue;
        }

        const auto& well = schedule_[reportStepNum].wells.get(wname);

        for (const auto& connection: well.getConnections()) {
            const std::size_t index = connection.global_index();

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                oilConnectionPressures_.emplace(index, 0.0);
            }

            if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                waterConnectionSaturations_.emplace(index, 0.0);
            }

            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                gasConnectionSaturations_.emplace(index, 0.0);
            }
        }
    }
}

template<class FluidSystem>
void RFTContainer<FluidSystem>::
assign(const unsigned cartesianIndex,
       const AssignmentFunc& oil,
       const AssignmentFunc& water,
       const AssignmentFunc& gas)
{
    auto cond_assign = [](auto& map, unsigned idx, const auto& func)
    {
        auto it = map.find(idx);
        if (it != map.end()) {
            it->second = func();
        }
    };

    cond_assign(oilConnectionPressures_, cartesianIndex, oil);
    cond_assign(waterConnectionSaturations_, cartesianIndex, water);
    cond_assign(gasConnectionSaturations_, cartesianIndex, gas);
}

template<class T> using FS = BlackOilFluidSystem<T, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE_TYPE(T) \
    template class RFTContainer<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP_THREEPHASE(NUM) \
    template<class T> using FS##NUM = GenericOilGasWaterFluidSystem<T, NUM, true>; \
    template class RFTContainer<FS##NUM<double>>;

#define INSTANTIATE_COMP_TWOPHASE(NUM) \
    template<class T> using GFS##NUM = GenericOilGasWaterFluidSystem<T, NUM, false>; \
    template class RFTContainer<GFS##NUM<double>>;

#define INSTANTIATE_COMP(NUM) \
    INSTANTIATE_COMP_THREEPHASE(NUM) \
    INSTANTIATE_COMP_TWOPHASE(NUM)

INSTANTIATE_COMP_THREEPHASE(0)  // \Note: to register the parameter ForceDisableFluidInPlaceOutput
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
