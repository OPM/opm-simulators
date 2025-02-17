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

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

namespace Opm {

template<class FluidSystem>
void RFTContainer<FluidSystem>::
allocate(const std::size_t reportStepNum,
         const WellQueryFunc& wellQuery)
{
    // Well RFT data
    const auto& rft_config = schedule_[reportStepNum].rft_config();
    for (const auto& well: schedule_.getWells(reportStepNum)) {

        // don't bother with wells not on this process
        if (wellQuery(well.name())) {
            continue;
        }

        if (!rft_config.active()) {
            continue;
        }

        for (const auto& connection: well.getConnections()) {
            const std::size_t i = std::size_t(connection.getI());
            const std::size_t j = std::size_t(connection.getJ());
            const std::size_t k = std::size_t(connection.getK());
            const std::size_t index = eclState_.gridDims().getGlobalIndex(i, j, k);

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                oilConnectionPressures_.emplace(std::make_pair(index, 0.0));
            }

            if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                waterConnectionSaturations_.emplace(std::make_pair(index, 0.0));
            }

            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                gasConnectionSaturations_.emplace(std::make_pair(index, 0.0));
            }
        }
    }
}

template<class T> using FS = BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>;

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
