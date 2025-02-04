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
#include <opm/simulators/flow/FlowsContainer.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <algorithm>

namespace Opm {

template<class FluidSystem>
FlowsContainer<FluidSystem>::
FlowsContainer(const Schedule& schedule,
               const SummaryConfig& summaryConfig)
{
    // Check for any BFLOW[I|J|K] summary keys
    blockFlows_ = summaryConfig.keywords("BFLOW*").size() > 0;

    // Check if FLORES/FLOWS is set in any RPTRST in the schedule
    enableFlores_ = false;  // Used for the output of i+, j+, k+
    enableFloresn_ = false; // Used for the special case of nnc
    enableFlows_ = false;
    enableFlowsn_ = false;

    anyFlores_ = std::any_of(schedule.begin(), schedule.end(),
                             [](const auto& block)
                             {
                                const auto& rstkw = block.rst_config().keywords;
                                return rstkw.find("FLORES") != rstkw.end();
                             });
    anyFlows_ = std::any_of(schedule.begin(), schedule.end(),
                            [](const auto& block)
                            {
                                const auto& rstkw = block.rst_config().keywords;
                                return rstkw.find("FLOWS") != rstkw.end();
                            });
}

template<class T> using FS = BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>;

#define INSTANTIATE_TYPE(T) \
    template class FlowsContainer<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP_THREEPHASE(NUM) \
    template<class T> using FS##NUM = GenericOilGasWaterFluidSystem<T, NUM, true>; \
    template class FlowsContainer<FS##NUM<double>>;

#define INSTANTIATE_COMP_TWOPHASE(NUM) \
    template<class T> using GFS##NUM = GenericOilGasWaterFluidSystem<T, NUM, false>; \
    template class FlowsContainer<GFS##NUM<double>>;

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
