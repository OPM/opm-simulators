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

template<class FluidSystem>
void FlowsContainer<FluidSystem>::
allocate(const std::size_t bufferSize,
         const unsigned numOutputNnc,
         const bool allocRestart,
         std::map<std::string, int>& rstKeywords)
{
    using Dir = FaceDir::DirEnum;

    // Flows may need to be allocated even when there is no restart due to BFLOW* summary keywords
    if (blockFlows_ ) {
        const std::array<int, 3> phaseIdxs { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs { gasCompIdx, oilCompIdx, waterCompIdx };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
            }
        }
    }

    if (!allocRestart) {
        return ;
    }

    enableFlows_ = false;
    enableFlowsn_ = false;
    const bool rstFlows = (rstKeywords["FLOWS"] > 0);
    if (rstFlows) {
        rstKeywords["FLOWS"] = 0;
        enableFlows_ = true;

        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto rstName = std::array { "FLOGASN+", "FLOOILN+", "FLOWATN+" };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                if (!blockFlows_) { // Already allocated if summary vectors requested
                    flows_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                }

                if (rstKeywords["FLOWS-"] > 0) {
                    flows_[FaceDir::ToIntersectionIndex(Dir::XMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::YMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flows_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                }

                if (numOutputNnc > 0) {
                    enableFlowsn_ = true;

                    flowsn_[compIdxs[ii]].name = rstName[ii];
                    flowsn_[compIdxs[ii]].indices.resize(numOutputNnc, -1);
                    flowsn_[compIdxs[ii]].values.resize(numOutputNnc, 0.0);
                }
            }
        }
        if (rstKeywords["FLOWS-"] > 0) {
            rstKeywords["FLOWS-"] = 0;
        }
    }

    enableFlores_ = false;
    enableFloresn_ = false;
    if (rstKeywords["FLORES"] > 0) {
        rstKeywords["FLORES"] = 0;
        enableFlores_ = true;

        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto rstName = std::array{ "FLRGASN+", "FLROILN+", "FLRWATN+" };

        for (unsigned ii = 0; ii < phaseIdxs.size(); ++ii) {
            if (FluidSystem::phaseIsActive(phaseIdxs[ii])) {
                flores_[FaceDir::ToIntersectionIndex(Dir::XPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flores_[FaceDir::ToIntersectionIndex(Dir::YPlus)][compIdxs[ii]].resize(bufferSize, 0.0);
                flores_[FaceDir::ToIntersectionIndex(Dir::ZPlus)][compIdxs[ii]].resize(bufferSize, 0.0);

                if (rstKeywords["FLORES-"] > 0) {
                    flores_[FaceDir::ToIntersectionIndex(Dir::XMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flores_[FaceDir::ToIntersectionIndex(Dir::YMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                    flores_[FaceDir::ToIntersectionIndex(Dir::ZMinus)][compIdxs[ii]].resize(bufferSize, 0.0);
                }

                if (numOutputNnc > 0) {
                    enableFloresn_ = true;

                    floresn_[compIdxs[ii]].name = rstName[ii];
                    floresn_[compIdxs[ii]].indices.resize(numOutputNnc, -1);
                    floresn_[compIdxs[ii]].values.resize(numOutputNnc, 0.0);
                }
            }
        }
        if (rstKeywords["FLORES-"] > 0) {
            rstKeywords["FLORES-"] = 0;
        }
    }
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
