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
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <tuple>

namespace {

    template<class Scalar>
    using DataEntry = std::tuple<std::string,
                                 Opm::UnitSystem::measure,
                                 std::array<std::vector<Scalar>,6>&>;

    template<int idx, class Array, class Scalar>
    void addEntry(std::vector<DataEntry<Scalar>>& container,
                  const std::string& name,
                  Opm::UnitSystem::measure measure,
                  Array& flowArray)
    {
        if constexpr (idx >= 0) {  // Only add if index is valid
            container.emplace_back(name, measure, flowArray[idx]);
        }
    }

    template<int idx, class Array, class Scalar>
    void assignToVec(Array& array,
                     const unsigned faceId,
                     const unsigned globalDofIdx,
                     const Scalar value)
    {
        if constexpr (idx != -1) {
            if (!array[idx][faceId].empty()) {
                array[idx][faceId][globalDofIdx] = value;
            }
        }
    }

    template<int idx, class Array, class Scalar>
    void assignToNnc(Array& array,
                     unsigned nncId,
                     const Scalar value)
    {
        if constexpr (idx != -1) {
            if (!array[idx].indices.empty()) {
                array[idx].indices[nncId] = nncId;
                array[idx].values[nncId] = value;
            }
        }
    }

}

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
                flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::XPlus)].resize(bufferSize, 0.0);
                flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::YPlus)].resize(bufferSize, 0.0);
                flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::ZPlus)].resize(bufferSize, 0.0);
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
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::XPlus)].resize(bufferSize, 0.0);
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::YPlus)].resize(bufferSize, 0.0);
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::ZPlus)].resize(bufferSize, 0.0);
                }

                if (rstKeywords["FLOWS-"] > 0) {
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::XMinus)].resize(bufferSize, 0.0);
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::YMinus)].resize(bufferSize, 0.0);
                    flows_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::ZMinus)].resize(bufferSize, 0.0);
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
                flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::XPlus)].resize(bufferSize, 0.0);
                flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::YPlus)].resize(bufferSize, 0.0);
                flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::ZPlus)].resize(bufferSize, 0.0);

                if (rstKeywords["FLORES-"] > 0) {
                    flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::XMinus)].resize(bufferSize, 0.0);
                    flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::YMinus)].resize(bufferSize, 0.0);
                    flores_[compIdxs[ii]][FaceDir::ToIntersectionIndex(Dir::ZMinus)].resize(bufferSize, 0.0);
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

template<class FluidSystem>
void FlowsContainer<FluidSystem>::
assignFlores(const unsigned globalDofIdx,
             const int faceId,
             const unsigned nncId,
             const Scalar gas,
             const Scalar oil,
             const Scalar water)
{
    if (faceId >= 0) {
        assignToVec<gasCompIdx>(this->flores_, faceId, globalDofIdx, gas);
        assignToVec<oilCompIdx>(this->flores_, faceId, globalDofIdx, oil);
        assignToVec<waterCompIdx>(this->flores_, faceId, globalDofIdx, water);
    }
    else if (faceId == -2) {
        assignToNnc<gasCompIdx>(this->floresn_, nncId, gas);
        assignToNnc<oilCompIdx>(this->floresn_, nncId, oil);
        assignToNnc<waterCompIdx>(this->floresn_, nncId, water);
    }
}

template<class FluidSystem>
void FlowsContainer<FluidSystem>::
assignFlows(const unsigned globalDofIdx,
            const int faceId,
            const unsigned nncId,
            const Scalar gas,
            const Scalar oil,
            const Scalar water)
{
    if (faceId >= 0) {
        assignToVec<gasCompIdx>(this->flows_, faceId, globalDofIdx, gas);
        assignToVec<oilCompIdx>(this->flows_, faceId, globalDofIdx, oil);
        assignToVec<waterCompIdx>(this->flows_, faceId, globalDofIdx, water);
    }
    else if (faceId == -2) {
        assignToNnc<gasCompIdx>(this->flowsn_, nncId, gas);
        assignToNnc<oilCompIdx>(this->flowsn_, nncId, oil);
        assignToNnc<waterCompIdx>(this->flowsn_, nncId, water);
    }
}

template<class FluidSystem>
void FlowsContainer<FluidSystem>::
outputRestart(data::Solution& sol)
{
    auto doInsert = [&sol](ScalarBuffer& value,
                           const std::string& name,
                           UnitSystem::measure measure)
    {
        if (!value.empty()) {
            sol.insert(name, measure, std::move(value),
                       data::TargetType::RESTART_SOLUTION);
        }
    };

    using Dir = FaceDir::DirEnum;
    std::vector<DataEntry<Scalar>> entries;
    if (this->enableFlores_) {
        addEntry<gasCompIdx>  (entries, "FLRGAS", UnitSystem::measure::rate,                flores_);
        addEntry<oilCompIdx>  (entries, "FLROIL", UnitSystem::measure::rate,                flores_);
        addEntry<waterCompIdx>(entries, "FLRWAT", UnitSystem::measure::rate,                flores_);
    }
    if (this->enableFlows_) {
        addEntry<gasCompIdx>  (entries, "FLOGAS", UnitSystem::measure::gas_surface_rate,    flows_);
        addEntry<oilCompIdx>  (entries, "FLOOIL", UnitSystem::measure::liquid_surface_rate, flows_);
        addEntry<waterCompIdx>(entries, "FLOWAT", UnitSystem::measure::liquid_surface_rate, flows_);
    }

    std::for_each(entries.begin(), entries.end(),
                  [&doInsert](auto& array)
                  {
                      static const auto dirs = std::array{
                          std::pair{FaceDir::ToIntersectionIndex(Dir::XMinus), "I-"},
                          std::pair{FaceDir::ToIntersectionIndex(Dir::XPlus), "I+"},
                          std::pair{FaceDir::ToIntersectionIndex(Dir::YMinus), "J-"},
                          std::pair{FaceDir::ToIntersectionIndex(Dir::YPlus), "J+"},
                          std::pair{FaceDir::ToIntersectionIndex(Dir::ZMinus), "K-"},
                          std::pair{FaceDir::ToIntersectionIndex(Dir::ZPlus), "K+"},
                      };
                      const auto& name = std::get<0>(array);
                      const auto& measure = std::get<1>(array);
                      auto& value = std::get<2>(array);
                      for (const auto& [index, postfix] : dirs) {
                          doInsert(value[index], name + postfix, measure);
                      }
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
