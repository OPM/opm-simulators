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
#ifndef OPM_FLOWS_CONTAINER_HPP
#define OPM_FLOWS_CONTAINER_HPP

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/simulators/flow/FlowsData.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Opm {

namespace data { class Solution; }
class Schedule;
class SummaryConfig;

template<class FluidSystem>
class FlowsContainer
{
    using Scalar = typename FluidSystem::Scalar;
    using ScalarBuffer = std::vector<Scalar>;

    static constexpr auto numPhases = FluidSystem::numPhases;
    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr auto gasCompIdx = FluidSystem::gasCompIdx;
    static constexpr auto oilCompIdx = FluidSystem::oilCompIdx;
    static constexpr auto waterCompIdx = FluidSystem::waterCompIdx;

public:
    FlowsContainer(const Schedule& schedule,
                   const SummaryConfig& summaryConfig,
                   std::function<bool(const int)> );

    void allocate(const std::size_t bufferSize,
                  const SummaryConfig& summaryConfig,
                  const unsigned numOutputNnc,
                  const bool allocRestart,
                  std::map<std::string, int>& rstKeywords);

    void assignFlores(const unsigned globalDofIdx,
                      const int faceId,
                      const unsigned nncId,
                      const Scalar gas,
                      const Scalar oil,
                      const Scalar water);

    void assignBlockVelocity(const unsigned globalDofIdx,
                             const int faceId,
                             const int comp_idx,
                             const Scalar velocity);

    void assignBlockFlows(const unsigned globalDofIdx,
                          const int faceId,
                          const int comp_idx,
                          const Scalar flow);

    void assignFlows(const unsigned globalDofIdx,
                     const int faceId,
                     const unsigned nncId,
                     const Scalar gas,
                     const Scalar oil,
                     const Scalar water);

    void outputRestart(data::Solution& sol);

    const std::array<FlowsData<double>, 3>& getFlowsn() const
    { return this->flowsn_; }

    bool hasFlowsn() const
    { return enableFlowsn_; }

    bool hasFlows() const
    { return enableFlows_; }

    const std::vector<int> blockVelocity() const
    { return blockVelocityAllIds_; }

    const std::vector<int> blockFlows() const
    { return blockFlowsAllIds_; }

    unsigned blockFlowsIds(const unsigned globalDofIdx,
                           const int dir,
                           const int comp_idx) const
    { 
        const auto& blockIdxs = blockFlowsIds_[comp_idx][dir];
        auto it = std::lower_bound(blockIdxs.begin(), blockIdxs.end(), globalDofIdx);
        return std::distance(blockIdxs.begin(), it); 
    }

    unsigned blockVelocityIds(const unsigned globalDofIdx,
                              const int dir,
                              const int comp_idx) const
    { 
        const auto& blockIdxs = blockVelocityIds_[comp_idx][dir];
        auto it = std::lower_bound(blockIdxs.begin(), blockIdxs.end(), globalDofIdx);
        return std::distance(blockIdxs.begin(), it); 
    }

    bool anyFlows() const
    { return anyFlows_; }

    const std::array<FlowsData<double>, 3>& getFloresn() const
    { return this->floresn_; }

    bool hasFloresn() const
    { return enableFloresn_; }

    bool hasFlores() const
    { return enableFlores_; }

    bool anyFlores() const
    { return anyFlores_; }

    bool hasBlockVelocityValue(const unsigned globalDofIdx,
                               const int dir,
                               const int comp_idx) const
    {   
        const auto& blockIdxs = blockVelocityIds_[comp_idx][dir];
        return std::ranges::binary_search(blockIdxs, globalDofIdx);
    }

    bool hasBlockFlowValue(const unsigned globalDofIdx,
                           const int dir,
                           const int comp_idx) const
    {   
        const auto& blockIdxs = blockFlowsIds_[comp_idx][dir];
        return std::ranges::binary_search(blockIdxs, globalDofIdx);
    }

    Scalar getVelocity(const unsigned globalDofIdx,
                       const FaceDir::DirEnum dir,
                       const int comp_idx) const
    { return velocity_[comp_idx][FaceDir::ToIntersectionIndex(dir)][globalDofIdx]; }

    Scalar getFlow(const unsigned globalDofIdx,
                   const FaceDir::DirEnum dir,
                   const int comp_idx) const
    { return flows_[comp_idx][FaceDir::ToIntersectionIndex(dir)][globalDofIdx]; }

private:
    bool anyFlows_{false};
    bool anyFlores_{false};
    bool enableFlows_{false};
    bool enableFlores_{false};
    bool enableFlowsn_{false};
    bool enableFloresn_{false};

    std::array<std::array<ScalarBuffer, 6>, numPhases> flows_;
    std::array<std::array<ScalarBuffer, 6>, numPhases> flores_;
    std::array<std::array<ScalarBuffer, 6>, numPhases> velocity_;

    std::array<FlowsData<double>, 3> floresn_;
    std::array<FlowsData<double>, 3> flowsn_;

    std::vector<int> blockFlowsAllIds_;
    std::vector<int> blockVelocityAllIds_;
    std::array<std::array<std::vector<int>, 6>, numPhases> blockFlowsIds_;
    std::array<std::array<std::vector<int>, 6>, numPhases> blockVelocityIds_;
};

} // namespace Opm

#endif // OPM_FLOWS_CONTAINER_HPP
