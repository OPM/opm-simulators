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
/*!
 * \file
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_RFT_CONTAINER_HPP
#define OPM_RFT_CONTAINER_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Opm {

namespace data { class Wells; }
class EclipseState;
class Schedule;

template<class FluidSystem>
class RFTContainer {
    using Scalar = typename FluidSystem::Scalar;
    using ScalarBuffer = std::vector<Scalar>;

    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

public:
    using AssignmentFunc = std::function<Scalar()>;
    using WellQueryFunc = std::function<bool(const std::string&)>;

    RFTContainer(const EclipseState& eclState,
                 const Schedule& schedule,
                 const WellQueryFunc& wellQuery)
        : eclState_(eclState)
        , schedule_(schedule)
        , wellQuery_(wellQuery)
    {}

   void addToWells(data::Wells& wellDatas,
                   const std::size_t reportStepNum,
                   const Parallel::Communication& comm);

    void allocate(const std::size_t reportStepNum);

   void assign(const unsigned cartesianIndex,
                const AssignmentFunc& oil,
                const AssignmentFunc& water,
                const AssignmentFunc& gas);

private:
    const EclipseState& eclState_;
    const Schedule& schedule_;
    WellQueryFunc wellQuery_;

    std::map<std::size_t, Scalar> oilConnectionPressures_;
    std::map<std::size_t, Scalar> waterConnectionSaturations_;
    std::map<std::size_t, Scalar> gasConnectionSaturations_;
};

} // namespace Opm

#endif // OPM_RFT_CONTAINER_HPP
