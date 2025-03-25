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
#ifndef OPM_TRACER_CONTAINER_HPP
#define OPM_TRACER_CONTAINER_HPP

#include <functional>
#include <vector>

namespace Opm {

namespace data { class Solution; }
class TracerConfig;

template<class FluidSystem>
class TracerContainer
{
    using Scalar = typename FluidSystem::Scalar;
    using ScalarBuffer = std::vector<Scalar>;

public:
    void allocate(const unsigned bufferSize,
                  const TracerConfig& tracers);

    using AssignFunction = std::function<Scalar(const unsigned)>;

    void assignFreeConcentrations(const unsigned globalDofIdx,
                                  const AssignFunction& concentration);

    void assignSolConcentrations(const unsigned globalDofIdx,
                                 const AssignFunction& concentration);

    void outputRestart(data::Solution& sol,
                       const TracerConfig& tracers);

private:
    std::vector<ScalarBuffer> freeConcentrations_{};
    std::vector<ScalarBuffer> solConcentrations_{};
    bool allocated_{false};
};

} // namespace Opm

#endif // OPM_TRACER_CONTAINER_HPP
