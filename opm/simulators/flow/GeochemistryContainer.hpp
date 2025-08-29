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
#ifndef GEOCHEMISTRY_CONTAINER_HPP
#define GEOCHEMISTRY_CONTAINER_HPP

#include <functional>
#include <vector>

namespace Opm {

namespace data { class Solution; }
class SpeciesConfig;
class MineralConfig;

template<class Scalar>
class GeochemistryContainer
{
    using ScalarBuffer = std::vector<Scalar>;

public:
    void allocate(const unsigned bufferSize,
                  const SpeciesConfig& species,
                  const MineralConfig& minerals);

    using AssignFunction = std::function<Scalar(const unsigned)>;

    void assignSpeciesConcentrations(const unsigned globalDofIdx,
                                     const AssignFunction& concentration);
    void assignMineralConcentrations(const unsigned globalDofIdx,
                                     const AssignFunction& concentration);
    void assignPH(const unsigned globalDofIdx, const Scalar ph);

    void outputRestart(data::Solution& sol,
                       const SpeciesConfig& species,
                       const MineralConfig& minerals);

    bool allocated() const
    { return allocated_; }

private:
    std::vector<ScalarBuffer> speciesConcentrations_{};
    std::vector<ScalarBuffer> mineralConcentrations_{};
    ScalarBuffer pH_;
    bool allocated_{false};
};  // class GeochemistryContainer

} // namespace Opm

#endif
