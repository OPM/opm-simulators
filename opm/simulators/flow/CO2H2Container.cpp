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
#include <opm/simulators/flow/CO2H2Container.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <opm/simulators/flow/SolutionContainers.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <tuple>

namespace Opm {

template<class Scalar>
void CO2H2Container<Scalar>::
allocate(const unsigned bufferSize, const bool isCO2)
{
    if (isCO2) {
        cXmfCO2_.resize(bufferSize, 0.0);
    }
    else {
        cXmfH2_.resize(bufferSize, 0.0);
    }
    cYmfwat_.resize(bufferSize, 0.0);
    allocated_ = true;
}

template<class Scalar>
void CO2H2Container<Scalar>::
assign(const unsigned globalDofIdx,
       const Scalar xfrac,
       const Scalar yfrac,
       const bool isCO2)
{
    if (isCO2) {
        cXmfCO2_[globalDofIdx] = xfrac;
    }
    else {
        cXmfH2_[globalDofIdx] = xfrac;
    }
    cYmfwat_[globalDofIdx] = yfrac;
}

template<class Scalar>
CO2H2SolutionContainer<Scalar>
CO2H2Container<Scalar>::
getSolution() const
{
    return {
        cXmfCO2_,
        cXmfH2_,
        cYmfwat_
    };
}

template<class Scalar>
void CO2H2Container<Scalar>::
outputRestart(data::Solution& sol)
{
    if (!this->allocated_) {
        return;
    }

    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, std::vector<Scalar>&>;

    auto insert = [&sol](auto& entry)

    {
        if (!std::get<2>(entry).empty()) {
            sol.insert(std::get<std::string>(entry),
            std::get<UnitSystem::measure>(entry),
            std::move(std::get<2>(entry)),
            data::TargetType::RESTART_OPM_EXTENDED);
        }
    };

    auto solutionXCO2 = DataEntry{"XMFCO2", UnitSystem::measure::identity, cXmfCO2_};
    insert(solutionXCO2);
    auto solutionXH2 = DataEntry{"XMFH2", UnitSystem::measure::identity, cXmfH2_};
    insert(solutionXH2);
    auto solutionYfrac = DataEntry{"YMFWAT", UnitSystem::measure::identity, cYmfwat_};
    insert(solutionYfrac);

    allocated_ = false;
}

template<class Scalar>
void CO2H2Container<Scalar>::
readRestart(const unsigned globalDofIdx,
            const unsigned elemIdx,
            const data::Solution& sol)
{
    if (this->allocated_) {
        return;
    }

    auto assign = [elemIdx, globalDofIdx, &sol](const std::string& name,
                                                ScalarBuffer& data)

    {
        if (!data.empty() && sol.has(name)) {
            data[elemIdx] = sol.data<double>(name)[globalDofIdx];
        }
    };

    assign("XMFCO2", *&cXmfCO2_);
    assign("XMFH2", *&cXmfH2_);
    assign("YMFWAT", *&cYmfwat_);
}

template class CO2H2Container<double>;

#if FLOW_INSTANTIATE_FLOAT
template class CO2H2Container<float>;
#endif

} // namespace Opm
