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
#include <opm/simulators/flow/BioeffectsContainer.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <opm/simulators/flow/SolutionContainers.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <tuple>

namespace Opm {

template<class Scalar>
void BioeffectsContainer<Scalar>::
allocate(const unsigned bufferSize, const bool isMICP)
{
    cMicrobes_.resize(bufferSize, 0.0);
    cBiofilm_.resize(bufferSize, 0.0);
    if (isMICP) {
        cOxygen_.resize(bufferSize, 0.0);
        cUrea_.resize(bufferSize, 0.0);
        cCalcite_.resize(bufferSize, 0.0);
    }
    allocated_ = true;
}

template<class Scalar>
void BioeffectsContainer<Scalar>::
assign(const unsigned globalDofIdx,
       const Scalar oxygenConcentration,
       const Scalar ureaConcentration,
       const Scalar calciteVolumeFraction)
{
    cOxygen_[globalDofIdx] = oxygenConcentration;
    cUrea_[globalDofIdx] = ureaConcentration;
    cCalcite_[globalDofIdx] = calciteVolumeFraction;
}

template<class Scalar>
void BioeffectsContainer<Scalar>::
assign(const unsigned globalDofIdx,
       const Scalar microbialConcentration,
       const Scalar biofilmVolumeFraction)
{
    cMicrobes_[globalDofIdx] = microbialConcentration;
    cBiofilm_[globalDofIdx] = biofilmVolumeFraction;
}

template<class Scalar>
BioeffectsSolutionContainer<Scalar>
BioeffectsContainer<Scalar>::
getSolution() const
{
    return {
        cMicrobes_,
        cOxygen_,
        cUrea_,
        cBiofilm_,
        cCalcite_
    };
}

template<class Scalar>
void BioeffectsContainer<Scalar>::
outputRestart(data::Solution& sol, const bool isMICP)
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

    auto solutionMicrobes = DataEntry{"MICROBES",  UnitSystem::measure::density, cMicrobes_};
    insert(solutionMicrobes);
    auto solutionBiofilm = DataEntry{"BIOFILM",  UnitSystem::measure::identity, cBiofilm_};
    insert(solutionBiofilm);

    if (isMICP) {
        auto solutionVectors = std::array {
            DataEntry{"CALCITE",  UnitSystem::measure::identity, cCalcite_},
            DataEntry{"OXYGEN",   UnitSystem::measure::density,  cOxygen_},
            DataEntry{"UREA",     UnitSystem::measure::density,  cUrea_},
        };
        std::for_each(solutionVectors.begin(), solutionVectors.end(),
            [&insert](auto& entry)
            { insert(entry); });
    }

    allocated_ = false;
}

template<class Scalar>
void BioeffectsContainer<Scalar>::
readRestart(const unsigned globalDofIdx,
            const unsigned elemIdx,
            const data::Solution& sol,
            const bool isMICP)
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

    assign("MICROBES", *&cMicrobes_);
    assign("BIOFILM", *&cBiofilm_);

    if (isMICP) {
        const auto fields = std::array{
            std::pair{"OXYGEN",   &cOxygen_},
            std::pair{"UREA",     &cUrea_},
            std::pair{"CALCITE",  &cCalcite_},
        };
        std::for_each(fields.begin(), fields.end(),
            [&assign](const auto& p)
            { assign(p.first, *p.second); });
    }
}

template class BioeffectsContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class BioeffectsContainer<float>;
#endif

} // namespace Opm
