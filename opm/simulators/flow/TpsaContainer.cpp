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

#include <opm/simulators/flow/TpsaContainer.hpp>

#include <opm/common/utility/Visitor.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Solution.hpp>

namespace Opm
{

template <class Scalar>
void
TpsaContainer<Scalar>::
allocate(const std::size_t bufferSize,
         std::map<std::string, int>& rstKeywords)
{
    this->solidPressure_.resize(bufferSize, 0.0);
    rstKeywords["SPRES"] = 0;

    std::ranges::for_each(rotation_,
                          [suffix = 'X', bufferSize, &rstKeywords](auto& v) mutable {
                              v.resize(bufferSize, 0.0);
                              rstKeywords[std::string{"ROT"} + suffix++] = 0;
                          });

    allocated_ = true;
}

template <class Scalar>
void
TpsaContainer<Scalar>::
assignRotation(const unsigned globalDofIdx,
               const Dune::FieldVector<Scalar, 3>& rotation)
{
    this->rotation_[0][globalDofIdx] = rotation[0];
    this->rotation_[1][globalDofIdx] = rotation[1];
    this->rotation_[2][globalDofIdx] = rotation[2];
}

template <class Scalar>
void
TpsaContainer<Scalar>::
assignSolidPressure(const unsigned globalDofIdx,
                    const Scalar solidPressure)
{
    solidPressure_[globalDofIdx] = solidPressure;
}

template <class Scalar>
void
TpsaContainer<Scalar>::
outputRestart(data::Solution& sol)
{
    if (!allocated_) {
        return;
    }
    using DataEntry = std::tuple<std::string,
                                 UnitSystem::measure,
                                 std::variant<std::vector<Scalar>*,
                                              std::array<std::vector<Scalar>, 3>*> >;

    auto doInsert = [&sol](const std::string& name,
                           const UnitSystem::measure& measure,
                           std::vector<Scalar>& entry) {
        if (!entry.empty()) {
            sol.insert(name, measure, std::move(entry), data::TargetType::RESTART_OPM_EXTENDED);
        }
    };

    const auto solutionVectors = std::array{
        DataEntry{"SPRES", UnitSystem::measure::pressure, &solidPressure_},
        DataEntry{"ROT", UnitSystem::measure::pressure, &rotation_},
    };

    std::ranges::for_each(solutionVectors,
                          [&doInsert](auto& array) {
                              std::visit(VisitorOverloadSet{
                                             [&array, &doInsert](std::vector<Scalar>* v) {
                                                 doInsert(std::get<0>(array),
                                                          std::get<1>(array),
                                                          *v);
                                             },
                                             [&array, &doInsert](
                                             std::array<std::vector<Scalar>, 3>* V) {
                                                 auto& v = *V;
                                                 const auto& name = std::get<0>(array);
                                                 const auto& measure = std::get<1>(array);
                                                 doInsert(name + "X", measure, v[0]);
                                                 doInsert(name + "Y", measure, v[1]);
                                                 doInsert(name + "Z", measure, v[2]);
                                             },
                                         },
                                         std::get<2>(array));
                          }
        );
    allocated_ = false;
}

template class TpsaContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class TpsaContainer<float>;
#endif

} // namespace Opm
