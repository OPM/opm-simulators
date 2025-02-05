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
#include <opm/simulators/flow/MechContainer.hpp>

#include <opm/common/utility/Visitor.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <array>
#include <tuple>
#include <variant>

namespace Opm {

template<class Scalar>
void MechContainer<Scalar>::
allocate(const std::size_t bufferSize,
         std::map<std::string, int>& rstKeywords)
{
    this->potentialForce_.resize(bufferSize, 0.0);
    rstKeywords["MECHPOTF"] = 0;
    this->potentialTempForce_.resize(bufferSize, 0.0);
    rstKeywords["TEMPPOTF"] = 0;
    this->potentialPressForce_.resize(bufferSize, 0.0);
    rstKeywords["PRESPOTF"] = 0;

    std::for_each(disp_.begin(), disp_.end(),
                  [suffix = "XYZ", is = 0, bufferSize, &rstKeywords](auto& v) mutable
                  {
                      v.resize(bufferSize, 0.0);
                      rstKeywords[std::string{"DISP"} + suffix[is++]] = 0;
                  });

    auto resizeAndRegister =
        [&rstKeywords,bufferSize](auto& tensor,
                                  const std::string& name)
    {
        static constexpr auto fields = std::array{
            "XX", "YY", "ZZ", "YZ", "XZ", "XY",
        };
        tensor.resize(bufferSize);
        for (const auto& f : fields) {
            rstKeywords[name + f] = 0;
        }
    };

    resizeAndRegister(stress_, "STRESS");

    this->strainXX_.resize(bufferSize, 0.0);
    rstKeywords["STRAINXX"] = 0;
    this->strainYY_.resize(bufferSize, 0.0);
    rstKeywords["STRAINYY"] = 0;
    this->strainZZ_.resize(bufferSize, 0.0);
    rstKeywords["STRAINZZ"] = 0;
    this->strainXY_.resize(bufferSize, 0.0);
    rstKeywords["STRAINXY"] = 0;
    this->strainXZ_.resize(bufferSize, 0.0);
    rstKeywords["STRAINXZ"] = 0;
    this->strainYZ_.resize(bufferSize, 0.0);
    rstKeywords["STRAINYZ"] = 0;

    this->delstressXX_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRXX"] = 0;
    this->delstressYY_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRYY"] = 0;
    this->delstressZZ_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRZZ"] = 0;
    this->delstressXY_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRXY"] = 0;
    this->delstressXZ_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRXZ"] = 0;
    this->delstressYZ_.resize(bufferSize, 0.0);
    rstKeywords["DELSTRYZ"] = 0;

    this->fracstressXX_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRXX"] = 0;
    this->fracstressYY_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRYY"] = 0;
    this->fracstressZZ_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRZZ"] = 0;
    this->fracstressXY_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRXY"] = 0;
    this->fracstressXZ_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRXZ"] = 0;
    this->fracstressYZ_.resize(bufferSize,0.0);
    rstKeywords["FRCSTRYZ"] = 0;

    this->linstressXX_.resize(bufferSize,0.0);
    rstKeywords["LINSTRXX"] = 0;
    this->linstressYY_.resize(bufferSize,0.0);
    rstKeywords["LINSTRYY"] = 0;
    this->linstressZZ_.resize(bufferSize,0.0);
    rstKeywords["LINSTRZZ"] = 0;
    this->linstressXY_.resize(bufferSize,0.0);
    rstKeywords["LINSTRXY"] = 0;
    this->linstressXZ_.resize(bufferSize,0.0);
    rstKeywords["LINSTRXZ"] = 0;
    this->linstressYZ_.resize(bufferSize,0.0);
    rstKeywords["LINSTRYZ"] = 0;

    allocated_ = true;
}

template<class Scalar>
void MechContainer<Scalar>::
assignDisplacement(const unsigned globalDofIdx,
                   const Dune::FieldVector<Scalar,3>& disp)
{
    this->disp_[0][globalDofIdx] = disp[0];
    this->disp_[1][globalDofIdx] = disp[1];
    this->disp_[2][globalDofIdx] = disp[2];
}

template<class Scalar>
void MechContainer<Scalar>::
assignPotentialForces(const unsigned globalDofIdx,
                      const Scalar force,
                      const Scalar pressForce,
                      const Scalar tempForce)
{
    potentialForce_[globalDofIdx] = force;
    potentialPressForce_[globalDofIdx] = pressForce;
    potentialTempForce_[globalDofIdx] = tempForce;
}

template<class Scalar>
void MechContainer<Scalar>::
assignStress(const unsigned globalDofIdx,
             const Dune::FieldVector<Scalar,6>& stress)
{
    this->stress_.assign(globalDofIdx, VoigtContainer<Scalar>(stress));
}

template<class Scalar>
void MechContainer<Scalar>::
outputRestart(data::Solution& sol) const
{
    if (!allocated_) {
        return;
    }
    using DataEntry = std::tuple<std::string,
                                 UnitSystem::measure,
                                 std::variant<const std::vector<Scalar>*,
                                              const std::array<std::vector<Scalar>,3>*,
                                              const VoigtArray<Scalar>*>>;

    auto doInsert = [&sol](const std::string& name,
                           const UnitSystem::measure& measure,
                           const std::vector<Scalar>& entry)
    {
        if (!entry.empty()) {
            sol.insert(name, measure, entry, data::TargetType::RESTART_OPM_EXTENDED);
        }
    };

    const auto solutionVectors = std::array{
        DataEntry{"DELSTRXX", UnitSystem::measure::pressure, &delstressXX_},
        DataEntry{"DELSTRYY", UnitSystem::measure::pressure, &delstressYY_},
        DataEntry{"DELSTRZZ", UnitSystem::measure::pressure, &delstressZZ_},
        DataEntry{"DELSTRXY", UnitSystem::measure::pressure, &delstressXY_},
        DataEntry{"DELSTRXZ", UnitSystem::measure::pressure, &delstressXZ_},
        DataEntry{"DELSTRYZ", UnitSystem::measure::pressure, &delstressYZ_},
        DataEntry{"DISP",     UnitSystem::measure::length,   &disp_},
        DataEntry{"FRCSTRXX", UnitSystem::measure::pressure, &fracstressXX_},
        DataEntry{"FRCSTRYY", UnitSystem::measure::pressure, &fracstressYY_},
        DataEntry{"FRCSTRZZ", UnitSystem::measure::pressure, &fracstressZZ_},
        DataEntry{"FRCSTRXY", UnitSystem::measure::pressure, &fracstressXY_},
        DataEntry{"FRCSTRXZ", UnitSystem::measure::pressure, &fracstressXZ_},
        DataEntry{"FRCSTRYZ", UnitSystem::measure::pressure, &fracstressYZ_},
        DataEntry{"LINSTRXX", UnitSystem::measure::pressure, &linstressXX_},
        DataEntry{"LINSTRYY", UnitSystem::measure::pressure, &linstressYY_},
        DataEntry{"LINSTRZZ", UnitSystem::measure::pressure, &linstressZZ_},
        DataEntry{"LINSTRXY", UnitSystem::measure::pressure, &linstressXY_},
        DataEntry{"LINSTRXZ", UnitSystem::measure::pressure, &linstressXZ_},
        DataEntry{"LINSTRYZ", UnitSystem::measure::pressure, &linstressYZ_},
        DataEntry{"MECHPOTF", UnitSystem::measure::pressure, &potentialForce_},
        DataEntry{"PRESPOTF", UnitSystem::measure::pressure, &potentialPressForce_},
        DataEntry{"STRAINXX", UnitSystem::measure::identity, &strainXX_},
        DataEntry{"STRAINYY", UnitSystem::measure::identity, &strainYY_},
        DataEntry{"STRAINZZ", UnitSystem::measure::identity, &strainZZ_},
        DataEntry{"STRAINXY", UnitSystem::measure::identity, &strainXY_},
        DataEntry{"STRAINXZ", UnitSystem::measure::identity, &strainXZ_},
        DataEntry{"STRAINYZ", UnitSystem::measure::identity, &strainYZ_},
        DataEntry{"STRESS",   UnitSystem::measure::length,   &stress_},
        DataEntry{"TEMPPOTF", UnitSystem::measure::pressure, &potentialTempForce_},
    };

    std::for_each(solutionVectors.begin(), solutionVectors.end(),
                  [&doInsert](const auto& array)
                  {
                      std::visit(VisitorOverloadSet{
                                    [&array, &doInsert](const std::vector<Scalar>* v)
                                    {
                                        doInsert(std::get<0>(array), std::get<1>(array), *v);
                                    },
                                    [&array, &doInsert](const std::array<std::vector<Scalar>,3>* V)
                                    {
                                        const auto& v = *V;
                                        const auto& name = std::get<0>(array);
                                        const auto& measure = std::get<1>(array);
                                        doInsert(name + "X", measure, v[0]);
                                        doInsert(name + "Y", measure, v[1]);
                                        doInsert(name + "Z", measure, v[2]);
                                    },
                                    [&array, &doInsert](const VoigtArray<Scalar>* V)
                                    {
                                        const auto& v = *V;
                                        const auto& name = std::get<0>(array);
                                        const auto& measure = std::get<1>(array);
                                        doInsert(name + "XX", measure, v[VoigtIndex::XX]);
                                        doInsert(name + "YY", measure, v[VoigtIndex::YY]);
                                        doInsert(name + "ZZ", measure, v[VoigtIndex::ZZ]);
                                        doInsert(name + "YZ", measure, v[VoigtIndex::YZ]);
                                        doInsert(name + "XZ", measure, v[VoigtIndex::XZ]);
                                        doInsert(name + "XY", measure, v[VoigtIndex::XY]);
                                    }
                                 }, std::get<2>(array));
                  }
    );
}

template class MechContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class MechContainer<float>;
#endif

} // namespace Opm
