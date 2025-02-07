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

    resizeAndRegister(delstress_, "DELSTR");
    resizeAndRegister(fracstress_, "FRCSTR");
    resizeAndRegister(linstress_, "LINSTR");
    resizeAndRegister(strain_, "STRAIN");
    resizeAndRegister(stress_, "STRESS");

    allocated_ = true;
}

template<class Scalar>
void MechContainer<Scalar>::
assignDelStress(const unsigned globalDofIdx,
                const Dune::FieldVector<Scalar,6>& delStress)
{
    this->delstress_.assign(globalDofIdx, VoigtContainer<Scalar>(delStress));
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
assignFracStress(const unsigned globalDofIdx,
                 const Dune::FieldVector<Scalar,6>& fracStress)
{
    this->fracstress_.assign(globalDofIdx, VoigtContainer<Scalar>(fracStress));
}

template<class Scalar>
void MechContainer<Scalar>::
assignLinStress(const unsigned globalDofIdx,
                const Dune::FieldVector<Scalar,6>& linStress)
{
    this->linstress_.assign(globalDofIdx, VoigtContainer<Scalar>(linStress));
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
assignStrain(const unsigned globalDofIdx,
             const Dune::FieldVector<Scalar,6>& strain)
{
    this->strain_.assign(globalDofIdx, VoigtContainer<Scalar>(strain));
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
outputRestart(data::Solution& sol)
{
    if (!allocated_) {
        return;
    }
    using DataEntry = std::tuple<std::string,
                                 UnitSystem::measure,
                                 std::variant<std::vector<Scalar>*,
                                              std::array<std::vector<Scalar>,3>*,
                                              VoigtArray<Scalar>*>>;

    auto doInsert = [&sol](const std::string& name,
                           const UnitSystem::measure& measure,
                           std::vector<Scalar>& entry)
    {
        if (!entry.empty()) {
            sol.insert(name, measure, std::move(entry), data::TargetType::RESTART_OPM_EXTENDED);
        }
    };

    const auto solutionVectors = std::array{
        DataEntry{"DELSTR",   UnitSystem::measure::pressure, &delstress_},
        DataEntry{"DISP",     UnitSystem::measure::length,   &disp_},
        DataEntry{"FRCSTR",   UnitSystem::measure::pressure, &fracstress_},
        DataEntry{"LINSTR",   UnitSystem::measure::pressure, &linstress_},
        DataEntry{"MECHPOTF", UnitSystem::measure::pressure, &potentialForce_},
        DataEntry{"PRESPOTF", UnitSystem::measure::pressure, &potentialPressForce_},
        DataEntry{"STRAIN",   UnitSystem::measure::identity, &strain_},
        DataEntry{"STRESS",   UnitSystem::measure::length,   &stress_},
        DataEntry{"TEMPPOTF", UnitSystem::measure::pressure, &potentialTempForce_},
    };

    std::for_each(solutionVectors.begin(), solutionVectors.end(),
                  [&doInsert](auto& array)
                  {
                      std::visit(VisitorOverloadSet{
                                    [&array, &doInsert](std::vector<Scalar>* v)
                                    {
                                        doInsert(std::get<0>(array), std::get<1>(array), *v);
                                    },
                                    [&array, &doInsert](std::array<std::vector<Scalar>,3>* V)
                                    {
                                        auto& v = *V;
                                        const auto& name = std::get<0>(array);
                                        const auto& measure = std::get<1>(array);
                                        doInsert(name + "X", measure, v[0]);
                                        doInsert(name + "Y", measure, v[1]);
                                        doInsert(name + "Z", measure, v[2]);
                                    },
                                    [&array, &doInsert](VoigtArray<Scalar>* V)
                                    {
                                        auto& v = *V;
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

    allocated_ = false;
}

template class MechContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class MechContainer<float>;
#endif

} // namespace Opm
