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

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <array>
#include <tuple>

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

    this->dispX_.resize(bufferSize, 0.0);
    rstKeywords["DISPX"] = 0;
    this->dispY_.resize(bufferSize, 0.0);
    rstKeywords["DISPY"] = 0;
    this->dispZ_.resize(bufferSize, 0.0);
    rstKeywords["DISPZ"] = 0;
    this->stressXX_.resize(bufferSize, 0.0);
    rstKeywords["STRESSXX"] = 0;
    this->stressYY_.resize(bufferSize, 0.0);
    rstKeywords["STRESSYY"] = 0;
    this->stressZZ_.resize(bufferSize, 0.0);
    rstKeywords["STRESSZZ"] = 0;
    this->stressXY_.resize(bufferSize, 0.0);
    rstKeywords["STRESSXY"] = 0;
    this->stressXZ_.resize(bufferSize, 0.0);
    rstKeywords["STRESSXZ"] = 0;
    this->stressYZ_.resize(bufferSize, 0.0);
    rstKeywords["STRESSYZ"] = 0;

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
    this->dispX_[globalDofIdx] = disp[0];
    this->dispY_[globalDofIdx] = disp[1];
    this->dispZ_[globalDofIdx] = disp[2];
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
outputRestart(data::Solution& sol) const
{
    if (!allocated_) {
        return;
    }
    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, const std::vector<Scalar>&>;

    auto doInsert = [&sol](const DataEntry&       entry,
                           const data::TargetType target)
    {
        if (std::get<2>(entry).empty()) {
            return;
        }

        sol.insert(std::get<std::string>(entry),
                   std::get<UnitSystem::measure>(entry),
                   std::move(std::get<2>(entry)),
                   target);
    };

    const auto solutionVectors = std::array{
        DataEntry{"DELSTRXX", UnitSystem::measure::pressure, delstressXX_},
        DataEntry{"DELSTRYY", UnitSystem::measure::pressure, delstressYY_},
        DataEntry{"DELSTRZZ", UnitSystem::measure::pressure, delstressZZ_},
        DataEntry{"DELSTRXY", UnitSystem::measure::pressure, delstressXY_},
        DataEntry{"DELSTRXZ", UnitSystem::measure::pressure, delstressXZ_},
        DataEntry{"DELSTRYZ", UnitSystem::measure::pressure, delstressYZ_},
        DataEntry{"DISPX",    UnitSystem::measure::length,   dispX_},
        DataEntry{"DISPY",    UnitSystem::measure::length,   dispY_},
        DataEntry{"DISPZ",    UnitSystem::measure::length,   dispZ_},
        DataEntry{"FRCSTRXX", UnitSystem::measure::pressure, fracstressXX_},
        DataEntry{"FRCSTRYY", UnitSystem::measure::pressure, fracstressYY_},
        DataEntry{"FRCSTRZZ", UnitSystem::measure::pressure, fracstressZZ_},
        DataEntry{"FRCSTRXY", UnitSystem::measure::pressure, fracstressXY_},
        DataEntry{"FRCSTRXZ", UnitSystem::measure::pressure, fracstressXZ_},
        DataEntry{"FRCSTRYZ", UnitSystem::measure::pressure, fracstressYZ_},
        DataEntry{"LINSTRXX", UnitSystem::measure::pressure, linstressXX_},
        DataEntry{"LINSTRYY", UnitSystem::measure::pressure, linstressYY_},
        DataEntry{"LINSTRZZ", UnitSystem::measure::pressure, linstressZZ_},
        DataEntry{"LINSTRXY", UnitSystem::measure::pressure, linstressXY_},
        DataEntry{"LINSTRXZ", UnitSystem::measure::pressure, linstressXZ_},
        DataEntry{"LINSTRYZ", UnitSystem::measure::pressure, linstressYZ_},
        DataEntry{"MECHPOTF", UnitSystem::measure::pressure, potentialForce_},
        DataEntry{"PRESPOTF", UnitSystem::measure::pressure, potentialPressForce_},
        DataEntry{"STRAINXX", UnitSystem::measure::identity, strainXX_},
        DataEntry{"STRAINYY", UnitSystem::measure::identity, strainYY_},
        DataEntry{"STRAINZZ", UnitSystem::measure::identity, strainZZ_},
        DataEntry{"STRAINXY", UnitSystem::measure::identity, strainXY_},
        DataEntry{"STRAINXZ", UnitSystem::measure::identity, strainXZ_},
        DataEntry{"STRAINYZ", UnitSystem::measure::identity, strainYZ_},
        DataEntry{"STRESSXX", UnitSystem::measure::length,   stressXX_},
        DataEntry{"STRESSYY", UnitSystem::measure::length,   stressYY_},
        DataEntry{"STRESSZZ", UnitSystem::measure::length,   stressZZ_},
        DataEntry{"STRESSXY", UnitSystem::measure::length,   stressXY_},
        DataEntry{"STRESSXZ", UnitSystem::measure::length,   stressXZ_},
        DataEntry{"STRESSYZ", UnitSystem::measure::length,   stressYZ_},
        DataEntry{"TEMPPOTF", UnitSystem::measure::pressure, potentialTempForce_},
    };

    std::for_each(solutionVectors.begin(), solutionVectors.end(),
                  [doInsert](const auto& array)
                  { doInsert(array, data::TargetType::RESTART_OPM_EXTENDED); });
}

template class MechContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class MechContainer<float>;
#endif

} // namespace Opm
