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
#include <opm/simulators/flow/ExtboContainer.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <tuple>

namespace Opm {

template<class Scalar>
void ExtboContainer<Scalar>::
allocate(const unsigned bufferSize)
{
    X_volume_.resize(bufferSize, 0.0);
    Y_volume_.resize(bufferSize, 0.0);
    Z_fraction_.resize(bufferSize, 0.0);
    mFracOil_.resize(bufferSize, 0.0);
    mFracGas_.resize(bufferSize, 0.0);
    mFracCo2_.resize(bufferSize, 0.0);

    allocated_ = true;
}

template<class Scalar>
void ExtboContainer<Scalar>::
assignMassFractions(const unsigned globalDofIdx,
                    const Scalar gas,
                    const Scalar oil,
                    const Scalar co2)
{
    this->mFracGas_[globalDofIdx] = gas;
    this->mFracOil_[globalDofIdx] = oil;
    this->mFracCo2_[globalDofIdx] = co2;
}

template<class Scalar>
void ExtboContainer<Scalar>::
assignVolumes(const unsigned globalDofIdx,
              const Scalar xVolume,
              const Scalar yVolume)
{
    X_volume_[globalDofIdx] = xVolume;
    Y_volume_[globalDofIdx] = yVolume;
}

template<class Scalar>
void ExtboContainer<Scalar>::
assignZFraction(const unsigned globalDofIdx,
                const Scalar zFraction)
{
    Z_fraction_[globalDofIdx] = zFraction;
}

template<class Scalar>
void ExtboContainer<Scalar>::
outputRestart(data::Solution& sol)
{
    if (!this->allocated_) {
        return;
    }

    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, std::vector<Scalar>&>;

    auto solutionArrays = std::array {
        DataEntry{"SS_X",     UnitSystem::measure::identity,           X_volume_},
        DataEntry{"SS_Y",     UnitSystem::measure::identity,           Y_volume_},
        DataEntry{"SS_Z",     UnitSystem::measure::identity,           Z_fraction_},
        DataEntry{"STD_CO2",  UnitSystem::measure::identity,           mFracCo2_},
        DataEntry{"STD_GAS",  UnitSystem::measure::identity,           mFracGas_},
        DataEntry{"STD_OIL",  UnitSystem::measure::identity,           mFracOil_},
    };

    std::for_each(solutionArrays.begin(), solutionArrays.end(),
                  [&sol](auto& entry)
                  {
                      if (!std::get<2>(entry).empty()) {
                          sol.insert(std::get<std::string>(entry),
                          std::get<UnitSystem::measure>(entry),
                          std::move(std::get<2>(entry)),
                          data::TargetType::RESTART_OPM_EXTENDED);
                      }
                  });

    this->allocated_ = false;
}

template class ExtboContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class ExtboContainer<float>;
#endif

} // namespace Opm
