/*
  Copyright 2024 Equinor AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp>

#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>

#include <cstddef>
#include <optional>
#include <unordered_set>
#include <vector>

template <typename Scalar>
std::optional<std::size_t>
Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<Scalar>::
pointID(const int cellIdx) const
{
    const auto& [regPos, inserted] =
        this->seen_.insert((*this->region_)[cellIdx]);

    // inserted == true  => new/unseen region,
    //          == false => previously seen region.
    return inserted
        ? std::optional { static_cast<std::size_t>(*regPos) }
        : std::nullopt;
}

template <typename Scalar>
void
Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<Scalar>::
populateCheckPoint(const int                        cellIdx,
                   EclEpsScalingPointsInfo<Scalar>& endPoints) const
{
    endPoints.extractUnscaled(*this->unscaledEndPoints_.rtep,
                              *this->unscaledEndPoints_.rfunc,
                              (*this->region_)[cellIdx] - this->regIdxOffset_);
}

// ===========================================================================
// Explicit Specialisations
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<float>;
template class Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<double>;
