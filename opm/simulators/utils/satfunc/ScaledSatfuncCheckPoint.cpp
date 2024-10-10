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

#include <opm/simulators/utils/satfunc/ScaledSatfuncCheckPoint.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>
#include <opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp>

template <typename Scalar>
void
Opm::Satfunc::PhaseChecks::ScaledSatfuncCheckPoint<Scalar>::
populateCheckPoint(const int                        cellIdx,
                   EclEpsScalingPointsInfo<Scalar>& endPoints) const
{
    this->unscaled_.populateCheckPoint(cellIdx, endPoints);

    endPoints.extractScaled(*this->eclipseState_,
                            *this->epsGridProps_, cellIdx);
}

// ===========================================================================
// Explicit Specialisations
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::ScaledSatfuncCheckPoint<float>;
template class Opm::Satfunc::PhaseChecks::ScaledSatfuncCheckPoint<double>;
