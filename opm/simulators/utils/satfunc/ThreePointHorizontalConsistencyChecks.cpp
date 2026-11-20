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

#include <opm/simulators/utils/satfunc/ThreePointHorizontalConsistencyChecks.hpp>

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_GO<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // SGCR < 1-SOGCR-SWL <= SGU

    this->swl_   = endPoints.Swl;
    this->sogcr_ = endPoints.Sogcr;
    this->sgcr_  = endPoints.Sgcr;
    this->sgu_   = endPoints.Sgu;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sogcr_) ||
        ! std::isfinite(this->sgcr_) ||
        ! std::isfinite(this->sgu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto sr = Scalar{1} - (this->sogcr_ + this->swl_);

    const auto low = ! (this->sgcr_ < sr + this->epsilon_);
    const auto high = sr > this->sgu_ + this->epsilon_;

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_OW<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // SWCR < 1-SOWCR-SGL <= SWU

    this->sgl_   = endPoints.Sgl;
    this->sowcr_ = endPoints.Sowcr;
    this->swcr_  = endPoints.Swcr;
    this->swu_   = endPoints.Swu;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sowcr_) ||
        ! std::isfinite(this->swcr_) ||
        ! std::isfinite(this->swu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto sr = Scalar{1} - (this->sowcr_ + this->sgl_);

    const auto low = ! (this->swcr_ < sr + this->epsilon_);
    const auto high = sr > this->swu_ + this->epsilon_;

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ===========================================================================
// Explicit Specialisations of Individual Check Templates
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_GO<float>;
template class Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_GO<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_OW<float>;
template class Opm::Satfunc::PhaseChecks::ThreePointHorizontal::DisplacingOil_OW<double>;
