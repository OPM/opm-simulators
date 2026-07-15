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

#include <opm/simulators/utils/satfunc/WaterPhaseConsistencyChecks.hpp>

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Water::SWmin<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // 0 <= SWL < 1

    this->swl_ = endPoints.Swl;

    if (! std::isfinite(this->swl_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->swl_ < Scalar{0} - this->epsilon_;
    const auto high = ! (this->swl_ < Scalar{1} + this->epsilon_);

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Water::SWmax<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // 0 < SWU <= 1

    this->swu_ = endPoints.Swu;

    if (! std::isfinite(this->swu_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = ! (this->swu_ > Scalar{0} - this->epsilon_);
    const auto high = this->swu_ > Scalar{1} + this->epsilon_;

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Water::SWcr<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // SWL <= SWCR < SWU

    this->swl_  = endPoints.Swl;
    this->swcr_ = endPoints.Swcr;
    this->swu_  = endPoints.Swu;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->swcr_) ||
        ! std::isfinite(this->swu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->swcr_ < this->swl_ - this->epsilon_;
    const auto high = ! (this->swcr_ < this->swu_ + this->epsilon_);

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

template class Opm::Satfunc::PhaseChecks::Water::SWmin<float>;
template class Opm::Satfunc::PhaseChecks::Water::SWmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Water::SWmax<float>;
template class Opm::Satfunc::PhaseChecks::Water::SWmax<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Water::SWcr<float>;
template class Opm::Satfunc::PhaseChecks::Water::SWcr<double>;
