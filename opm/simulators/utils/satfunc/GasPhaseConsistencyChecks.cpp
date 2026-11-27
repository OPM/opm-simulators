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

#include <opm/simulators/utils/satfunc/GasPhaseConsistencyChecks.hpp>

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Gas::SGmin<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sgl_ = endPoints.Sgl;

    if (! std::isfinite(this->sgl_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->sgl_ < Scalar{0} - this->epsilon_;
    const auto high = ! (this->sgl_ < Scalar{1} + this->epsilon_);

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Gas::SGmax<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    // 0 < SGU <= 1

    this->sgu_ = endPoints.Sgu;

    if (! std::isfinite(this->sgu_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = ! (this->sgu_ > Scalar{0} - this->epsilon_);
    const auto high = this->sgu_ > Scalar{1} + this->epsilon_;

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Gas::SGcr<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sgl_  = endPoints.Sgl;
    this->sgcr_ = endPoints.Sgcr;
    this->sgu_  = endPoints.Sgu;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sgcr_) ||
        ! std::isfinite(this->sgu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->sgcr_ < this->sgl_ - this->epsilon_;
    const auto high = ! (this->sgcr_ < this->sgu_ + this->epsilon_);

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

template class Opm::Satfunc::PhaseChecks::Gas::SGmin<float>;
template class Opm::Satfunc::PhaseChecks::Gas::SGmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Gas::SGmax<float>;
template class Opm::Satfunc::PhaseChecks::Gas::SGmax<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Gas::SGcr<float>;
template class Opm::Satfunc::PhaseChecks::Gas::SGcr<double>;
