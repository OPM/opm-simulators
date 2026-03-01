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

#include <opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.hpp>

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sogcr_ = endPoints.Sogcr;

    if (! std::isfinite(this->sogcr_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->sogcr_ < Scalar{0} - this->epsilon_;
    const auto high = !(this->sogcr_ < Scalar{1} + this->epsilon_);

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->swl_ = endPoints.Swl;
    this->sgu_ = endPoints.Sgu;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (this->swl_ + this->sgu_ > Scalar{1} + this->epsilon_) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgl_   = endPoints.Sgl;
    this->sogcr_ = endPoints.Sogcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sogcr_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (! (this->sogcr_ < Scalar{1} - this->swl_ - this->sgl_ + this->epsilon_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgcr_  = endPoints.Sgcr;
    this->sogcr_ = endPoints.Sogcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgcr_) ||
        ! std::isfinite(this->sogcr_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (! (this->sogcr_ < Scalar{1} - this->swl_ - this->sgcr_ + this->epsilon_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->sowcr_)) {
        this->setViolated();
        this->setCritical();

        return;
    }

    const auto low = this->sowcr_ < Scalar{0} - this->epsilon_;
    const auto high = ! (this->sowcr_ < Scalar{1} + this->epsilon_);

    if (low || high) {
        this->setViolated();
        this->setCritical();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sgl_ = endPoints.Sgl;
    this->swu_ = endPoints.Swu;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->swu_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (this->sgl_ + this->swu_ > Scalar{1} + this->epsilon_) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgl_   = endPoints.Sgl;
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sowcr_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (! (this->sowcr_ < Scalar{1} - this->swl_ - this->sgl_ + this->epsilon_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWcr<Scalar>::
testImpl(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->sgl_   = endPoints.Sgl;
    this->swcr_  = endPoints.Swcr;
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->swcr_) ||
        ! std::isfinite(this->sowcr_))
    {
        this->setViolated();
        this->setCritical();

        return;
    }

    if (! (this->sowcr_ < Scalar{1} - this->swcr_ - this->sgl_ + this->epsilon_)) {
        this->setViolated();
    }
}

// ===========================================================================
// Explicit Specialisations of Individual Check Templates
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWcr<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWcr<double>;
