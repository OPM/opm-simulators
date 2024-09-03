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

#include <opm/simulators/utils/satfunc/PhaseCheckBase.hpp>

namespace {
    constexpr auto one_bit = static_cast<unsigned char>(1);
    constexpr auto failed_bit = one_bit << 0u;
    constexpr auto critical_bit = one_bit << 1u;
}

// ===========================================================================

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::PhaseCheckBase<Scalar>::
test(const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->flags_ = static_cast<unsigned char>(0);

    this->testImpl(endPoints);
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::PhaseCheckBase<Scalar>::isViolated() const
{
    return (this->flags_ & failed_bit) != 0;
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::PhaseCheckBase<Scalar>::isCritical() const
{
    return (this->flags_ & critical_bit) != 0;
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::PhaseCheckBase<Scalar>::setViolated()
{
    this->flags_ |= failed_bit;
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::PhaseCheckBase<Scalar>::setCritical()
{
    this->flags_ |= critical_bit;
}

// ===========================================================================
// Explicit Specialisations of PhaseCheckBase Templates
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::PhaseCheckBase<float>;
template class Opm::Satfunc::PhaseChecks::PhaseCheckBase<double>;
