/*
  Copyright 2021 Equinor ASA

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/wells/ALQState.hpp>

#include <cstddef>
#include <stdexcept>

namespace Opm {

template<class Scalar>
ALQState<Scalar> ALQState<Scalar>::serializationTestObject()
{
    ALQState result;
    result.current_alq_ = {1.0};
    result.default_alq_ = {2.0};
    result.alq_increase_count_= {4};
    result.alq_decrease_count_= {5};
    return result;
}

template<class Scalar>
Scalar ALQState<Scalar>::get() const
{
    if (this->current_alq_)
        return *this->current_alq_;

    return this->default_alq_;
}

template<class Scalar>
void ALQState<Scalar>::update_default(Scalar value)
{
    if (default_alq_ != value) {
        this->default_alq_ = value;
        this->current_alq_ = value;
    }
}

template<class Scalar>
void ALQState<Scalar>::set(Scalar value)
{
    this->current_alq_ = value;
}

template<class Scalar>
bool ALQState<Scalar>::oscillation() const
{
    return this->alq_increase_count_ > 0 && this->alq_decrease_count_ >= 1;
}

template<class Scalar>
void ALQState<Scalar>::update_count(bool increase)
{
    if (increase)
        this->alq_increase_count_ += 1;
    else
        this->alq_decrease_count_ += 1;
}

template<class Scalar>
void ALQState<Scalar>::reset_count()
{
    this->alq_decrease_count_ = 0;
    this->alq_increase_count_ = 0;
}

template<class Scalar>
int ALQState<Scalar>::get_increment_count() const
{
    return this->alq_increase_count_;
}

template<class Scalar>
int ALQState<Scalar>::get_decrement_count() const
{
    return this->alq_decrease_count_;
}

template<class Scalar>
bool ALQState<Scalar>::operator==(const ALQState& rhs) const
{
    return this->current_alq_ == rhs.current_alq_ &&
           this->default_alq_ == rhs.default_alq_ &&
           this->alq_increase_count_ == rhs.alq_increase_count_ &&
           this->alq_decrease_count_ == rhs.alq_decrease_count_;
}

template class ALQState<double>;

#if FLOW_INSTANTIATE_FLOAT
template class ALQState<float>;
#endif

}
