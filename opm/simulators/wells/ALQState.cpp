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
    result.current_alq_ = {{"test1", 1.0}};
    result.default_alq_ = {{"test2", 2.0}, {"test3", 3.0}};
    result.alq_increase_count_= {{"test4", 4}};
    result.alq_decrease_count_= {{"test5", 5}};
    result.debug_counter_ = 6;

    return result;
}

template<class Scalar>
Scalar ALQState<Scalar>::get(const std::string& wname) const
{
    auto iter = this->current_alq_.find(wname);
    if (iter != this->current_alq_.end())
        return iter->second;

    auto default_iter = this->default_alq_.find(wname);
    if (default_iter != this->default_alq_.end())
        return default_iter->second;

    throw std::logic_error("No ALQ value registered for well: " + wname);
}

template<class Scalar>
void ALQState<Scalar>::update_default(const std::string& wname, Scalar value)
{
    auto default_iter = this->default_alq_.find(wname);
    if (default_iter == this->default_alq_.end() || default_iter->second != value) {
        this->default_alq_.insert_or_assign(wname, value);
        this->current_alq_.insert_or_assign(wname, value);
    }
}

template<class Scalar>
void ALQState<Scalar>::set(const std::string& wname, Scalar value)
{
    this->current_alq_[wname] = value;
}

template<class Scalar>
int ALQState<Scalar>::get_debug_counter()
{
    return this->debug_counter_;
}

template<class Scalar>
int ALQState<Scalar>::update_debug_counter()
{
    this->debug_counter_++;
    return this->debug_counter_;
}

template<class Scalar>
void ALQState<Scalar>::set_debug_counter(int value)
{
    this->debug_counter_ = value;
}

namespace {

int get_counter(const std::map<std::string, int>& count_map, const std::string& wname) {
    auto count_iter = count_map.find(wname);
    if (count_iter == count_map.end())
        return 0;
    return count_iter->second;
}

}

template<class Scalar>
bool ALQState<Scalar>::oscillation(const std::string& wname) const
{
    auto inc_count = get_counter(this->alq_increase_count_, wname);
    if (inc_count == 0)
        return false;

    auto dec_count = get_counter(this->alq_decrease_count_, wname);
    return dec_count >= 1;
}

template<class Scalar>
void ALQState<Scalar>::update_count(const std::string& wname, bool increase)
{
    if (increase)
        this->alq_increase_count_[wname] += 1;
    else
        this->alq_decrease_count_[wname] += 1;
}

template<class Scalar>
void ALQState<Scalar>::reset_count()
{
    this->alq_decrease_count_.clear();
    this->alq_increase_count_.clear();
}

template<class Scalar>
int ALQState<Scalar>::get_increment_count(const std::string& wname) const
{
    return get_counter(this->alq_increase_count_, wname);
}

template<class Scalar>
int ALQState<Scalar>::get_decrement_count(const std::string& wname) const
{
    return get_counter(this->alq_decrease_count_, wname);
}

template<class Scalar>
std::size_t ALQState<Scalar>::pack_size() const
{
    return this->current_alq_.size();
}

template<class Scalar>
std::size_t ALQState<Scalar>::pack_data(Scalar* data) const
{
    std::size_t index = 0;
    for (const auto& [_, value] : this->current_alq_) {
        (void)_;
        data[index++] = value;
    }
    return index;
}

template<class Scalar>
std::size_t ALQState<Scalar>::unpack_data(const Scalar* data)
{
    std::size_t index = 0;
    for (auto& [_, value] : this->current_alq_) {
        (void)_;
        value = data[index++];
    }
    return index;
}

template<class Scalar>
bool ALQState<Scalar>::operator==(const ALQState& rhs) const
{
    return this->current_alq_ == rhs.current_alq_ &&
           this->default_alq_ == rhs.default_alq_ &&
           this->alq_increase_count_ == rhs.alq_increase_count_ &&
           this->alq_decrease_count_ == rhs.alq_decrease_count_ &&
           this->debug_counter_ == rhs.debug_counter_;
}

template class ALQState<double>;

}
