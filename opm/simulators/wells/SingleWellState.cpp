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

#include <opm/simulators/wells/SingleWellState.hpp>

namespace Opm {

SingleWellState::SingleWellState(bool is_producer, std::size_t num_phases, double temp)
    : producer(is_producer)
    , temperature(temp)
    , well_potentials(num_phases)
    , productivity_index(num_phases)
    , surface_rates(num_phases)
{}


void SingleWellState::init_timestep(const SingleWellState& other) {
    if (this->producer != other.producer)
        return;

    if (this->status == Well::Status::SHUT)
        return;

    if (other.status == Well::Status::SHUT)
        return;

    this->bhp = other.bhp;
    this->thp = other.thp;
    this->temperature = other.temperature;
}


void SingleWellState::shut() {
    this->bhp = 0;
    this->thp = 0;
    this->status = Well::Status::SHUT;
    std::fill(this->surface_rates.begin(), this->surface_rates.end(), 0);
    std::fill(this->productivity_index.begin(), this->productivity_index.end(), 0);
}

void SingleWellState::stop() {
    this->thp = 0;
    this->status = Well::Status::STOP;
}

}



