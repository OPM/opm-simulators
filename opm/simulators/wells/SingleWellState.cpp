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
#include <opm/simulators/wells/PerforationData.hpp>

namespace Opm {

SingleWellState::SingleWellState(const std::string& name_, const ParallelWellInfo& pinfo, bool is_producer, const std::vector<PerforationData>& perf_input, std::size_t num_phases, double temp)
    : name(name_)
    , parallel_info(pinfo)
    , producer(is_producer)
    , temperature(temp)
    , well_potentials(num_phases)
    , productivity_index(num_phases)
    , surface_rates(num_phases)
    , reservoir_rates(num_phases)
    , perf_data(perf_input.size(), !is_producer, num_phases)
{
    for (std::size_t perf = 0; perf < perf_input.size(); perf++) {
        this->perf_data.cell_index[perf] = perf_input[perf].cell_index;
        this->perf_data.connection_transmissibility_factor[perf] = perf_input[perf].connection_transmissibility_factor;
        this->perf_data.satnum_id[perf] = perf_input[perf].satnum_id;
        this->perf_data.ecl_index[perf] = perf_input[perf].ecl_index;
    }
}


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
    std::fill(this->reservoir_rates.begin(), this->reservoir_rates.end(), 0);
    std::fill(this->productivity_index.begin(), this->productivity_index.end(), 0);

    auto& connpi = this->perf_data.prod_index;
    connpi.assign(connpi.size(), 0);
}

void SingleWellState::stop() {
    this->thp = 0;
    this->status = Well::Status::STOP;
}

void SingleWellState::open() {
    this->status = Well::Status::OPEN;
}

void SingleWellState::updateStatus(Well::Status new_status) {
    switch (new_status) {
    case Well::Status::OPEN:
        this->open();
        break;
    case Well::Status::SHUT:
        this->shut();
        break;
    case Well::Status::STOP:
        this->stop();
        break;
    default:
        throw std::logic_error("Invalid well status");
    }
}

double SingleWellState::sum_connection_rates(const std::vector<double>& connection_rates) const {
    return this->parallel_info.get().sumPerfValues(connection_rates.begin(), connection_rates.end());
}

double SingleWellState::sum_brine_rates() const {
    return this->sum_connection_rates(this->perf_data.brine_rates);
}


double SingleWellState::sum_polymer_rates() const {
    return this->sum_connection_rates(this->perf_data.polymer_rates);
}


double SingleWellState::sum_solvent_rates() const {
    return this->sum_connection_rates(this->perf_data.solvent_rates);
}

}



