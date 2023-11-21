/*
  Copyright 2021 Equinor ASA.


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

#include <opm/simulators/wells/ConnFiltrateData.hpp>
#include <opm/simulators/wells/PerfData.hpp>

namespace Opm {

PerfData::PerfData(std::size_t num_perf, double pressure_first_connection_, bool injector_, std::size_t num_phases)
    : injector(injector_)
    , pressure_first_connection(pressure_first_connection_)
    , pressure(num_perf)
    , rates(num_perf)
    , phase_rates(num_perf * num_phases)
    , phase_mixing_rates(num_perf)
    , solvent_rates(num_perf)
    , polymer_rates(num_perf)
    , brine_rates(num_perf)
    , prod_index(num_perf * num_phases)
    , micp_rates(num_perf)
    , cell_index(num_perf)
    , connection_transmissibility_factor(num_perf)
    , connection_d_factor(num_perf)
    , satnum_id(num_perf)
    , ecl_index(num_perf)
{
    if (injector) {
        this->water_throughput.resize(num_perf);
        this->skin_pressure.resize(num_perf);
        this->water_velocity.resize(num_perf);
        this->filtrate_data.resize(num_perf);
    }
}

PerfData PerfData::serializationTestObject()
{
    PerfData result;
    result.pressure_first_connection = 1.0;
    result.pressure = {2.0, 3.0, 4.0};
    result.rates = {5.0, 6.0};
    result.phase_rates = {7.0};
    result.phase_mixing_rates = { {1.0, 2.0, 3.0, 4.0}};
    result.solvent_rates = {8.0, 9.0};
    result.polymer_rates = {10.0, 11.0, 12.0};
    result.brine_rates = {13.0};
    result.prod_index = {14.0, 15.0};
    result.micp_rates = {16.0};
    result.cell_index = {17, 18, 19, 20};
    result.connection_transmissibility_factor = {21.0};
    result.connection_d_factor = {21.5};
    result.satnum_id = {22, 23};
    result.ecl_index = {24};
    result.water_throughput = {25.0, 26.0};
    result.skin_pressure = {27.0, 28.0};
    result.water_velocity = {29.0, 30.0};
    result.filtrate_data = ConnFiltrateData::serializationTestObject();

    return result;
}

std::size_t PerfData::size() const {
    return this->pressure.size();
}

bool PerfData::empty() const {
    return this->pressure.empty();
}

bool PerfData::try_assign(const PerfData& other) {
    if (this->size() != other.size())
        return false;

    if (this->injector != other.injector)
        return false;

    this->pressure_first_connection = other.pressure_first_connection;
    this->pressure = other.pressure;
    this->rates = other.rates;
    this->phase_rates = other.phase_rates;
    this->phase_mixing_rates = other.phase_mixing_rates;
    this->solvent_rates = other.solvent_rates;
    this->polymer_rates = other.polymer_rates;
    this->brine_rates = other.brine_rates;
    this->water_throughput = other.water_throughput;
    this->skin_pressure = other.skin_pressure;
    this->water_velocity = other.water_velocity;
    this->prod_index = other.prod_index;
    this->micp_rates = other.micp_rates;
    this->filtrate_data = other.filtrate_data;
    return true;
}

bool PerfData::operator==(const PerfData& rhs) const
{
    return this->pressure_first_connection == rhs.pressure_first_connection &&
           this->pressure == rhs.pressure &&
           this->rates == rhs.rates &&
           this->phase_rates == rhs.phase_rates &&
           this->phase_mixing_rates == rhs.phase_mixing_rates &&
           this->solvent_rates == rhs.solvent_rates &&
           this->polymer_rates == rhs.polymer_rates &&
           this->brine_rates == rhs.brine_rates &&
           this->prod_index == rhs.prod_index &&
           this->micp_rates == rhs.micp_rates &&
           this->cell_index == rhs.cell_index &&
           this->connection_transmissibility_factor == rhs.connection_transmissibility_factor &&
           this->connection_d_factor == rhs.connection_d_factor &&
           this->satnum_id == rhs.satnum_id &&
           this->ecl_index == rhs.ecl_index &&
           this->water_throughput == rhs.water_throughput &&
           this->skin_pressure == rhs.skin_pressure &&
           this->water_velocity == rhs.water_velocity &&
           this->filtrate_data == rhs.filtrate_data;
}

}
