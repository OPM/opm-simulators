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
#include <opm/simulators/wells/PerfData.hpp>

#include <opm/simulators/wells/ConnFiltrateData.hpp>

namespace Opm {

template<class Scalar>
PerfData<Scalar>::PerfData(const std::size_t num_perf,
                           const Scalar pressure_first_connection_,
                           const bool injector_,
                           const std::size_t num_phases)
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
    , microbial_rates(num_perf)
    , oxygen_rates(num_perf)
    , urea_rates(num_perf)
    , cell_index(num_perf)
    , connection_transmissibility_factor(num_perf)
    , connection_d_factor(num_perf)
    , connection_compaction_tmult(num_perf)
    , satnum_id(num_perf)
    , ecl_index(num_perf)
    , gas_mass_rates(num_perf)
{
    if (injector) {
        prepareInjectorContainers();
    }
}

template<class Scalar>
void PerfData<Scalar>::prepareInjectorContainers()
{
    auto num_perf = pressure.size();
    this->water_throughput.resize(num_perf);
    this->skin_pressure.resize(num_perf);
    this->water_velocity.resize(num_perf);
    this->filtrate_data.resize(num_perf);
}

template<class Scalar>
PerfData<Scalar> PerfData<Scalar>::serializationTestObject()
{
    PerfData result;
    result.injector = true;
    result.pressure_first_connection = 1.0;
    result.pressure = {2.0, 3.0, 4.0};
    result.rates = {5.0, 6.0};
    result.phase_rates = {7.0};
    result.phase_mixing_rates = { {1.0, 2.0, 3.0, 4.0}};
    result.solvent_rates = {8.0, 9.0};
    result.polymer_rates = {10.0, 11.0, 12.0};
    result.brine_rates = {13.0};
    result.prod_index = {14.0, 15.0};
    result.microbial_rates = {16.0};
    result.oxygen_rates = {16.0};
    result.urea_rates = {16.0};
    result.cell_index = {17, 18, 19, 20};
    result.connection_transmissibility_factor = {21.0};
    result.connection_d_factor = {21.5};
    result.connection_compaction_tmult.assign(1, 21.75);
    result.satnum_id = {22, 23};
    result.ecl_index = {24};
    result.water_throughput = {25.0, 26.0};
    result.skin_pressure = {27.0, 28.0};
    result.water_velocity = {29.0, 30.0};
    result.filtrate_data = ConnFiltrateData<Scalar>::serializationTestObject();
    result.connFracStatistics.assign(3, ConnFracStatistics<Scalar>::serializationTestObject());
    result.gas_mass_rates = {31.0};

    return result;
}

template<class Scalar>
std::size_t PerfData<Scalar>::size() const
{
    return this->pressure.size();
}

template<class Scalar>
bool PerfData<Scalar>::empty() const
{
    return this->pressure.empty();
}

template<class Scalar>
bool PerfData<Scalar>::try_assign(const PerfData& other)
{
    if (this->size() != other.size()) {
        return false;
    }

    if (this->injector != other.injector) {
        return false;
    }

    this->pressure_first_connection = other.pressure_first_connection;
    this->pressure = other.pressure;
    this->rates = other.rates;
    this->phase_rates = other.phase_rates;
    this->phase_mixing_rates = other.phase_mixing_rates;
    this->solvent_rates = other.solvent_rates;
    this->polymer_rates = other.polymer_rates;
    this->brine_rates = other.brine_rates;
    this->prod_index = other.prod_index;
    this->microbial_rates = other.microbial_rates;
    this->oxygen_rates = other.oxygen_rates;
    this->urea_rates = other.urea_rates;
    this->water_throughput = other.water_throughput;
    this->skin_pressure = other.skin_pressure;
    this->water_velocity = other.water_velocity;
    this->filtrate_data = other.filtrate_data;
    this->connFracStatistics = other.connFracStatistics;
    this->gas_mass_rates = other.gas_mass_rates;

    return true;
}

template<class Scalar>
bool PerfData<Scalar>::operator==(const PerfData& rhs) const
{
    return (this->injector == rhs.injector)
        && (this->pressure_first_connection == rhs.pressure_first_connection)
        && (this->pressure == rhs.pressure)
        && (this->rates == rhs.rates)
        && (this->phase_rates == rhs.phase_rates)
        && (this->phase_mixing_rates == rhs.phase_mixing_rates)
        && (this->solvent_rates == rhs.solvent_rates)
        && (this->polymer_rates == rhs.polymer_rates)
        && (this->brine_rates == rhs.brine_rates)
        && (this->prod_index == rhs.prod_index)
        && (this->microbial_rates == rhs.microbial_rates)
        && (this->oxygen_rates == rhs.oxygen_rates)
        && (this->urea_rates == rhs.urea_rates)
        && (this->cell_index == rhs.cell_index)
        && (this->connection_transmissibility_factor == rhs.connection_transmissibility_factor)
        && (this->connection_d_factor == rhs.connection_d_factor)
        && (this->connection_compaction_tmult == rhs.connection_compaction_tmult)
        && (this->satnum_id == rhs.satnum_id)
        && (this->ecl_index == rhs.ecl_index)
        && (this->water_throughput == rhs.water_throughput)
        && (this->skin_pressure == rhs.skin_pressure)
        && (this->water_velocity == rhs.water_velocity)
        && (this->filtrate_data == rhs.filtrate_data)
        && (this->connFracStatistics == rhs.connFracStatistics)
        && (this->gas_mass_rates == rhs.gas_mass_rates)
        ;
}

template class PerfData<double>;

#if FLOW_INSTANTIATE_FLOAT
template class PerfData<float>;
#endif

} // namespace Opm
