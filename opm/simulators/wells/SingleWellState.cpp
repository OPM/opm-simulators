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

#include <config.h>
#include <opm/simulators/wells/SingleWellState.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/wells/PerforationData.hpp>

namespace Opm {

SingleWellState::SingleWellState(const std::string& name_,
                                 const ParallelWellInfo& pinfo,
                                 bool is_producer,
                                 double pressure_first_connection,
                                 const std::vector<PerforationData>& perf_input,
                                 const PhaseUsage& pu_,
                                 double temp)
    : name(name_)
    , parallel_info(pinfo)
    , producer(is_producer)
    , pu(pu_)
    , temperature(temp)
    , well_potentials(pu_.num_phases)
    , productivity_index(pu_.num_phases)
    , surface_rates(pu_.num_phases)
    , reservoir_rates(pu_.num_phases)
    , prev_surface_rates(pu_.num_phases)
    , perf_data(perf_input.size(), pressure_first_connection, !is_producer, pu_.num_phases)
    , trivial_target(false)
{
    for (std::size_t perf = 0; perf < perf_input.size(); perf++) {
        this->perf_data.cell_index[perf] = perf_input[perf].cell_index;
        this->perf_data.connection_transmissibility_factor[perf] = perf_input[perf].connection_transmissibility_factor;
        this->perf_data.connection_d_factor[perf] = perf_input[perf].connection_d_factor;
        this->perf_data.satnum_id[perf] = perf_input[perf].satnum_id;
        this->perf_data.ecl_index[perf] = perf_input[perf].ecl_index;
    }
}

SingleWellState SingleWellState::serializationTestObject(const ParallelWellInfo& pinfo)
{
    SingleWellState result("testing", pinfo, true, 1.0, {}, PhaseUsage{}, 2.0);
    result.perf_data = PerfData::serializationTestObject();

    return result;
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
    std::fill(this->prev_surface_rates.begin(), this->prev_surface_rates.end(), 0);
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

void SingleWellState::reset_connection_factors(const std::vector<PerforationData>& new_perf_data) {
    if (this->perf_data.size() != new_perf_data.size()) {
        throw std::invalid_argument {
            "Size mismatch for perforation data in well " + this->name
        };
    }

    for (std::size_t conn_index = 0; conn_index < new_perf_data.size(); conn_index++) {
        if (this->perf_data.cell_index[conn_index] != static_cast<std::size_t>(new_perf_data[conn_index].cell_index)) {
            throw std::invalid_argument {
                "Cell index mismatch in connection "
                    + std::to_string(conn_index)
                    + " of well "
                    + this->name
            };
        }

        if (this->perf_data.satnum_id[conn_index] != new_perf_data[conn_index].satnum_id) {
            throw std::invalid_argument {
                "Saturation function table mismatch in connection "
                + std::to_string(conn_index)
                        + " of well "
                        + this->name
            };
        }
        this->perf_data.connection_transmissibility_factor[conn_index] = new_perf_data[conn_index].connection_transmissibility_factor;
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

double SingleWellState::sum_filtrate_rate() const {
    if (this->producer) return 0.;

    return this->sum_connection_rates(this->perf_data.filtrate_data.rates);
}

double SingleWellState::sum_filtrate_total() const {
    if (this->producer) return 0.;

    return this->sum_connection_rates(this->perf_data.filtrate_data.total);
}

void SingleWellState::update_producer_targets(const Well& ecl_well, const SummaryState& st) {
    const double bhp_safety_factor = 0.99;
    const auto& prod_controls = ecl_well.productionControls(st);

    auto cmode_is_bhp = (prod_controls.cmode == Well::ProducerCMode::BHP);
    auto bhp_limit = prod_controls.bhp_limit;

    if (ecl_well.getStatus() == Well::Status::STOP) {
        if (cmode_is_bhp)
            this->bhp = bhp_limit;
        else
            this->bhp = this->perf_data.pressure_first_connection;

        return;
    }

    switch (prod_controls.cmode) {
    case Well::ProducerCMode::ORAT:
        assert(this->pu.phase_used[BlackoilPhases::Liquid]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
        break;
    case Well::ProducerCMode::WRAT:
        assert(this->pu.phase_used[BlackoilPhases::Aqua]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
        break;
    case Well::ProducerCMode::GRAT:
        assert(this->pu.phase_used[BlackoilPhases::Vapour]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
        break;
    case Well::ProducerCMode::GRUP:
    case Well::ProducerCMode::THP:
    case Well::ProducerCMode::BHP:
        if (this->pu.phase_used[BlackoilPhases::Liquid]) {
            this->surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = -1000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
        }
        if (this->pu.phase_used[BlackoilPhases::Aqua]) {
            this->surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = -1000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
        }
        if (this->pu.phase_used[BlackoilPhases::Vapour]){
            this->surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = -100000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
        }
        break;

    default:
        // Keep zero init.
        break;
    }

    if (prod_controls.cmode == Well::ProducerCMode::THP) {
        this->thp = prod_controls.thp_limit;
    }

    if (cmode_is_bhp)
        this->bhp = bhp_limit;
    else
        this->bhp = this->perf_data.pressure_first_connection * bhp_safety_factor;

}

void SingleWellState::update_injector_targets(const Well& ecl_well, const SummaryState& st) {
    const double bhp_safety_factor = 1.01;
    const auto& inj_controls = ecl_well.injectionControls(st);

    if (inj_controls.hasControl(Well::InjectorCMode::THP))
        this->thp = inj_controls.thp_limit;

    auto cmode_is_bhp = (inj_controls.cmode == Well::InjectorCMode::BHP);
    auto bhp_limit = inj_controls.bhp_limit;

    if (ecl_well.getStatus() == Well::Status::STOP) {
        if (cmode_is_bhp)
            this->bhp = bhp_limit;
        else
            this->bhp = this->perf_data.pressure_first_connection;

        return;
    }

    // we initialize all open wells with a rate to avoid singularities
    double inj_surf_rate = 10.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
        inj_surf_rate = inj_controls.surface_rate;
    }

    switch (inj_controls.injector_type) {
    case InjectorType::WATER:
        assert(pu.phase_used[BlackoilPhases::Aqua]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
        break;
    case InjectorType::GAS:
        assert(pu.phase_used[BlackoilPhases::Vapour]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
        break;
    case InjectorType::OIL:
        assert(pu.phase_used[BlackoilPhases::Liquid]);
        this->surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
        break;
    case InjectorType::MULTI:
        // Not currently handled, keep zero init.
        break;
    }

    if (cmode_is_bhp)
        this->bhp = bhp_limit;
    else
        this->bhp = this->perf_data.pressure_first_connection * bhp_safety_factor;
}

void SingleWellState::update_targets(const Well& ecl_well, const SummaryState& st) {
    if (this->producer)
        this->update_producer_targets(ecl_well, st);
    else
        this->update_injector_targets(ecl_well, st);
}

bool SingleWellState::operator==(const SingleWellState& rhs) const
{
    return this->name == rhs.name &&
           this->status == rhs.status &&
           this->producer == rhs.producer &&
           this->bhp == rhs.bhp &&
           this->thp == rhs.thp &&
           this->temperature == rhs.temperature &&
           this->phase_mixing_rates == rhs.phase_mixing_rates &&
           this->well_potentials == rhs.well_potentials &&
           this->productivity_index == rhs.productivity_index &&
           this->surface_rates == rhs.surface_rates &&
           this->reservoir_rates == rhs.reservoir_rates &&
           this->prev_surface_rates == rhs.prev_surface_rates &&
           this->perf_data == rhs.perf_data &&
           this->filtrate_conc == rhs.filtrate_conc &&
           this->trivial_target == rhs.trivial_target &&
           this->segments == rhs.segments &&
           this->events == rhs.events &&
           this->injection_cmode == rhs.injection_cmode &&
           this->production_cmode == rhs.production_cmode;
}

}
