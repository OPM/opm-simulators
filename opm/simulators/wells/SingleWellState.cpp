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

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/wells/PerforationData.hpp>

#include <fmt/format.h>

namespace Opm {

template<typename Scalar, typename IndexTraits>
SingleWellState<Scalar, IndexTraits>::
SingleWellState(const std::string& name_,
                const ParallelWellInfo<Scalar>& pinfo,
                const PhaseUsageInfo<IndexTraits>&  pu_arg,
                bool is_producer,
                Scalar pressure_first_connection_,
                const std::vector<PerforationData<Scalar>>& perf_input,
                Scalar temp)
    : name(name_)
    , parallel_info(pinfo)
    , producer(is_producer)
    , pu(pu_arg)
    , pressure_first_connection(pressure_first_connection_)
    , temperature(temp)
    , well_potentials(pu.numActivePhases())
    , productivity_index(pu.numActivePhases())
    , implicit_ipr_a(pu.numActivePhases())
    , implicit_ipr_b(pu.numActivePhases())
    , surface_rates(pu.numActivePhases())
    , reservoir_rates(pu.numActivePhases())
    , prev_surface_rates(pu.numActivePhases())
    , perf_data(perf_input.size(), !is_producer, pu.numActivePhases())
    , trivial_group_target(false)
{
    for (std::size_t perf = 0; perf < perf_input.size(); perf++) {
        this->perf_data.cell_index[perf] = perf_input[perf].cell_index;
        this->perf_data.connection_transmissibility_factor[perf] = perf_input[perf].connection_transmissibility_factor;
        this->perf_data.connection_d_factor[perf] = perf_input[perf].connection_d_factor;
        this->perf_data.satnum_id[perf] = perf_input[perf].satnum_id;
        this->perf_data.ecl_index[perf] = perf_input[perf].ecl_index;
    }
}

template<typename Scalar, typename IndexTraits>
SingleWellState<Scalar, IndexTraits> SingleWellState<Scalar, IndexTraits>::
serializationTestObject(const ParallelWellInfo<Scalar>& pinfo)
{
    SingleWellState<Scalar, IndexTraits> result("testing", pinfo, PhaseUsageInfo<IndexTraits>{},
                                                true, 1.0, {}, 2.0);
    result.perf_data = PerfData<Scalar>::serializationTestObject();

    return result;
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::init_timestep(const SingleWellState& other)
{
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

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::shut()
{
    this->bhp = 0;
    this->thp = 0;
    this->status = Well::Status::SHUT;
    std::fill(this->surface_rates.begin(), this->surface_rates.end(), 0);
    std::fill(this->prev_surface_rates.begin(), this->prev_surface_rates.end(), 0);
    std::fill(this->reservoir_rates.begin(), this->reservoir_rates.end(), 0);
    std::fill(this->productivity_index.begin(), this->productivity_index.end(), 0);
    std::fill(this->implicit_ipr_a.begin(), this->implicit_ipr_a.end(), 0);
    std::fill(this->implicit_ipr_b.begin(), this->implicit_ipr_b.end(), 0);

    auto& connpi = this->perf_data.prod_index;
    connpi.assign(connpi.size(), 0);
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::stop()
{
    this->thp = 0;
    this->status = Well::Status::STOP;
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::open()
{
    this->status = Well::Status::OPEN;
}

template<typename Scalar, typename IndexTraits>
bool SingleWellState<Scalar, IndexTraits>::updateStatus(Well::Status new_status)
{
    if (this->status == new_status) {
        return false;
    }

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
    return true;
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::
reset_connection_factors(const std::vector<PerforationData<Scalar>>& new_perf_data)
{
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

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::
sum_connection_rates(const std::vector<Scalar>& connection_rates) const
{
    return this->parallel_info.get().sumPerfValues(connection_rates.begin(), connection_rates.end());
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_brine_rates() const
{
    return this->sum_connection_rates(this->perf_data.brine_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_polymer_rates() const
{
    return this->sum_connection_rates(this->perf_data.polymer_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_solvent_rates() const
{
    return this->sum_connection_rates(this->perf_data.solvent_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_microbial_rates() const
{
    return this->sum_connection_rates(this->perf_data.microbial_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_oxygen_rates() const
{
    return this->sum_connection_rates(this->perf_data.oxygen_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_urea_rates() const
{
    return this->sum_connection_rates(this->perf_data.urea_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_wat_mass_rates() const
{
    return this->sum_connection_rates(this->perf_data.wat_mass_rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_filtrate_rate() const
{
    if (this->producer) return 0.;

    return this->sum_connection_rates(this->perf_data.filtrate_data.rates);
}

template<typename Scalar, typename IndexTraits>
Scalar SingleWellState<Scalar, IndexTraits>::sum_filtrate_total() const
{
    if (this->producer) return 0.;

    return this->sum_connection_rates(this->perf_data.filtrate_data.total);
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::
update_producer_targets(const Well& ecl_well, const SummaryState& st)
{
    constexpr Scalar bhp_safety_factor = 0.99;
    const auto& prod_controls = ecl_well.productionControls(st);

    auto cmode_is_bhp = (prod_controls.cmode == Well::ProducerCMode::BHP);
    auto bhp_limit = prod_controls.bhp_limit;

    if (ecl_well.getStatus() == Well::Status::STOP) {
        if (cmode_is_bhp)
            this->bhp = bhp_limit;
        else
            this->bhp = this->pressure_first_connection;

        return;
    }

    std::fill(this->surface_rates.begin(), this->surface_rates.end(), 0.0);
    switch (prod_controls.cmode) {
    case Well::ProducerCMode::ORAT:
        assert(pu.phaseIsActive(oilPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(oilPhaseIdx)] = -prod_controls.oil_rate;
        break;
    case Well::ProducerCMode::WRAT:
        assert(pu.phaseIsActive(waterPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(waterPhaseIdx)] = -prod_controls.water_rate;
        break;
    case Well::ProducerCMode::GRAT:
        assert(pu.phaseIsActive(gasPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(gasPhaseIdx)] = -prod_controls.gas_rate;
        break;
    case Well::ProducerCMode::GRUP:
    case Well::ProducerCMode::THP:
    case Well::ProducerCMode::BHP:
        // Keeping all rates at zero, they will be initialized properly in
        // a call to WellInterface::initializeProducerWellState() later, which will
        // use the reservoir state to find a better initial value.
        // This also applies to the ORAT/WRAT/GRAT above, but then only the
        // rate not set in the above will be modified in
        // WellInterface::initializeProducerWellState().
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
        this->bhp = this->pressure_first_connection * bhp_safety_factor;
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::
update_injector_targets(const Well& ecl_well, const SummaryState& st)
{
    const Scalar bhp_safety_factor = 1.01;
    const auto& inj_controls = ecl_well.injectionControls(st);

    if (inj_controls.hasControl(Well::InjectorCMode::THP))
        this->thp = inj_controls.thp_limit;

    auto cmode_is_bhp = (inj_controls.cmode == Well::InjectorCMode::BHP);
    auto bhp_limit = inj_controls.bhp_limit;

    if (ecl_well.getStatus() == Well::Status::STOP) {
        if (cmode_is_bhp)
            this->bhp = bhp_limit;
        else
            this->bhp = this->pressure_first_connection;

        return;
    }

    // we initialize all open wells with a rate to avoid singularities
    Scalar inj_surf_rate = 10.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
        inj_surf_rate = inj_controls.surface_rate;
    }

    switch (inj_controls.injector_type) {
    case InjectorType::WATER:
        assert(pu.phaseIsActive(waterPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(waterPhaseIdx)] = inj_surf_rate;
        break;
    case InjectorType::GAS:
        assert(pu.phaseIsActive(gasPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(gasPhaseIdx)] = inj_surf_rate;
        break;
    case InjectorType::OIL:
        assert(pu.phaseIsActive(oilPhaseIdx));
        this->surface_rates[pu.canonicalToActivePhaseIdx(oilPhaseIdx)] = inj_surf_rate;
        break;
    case InjectorType::MULTI:
        // Not currently handled, keep zero init.
        break;
    }

    if (cmode_is_bhp)
        this->bhp = bhp_limit;
    else
        this->bhp = this->pressure_first_connection * bhp_safety_factor;
}

template<typename Scalar, typename IndexTraits>
void SingleWellState<Scalar, IndexTraits>::
update_type_and_targets(const Well& ecl_well, const SummaryState& st)
{
    if (this->producer != ecl_well.isProducer()) {
        // type has changed due to ACTIONX
        // Make sure that we are consistent with the ecl_well
        const bool switchedToProducer = this->producer = ecl_well.isProducer();
        if (switchedToProducer) {
            this->production_cmode = ecl_well.productionControls(st).cmode;
            // clear injection rates (those are positive)
            std::transform(this->surface_rates.begin(), this->surface_rates.end(),
                           this->surface_rates.begin(),
                           [](const Scalar& val){ return std::min(Scalar(), val);});
        } else {
            perf_data.prepareInjectorContainers();
            this->injection_cmode = ecl_well.injectionControls(st).cmode;
            // clear production rates (those are negative)
            std::transform(this->surface_rates.begin(), this->surface_rates.end(),
                           this->surface_rates.begin(),
                           [](const Scalar& val){ return std::max(Scalar(), val);});
        }
    }

    if (this->producer)
        this->update_producer_targets(ecl_well, st);
    else
        this->update_injector_targets(ecl_well, st);

}

template<typename Scalar, typename IndexTraits>
bool SingleWellState<Scalar, IndexTraits>::operator==(const SingleWellState& rhs) const
{
    // TODO: we are missing pu, while it was not there in the first place
    return this->name == rhs.name &&
           this->status == rhs.status &&
           this->producer == rhs.producer &&
           this->bhp == rhs.bhp &&
           this->thp == rhs.thp &&
           this->pressure_first_connection == rhs.pressure_first_connection &&
           this->temperature == rhs.temperature &&
           this->phase_mixing_rates == rhs.phase_mixing_rates &&
           this->well_potentials == rhs.well_potentials &&
           this->productivity_index == rhs.productivity_index &&
           this->implicit_ipr_a == rhs.implicit_ipr_a &&
           this->implicit_ipr_b == rhs.implicit_ipr_b &&
           this->surface_rates == rhs.surface_rates &&
           this->reservoir_rates == rhs.reservoir_rates &&
           this->prev_surface_rates == rhs.prev_surface_rates &&
           this->perf_data == rhs.perf_data &&
           this->filtrate_conc == rhs.filtrate_conc &&
           this->trivial_group_target == rhs.trivial_group_target &&
           this->segments == rhs.segments &&
           this->events == rhs.events &&
           this->injection_cmode == rhs.injection_cmode &&
           this->production_cmode == rhs.production_cmode &&
           this->alq_state == rhs.alq_state &&
           this->primaryvar == rhs.primaryvar &&
           this->group_target == rhs.group_target;
}

template<typename Scalar, typename IndexTraits>
std::string
SingleWellState<Scalar, IndexTraits>::
debugInfo() const
{
    std::string info = "Well name: " + this->name + " well state:\n";
    fmt::format_to(std::back_inserter(info), " type: {}, staus: {}, control type: {}, BHP: {:8.2e} Pa, THP: {:8.2e} Pa\n",
                   this->producer? "Producer" : " Injector", WellStatus2String(this->status),
                   this->producer ? WellProducerCMode2String(this->production_cmode)
                                  : WellInjectorCMode2String(this->injection_cmode),
                   this->bhp / Opm::unit::barsa,
                   this->thp / Opm::unit::barsa);

    fmt::format_to(std::back_inserter(info), "  Surface rates (m3/s): ");
    for (const auto& r : this->surface_rates) {
        fmt::format_to(std::back_inserter(info), " {:8.2e}", r);
    }
    fmt::format_to(std::back_inserter(info), "\n");
    fmt::format_to(std::back_inserter(info), "  Connection pressures and connection rates:\n");
    for (std::size_t perf = 0; perf < this->perf_data.size(); perf++) {
        fmt::format_to(std::back_inserter(info),
                       "  Connection {:4}: Pressure: {:8.2e} Pa, Rate:",
                       perf,
                       this->perf_data.pressure[perf]);
        const auto& connection_rates = this->perf_data.phase_rates;
        const std::size_t num_phases = this->pu.numActivePhases();

        for (std::size_t p = 0; p < num_phases; ++p) {
            fmt::format_to(std::back_inserter(info), " {: 8.2e} m3/s", connection_rates[perf * num_phases + p]);
        }
        fmt::format_to(std::back_inserter(info), "\n");
    }

    info += this->segments.debugInfo();

    return info;
}

template class SingleWellState<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class SingleWellState<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
