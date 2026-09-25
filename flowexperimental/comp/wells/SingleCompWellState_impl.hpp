/*
  Copyright 2024, SINTEF Digital

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

namespace Opm {

template <class Scalar>
CompConnectionData<Scalar>::
CompConnectionData(std::size_t num_connection,
                   std::size_t num_phases,
                   std::size_t num_components)
  : pressure(num_connection)
  , surface_phase_rates(num_connection * num_phases)
  , reservoir_phase_rates(num_connection * num_phases)
  , total_molar_fractions(num_connection * num_components)
  , transmissibility_factor(num_connection)
  , satnum_id(num_connection)
  , ecl_index(num_connection)
{
}


template <class Scalar>
CompConnectionData<Scalar>::
CompConnectionData(const std::vector<PerforationData<Scalar>>& connections,
                   const std::size_t num_phases,
                   const CompositionalConfig& comp_config)
  : CompConnectionData(connections.size(), num_phases, comp_config.numComps())
{
    for (std::size_t con = 0; con < connections.size(); ++con) {
        this->transmissibility_factor[con] = connections[con].connection_transmissibility_factor;
        this->satnum_id[con] = connections[con].satnum_id;
        this->ecl_index[con] = connections[con].cell_index;
    }
}

template <typename FluidSystem>
SingleCompWellState<FluidSystem>::
SingleCompWellState(const std::string& well_name,
                    const CompositionalConfig& comp_config,
                    const Scalar temperature_arg,
                    const std::vector<PerforationData<Scalar>>& connections,
                    bool is_producer)
   : name(well_name)
   , producer(is_producer)
   , temperature(temperature_arg)
   , surface_phase_rates(FluidSystem::numPhases)
   , phase_fractions(FluidSystem::numPhases)
   , reservoir_phase_rates(FluidSystem::numPhases)
   , total_molar_fractions(comp_config.numComps())
   , connection_data(connections, FluidSystem::numPhases, comp_config)
{
}

template <typename FluidSystem>
void SingleCompWellState<FluidSystem>::
update_injector_targets(const Well& well,
                        const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                        const SummaryState& st)
{
    const auto& inj_controls = well.injectionControls(st);
    const auto& injection_properties = well.getInjectionProperties();

    // Report this the way the black-oil model does in WellAssemble.
    if (inj_controls.cmode == Well::InjectorCMode::CMODE_UNDEFINED) {
        OPM_THROW(std::runtime_error,
                  "Well control must be specified for well " + this->name);
    }

    this->bhp = inj_controls.bhp_limit;
    this->injection_cmode = inj_controls.cmode;
    const auto injector_type = inj_controls.injector_type;
    if (injector_type == InjectorType::WATER && !FluidSystem::waterEnabled) {
        OPM_THROW(std::runtime_error,
                  "The water injector " + this->name + " needs an active water phase");
    }
    if (injector_type == InjectorType::WATER) {
        // The wellbore holds water alone. It still needs a hydrocarbon
        // composition the flash accepts, and a water injector has no stream
        // to take one from. A well without open local connections has no
        // well equations, so it can wait for a connection to open.
        if (!this->connection_data.ecl_index.empty()) {
            this->total_molar_fractions =
                cell_mole_fractions[this->connection_data.ecl_index.front()];
        }
        this->wellbore_water_volume_fraction = 1.;
    } else if (injector_type == InjectorType::GAS) {
        const auto& inj_composition = injection_properties.gasInjComposition();
        assert(this->total_molar_fractions.size() == inj_composition.size());
        // TODO: this might not be correct when crossing flow is involved
        this->total_molar_fractions = inj_composition;
    } else {
        OPM_THROW(std::runtime_error,
                  "Only gas and water injection is supported for well " + this->name);
    }

    // we initialize all open wells with a rate to avoid singularities
    Scalar inj_surf_rate = 10.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
        inj_surf_rate = inj_controls.surface_rate;
    }

    const auto injected_phase = injector_type == InjectorType::WATER ? FluidSystem::waterPhaseIdx
                                                                     : FluidSystem::gasPhaseIdx;
    this->surface_phase_rates[injected_phase] = inj_surf_rate;
}

template <typename FluidSystem>
void SingleCompWellState<FluidSystem>::
update_producer_targets(const Well& well,
                        const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                        const SummaryState& st)
{
    const auto& prod_controls = well.productionControls(st);

    // Report this the way the black-oil model does in WellAssemble.
    if (prod_controls.cmode == Well::ProducerCMode::CMODE_UNDEFINED) {
        OPM_THROW(std::runtime_error,
                  "Well control must be specified for well " + this->name);
    }

    // Read the reservoir composition only when an open perforation is
    // available. Initialize schedule targets below even without perforations.
    if (!this->connection_data.ecl_index.empty()) {
        this->total_molar_fractions = cell_mole_fractions[this->connection_data.ecl_index[0]];
    }

    this->bhp = prod_controls.bhp_limit;
    this->production_cmode = prod_controls.cmode;

    // Start BHP-controlled wells with a set of rates and rate-controlled wells
    // at their target. With a zero total rate, a rate control equation does
    // not depend on any primary variable and the well matrix is singular.
    const Scalar production_rate = -1000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    switch (prod_controls.cmode) {
    case Well::ProducerCMode::BHP:
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::oilPhaseIdx] = production_rate;
        }
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::waterPhaseIdx] = production_rate;
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::gasPhaseIdx] = 100. * production_rate;
        }
        break;
    case Well::ProducerCMode::ORAT:
        this->surface_phase_rates[FluidSystem::oilPhaseIdx] = -prod_controls.oil_rate;
        break;
    case Well::ProducerCMode::WRAT:
        if constexpr (FluidSystem::waterEnabled) {
            this->surface_phase_rates[FluidSystem::waterPhaseIdx] = -prod_controls.water_rate;
        }
        break;
    case Well::ProducerCMode::GRAT:
        this->surface_phase_rates[FluidSystem::gasPhaseIdx] = -prod_controls.gas_rate;
        break;
    case Well::ProducerCMode::LRAT:
        this->surface_phase_rates[FluidSystem::oilPhaseIdx] = -prod_controls.liquid_rate;
        break;
    default:
        break;
    }
}

template <typename FluidSystem>
void SingleCompWellState<FluidSystem>::
copyRuntimeStateFrom(const SingleCompWellState& other)
{
    // Keep the freshly initialized schedule-derived status, controls and
    // targets from base_init() and carry over the wellbore inventory, which a
    // shut well, or one that switched between producer and injector, lacks.
    if (status == WellStatus::SHUT || other.status == WellStatus::SHUT
        || producer != other.producer) {
        return;
    }
    bhp = other.bhp;
    wellbore_water_volume_fraction = other.wellbore_water_volume_fraction;
    total_molar_fractions = other.total_molar_fractions;
}

template <typename FluidSystem>
typename SingleCompWellState<FluidSystem>::Scalar
SingleCompWellState<FluidSystem>::
get_total_surface_rate() const
{
    return std::accumulate(surface_phase_rates.begin(), surface_phase_rates.end(), Scalar(0));
}


} // namespace Opm
