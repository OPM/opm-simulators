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

template <typename FluidSystem, class Scalar>
SingleCompWellState<FluidSystem, Scalar>::
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

template <typename FluidSystem, class Scalar>
void SingleCompWellState<FluidSystem, Scalar>::
update_injector_targets(const Well& well,
                        const SummaryState& st)
{
    const auto& inj_controls = well.injectionControls(st);
    const bool cmode_is_undefined = (inj_controls.cmode == Well::InjectorCMode::CMODE_UNDEFINED);
    assert(!cmode_is_undefined && "control types should be specified");
    const auto& injection_properties = well.getInjectionProperties();
    {
        const auto injection_type = injection_properties.injectorType;
        const bool is_gas_injecting = (injection_type == InjectorType::GAS);
        assert(is_gas_injecting && "Only gas injection is supported for now");
    }
    this->bhp = inj_controls.bhp_limit;
    this->injection_cmode = inj_controls.cmode;
    const auto& inj_composition = injection_properties.gasInjComposition();
    assert(this->total_molar_fractions.size() == inj_composition.size());
    // TODO: this might not be correct when crossing flow is involved
    this->total_molar_fractions = inj_composition;

    // we initialize all open wells with a rate to avoid singularities
    Scalar inj_surf_rate = 10.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
        inj_surf_rate = inj_controls.surface_rate;
    }

    switch (inj_controls.injector_type) {
        case InjectorType::WATER:
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            this->surface_phase_rates[FluidSystem::waterPhaseIdx] = inj_surf_rate;
            break;
        case InjectorType::GAS:
            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
            this->surface_phase_rates[FluidSystem::gasPhaseIdx] = inj_surf_rate;
            break;
        case InjectorType::OIL:
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            this->surface_phase_rates[FluidSystem::oilPhaseIdx] = inj_surf_rate;
            break;
        case InjectorType::MULTI:
            // Not currently handled, keep zero init.
            break;
    }
}

template <typename FluidSystem, class Scalar>
void SingleCompWellState<FluidSystem, Scalar>::
update_producer_targets(const Well& well,
                        const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                        const SummaryState& st)
{
    const auto& prod_controls = well.productionControls(st);

    const auto cmode_is_undefined = (prod_controls.cmode == Well::ProducerCMode::CMODE_UNDEFINED);
    assert(!cmode_is_undefined && "control types should be specified");

    this->total_molar_fractions = cell_mole_fractions[this->connection_data.ecl_index[0]];

    this->bhp = prod_controls.bhp_limit;
    this->production_cmode = prod_controls.cmode;

    // we give a set of rates for BHP-controlled wells for initialization
    const Scalar production_rate = -1000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    if (prod_controls.cmode == Well::ProducerCMode::BHP) {
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::oilPhaseIdx] = production_rate;
        }
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::waterPhaseIdx] = production_rate;
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->surface_phase_rates[FluidSystem::gasPhaseIdx] = 100. * production_rate;
        }
    }
}

template <typename FluidSystem, class Scalar>
Scalar
SingleCompWellState<FluidSystem, Scalar>::
get_total_surface_rate() const
{
    return std::accumulate(surface_phase_rates.begin(), surface_phase_rates.end(), Scalar(0));
}


} // namespace Opm
