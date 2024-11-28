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

namespace Opm
{

template <class Scalar>
ConnectionData<Scalar>::
ConnectionData(std::size_t num_connection,
               std::size_t num_phases,
               std::size_t num_components)
  : pressure(num_connection)
  , surface_phase_rates(num_connection * num_phases)
  , reservoir_phase_rates(num_connection * num_phases)
  , total_molar_fractions(num_connection * num_components)
  , tranmissibility_factor(num_connection)
  , satnum_id(num_connection)
  , ecl_index(num_connection)
{
}


template <class Scalar>
ConnectionData<Scalar>::
ConnectionData(const std::vector<CompConnectionData<Scalar>>& connections,
               const PhaseUsage& phase_usage,
               const CompositionalConfig& comp_config)
  : ConnectionData(connections.size(), phase_usage.num_phases, comp_config.numComps())
{
    for (std::size_t con = 0; con < connections.size(); ++con) {
        this->tranmissibility_factor[con] = connections[con].connection_transmissibility_factor;
        this->satnum_id[con] = connections[con].satnum_id;
        this->ecl_index[con] = connections[con].cell_index;
    }
}

template <class Scalar>
SingleCompWellState<Scalar>::
SingleCompWellState(const std::string& well_name,
                    const CompositionalConfig& comp_config,
                    const PhaseUsage& phase_usage,
                    const std::vector<CompConnectionData<Scalar>>& connections,
                    bool is_producer)
   : name(well_name)
   , phase_usage(phase_usage)
   , producer(is_producer)
   , surface_phase_rates(phase_usage.num_phases)
   , phase_fractions(phase_usage.num_phases)
   , reservoir_phase_rates(phase_usage.num_phases)
   , total_molar_fractions(comp_config.numComps())
   , connection_data(connections, phase_usage, comp_config)
{

}

template <typename Scalar>
void SingleCompWellState<Scalar>::
update_injector_targets(const Well& well,
                        const SummaryState& st) {
    const auto& inj_controls = well.injectionControls(st);
    const auto cmode_is_bhp = (inj_controls.cmode == Well::InjectorCMode::BHP);
    assert(cmode_is_bhp && "Only BHP control mode is supported for now");
    const auto& injection_properties = well.getInjectionProperties();
    {
        const auto inection_type = injection_properties.injectorType;
        const bool is_gas_injecting = (inection_type == InjectorType::GAS);
        assert(is_gas_injecting && "Only gas injection is supported for now");
    }
    this->bhp = inj_controls.bhp_limit;
    this->injection_cmode = Well::InjectorCMode::BHP;
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
            assert(phase_usage.phase_used[BlackoilPhases::Aqua]);
            this->surface_phase_rates[phase_usage.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
            break;
        case InjectorType::GAS:
            assert(phase_usage.phase_used[BlackoilPhases::Vapour]);
            this->surface_phase_rates[phase_usage.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
            break;
        case InjectorType::OIL:
            assert(phase_usage.phase_used[BlackoilPhases::Liquid]);
            this->surface_phase_rates[phase_usage.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
            break;
        case InjectorType::MULTI:
            // Not currently handled, keep zero init.
            break;
    }
}

template <typename Scalar>
void SingleCompWellState<Scalar>::
update_producer_targets(const Well& well,
                        const SummaryState& st)
{
    const auto& prod_controls = well.productionControls(st);

    const auto cmode_is_undefined = (prod_controls.cmode == Well::ProducerCMode::CMODE_UNDEFINED);
    const auto cmode_is_bhp = (prod_controls.cmode == Well::ProducerCMode::BHP);
    assert((cmode_is_undefined || cmode_is_bhp) && "Only BHP control mode is supported for now");

    this->bhp = prod_controls.bhp_limit;
    this->production_cmode = Well::ProducerCMode::BHP;

    // we give a random rates for BHP controlled wells
    const Scalar production_rate = -1000.0 * Opm::unit::cubic(Opm::unit::meter) / Opm::unit::day;
    // different from the black oil situation, surface rate is not the primary variable, so we need to
    // check the how to convert it.
    // TODO: we should use our own phase index system
    // TODO: the following must be changed
    if (prod_controls.cmode == Well::ProducerCMode::BHP) {
        if (this->phase_usage.phase_used[BlackoilPhases::Liquid]) {
            this->surface_phase_rates[BlackoilPhases::Liquid] = production_rate;
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Aqua]) {
            this->surface_phase_rates[BlackoilPhases::Aqua] = production_rate;
        }
        if (this->phase_usage.phase_used[BlackoilPhases::Vapour]) {
            this->surface_phase_rates[BlackoilPhases::Vapour] = 100. * production_rate;
        }
    }


}

template <typename Scalar>
Scalar
SingleCompWellState<Scalar>::
get_total_surface_rate() const
{
    return std::accumulate(surface_phase_rates.begin(), surface_phase_rates.end(), Scalar(0));
}


} // namespace Opm