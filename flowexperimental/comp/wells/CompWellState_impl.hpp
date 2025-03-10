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

template <typename Scalar>
CompWellState<Scalar>::
CompWellState(const PhaseUsage& phase_usage,
              const CompositionalConfig& comp_config)
    : phase_usage_(phase_usage)
    , comp_config_(comp_config)
{
}


template <typename Scalar>
void CompWellState<Scalar>::
init(const std::vector<Well>& wells_ecl,
     const std::vector<Scalar>& cell_pressures,
     const std::vector<std::vector<Scalar>>& cell_mole_fractions,
     const std::vector<std::vector<CompConnectionData<Scalar> > >& well_connection_data,
     const SummaryState& summary_state,
     const CompWellState* /*prev_well_state*/)
{
    this->base_init(wells_ecl, cell_pressures, cell_mole_fractions, well_connection_data, summary_state);

    // TODO: for the simple case we have, I think we can just copy from the prev_well_state
    // let us see how we gonna use it though
}

template <typename Scalar>
void CompWellState<Scalar>::
base_init(const std::vector<Well>& wells_ecl,
          const std::vector<Scalar>& cell_pressures,
          const std::vector<std::vector<Scalar>>& cell_mole_fractions,
          const std::vector<std::vector<CompConnectionData<Scalar> > >& well_connection_data,
          const SummaryState& summary_state)
{
    this->wells_.clear();

    const auto num_wells = wells_ecl.size();

    for (auto w = 0*num_wells; w < num_wells; ++w) {
        const Well& well = wells_ecl[w];
        const auto& conn_data = well_connection_data[w];
        initSingleWell(well, cell_pressures, cell_mole_fractions, conn_data, summary_state);
    }

}

template <typename Scalar>
void CompWellState<Scalar>::
initSingleWell(const Well& well,
               const std::vector<Scalar>& cell_pressures,
               const std::vector<std::vector<Scalar>>& cell_mole_fractions,
               const std::vector<CompConnectionData<Scalar> >& conn_data,
               const SummaryState& summary_state)
{
    if (well.isInjector()) {
        initSingleInjector(well, cell_pressures, conn_data, summary_state);
    } else {
        initSingleProducer(well, cell_pressures, cell_mole_fractions, conn_data, summary_state);
    }

}

template <typename Scalar>
void CompWellState<Scalar>::
initSingleInjector(const Well& well,
                   const std::vector<Scalar>& /* cell_pressures */,
                   const std::vector<CompConnectionData<Scalar> >& conn_data,
                   const SummaryState& summary_state)
{
    auto& ws = this->wells_.add(well.name(),
                                SingleCompWellState<Scalar>(well.name(),
                                    this->comp_config_,
                                    this->phase_usage_,
                                    conn_data,
                                    false) );
    ws.update_injector_targets(well, summary_state);
}

template <typename Scalar>
void CompWellState<Scalar>::
initSingleProducer(const Well& well,
                   const std::vector<Scalar>& /* cell_pressures */,
                   const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                   const std::vector<CompConnectionData<Scalar> >& conn_data,
                   const SummaryState& summary_state)
{
    auto& ws = this->wells_.add(well.name(),
                                SingleCompWellState<Scalar>(well.name(),
                                    this->comp_config_,
                                    this->phase_usage_,
                                    conn_data,
                                    true) );
    ws.update_producer_targets(well, cell_mole_fractions, summary_state);
}

template <typename Scalar>
const SingleCompWellState<Scalar>&
CompWellState<Scalar>::
operator[](const std::string& well_name) const
{
    return this->wells_[well_name];
}

template <typename Scalar>
SingleCompWellState<Scalar>&
CompWellState<Scalar>::
operator[](const std::string& well_name)
{
    return this->wells_[well_name];
}

template <typename Scalar>
data::Wells
CompWellState<Scalar>::
report() const
{
    if (this->wells_.empty()) {
        return {};
    }
    using rt = data::Rates::opt;
    const auto& pu = this->phase_usage_;

    data::Wells res;
    for (std::size_t w = 0; w < this->wells_.size(); ++w) {
        const auto& ws = this->wells_[w];
        if (ws.status == Well::Status::SHUT) {
            continue;
        }
        auto& well = res[ws.name];
        well.bhp = ws.bhp;
        well.temperature = ws.temperature;
        const auto& surface_rates = ws.surface_phase_rates;
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            well.rates.set(rt::wat, surface_rates[pu.phase_pos[BlackoilPhases::Aqua]]);
        }
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            well.rates.set(rt::oil, surface_rates[pu.phase_pos[BlackoilPhases::Liquid]]);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            well.rates.set(rt::gas, surface_rates[pu.phase_pos[BlackoilPhases::Vapour]]);
        }
    }
    return res;
}

} // end of namespace Opm