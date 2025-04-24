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

#ifndef OPM_COMP_WELL_STATES_HPP
#define OPM_COMP_WELL_STATES_HPP

// TODO: evaluate whether we can continue use this phase usage
// if yes, we might want to rename it
// TODO: we also the informaton regarding the components
// #include <opm/simulators/utils/BlackoilPhases.hpp>
#include <opm/simulators/wells/WellContainer.hpp>
#include <opm/simulators/wells/PerforationData.hpp>

#include <opm/output/data/Wells.hpp>


#include <flowexperimental/comp/wells/SingleCompWellState.hpp>


namespace Opm {

template <typename FluidSystem, typename Scalar>
class CompWellState
{
public:
    using CompConnectionData = PerforationData<Scalar>;
    using SingleWellState = SingleCompWellState<FluidSystem, Scalar>;

    explicit CompWellState(const CompositionalConfig& comp_config);

    void init(const std::vector<Well>& wells_ecl,
              const std::vector<Scalar>& cell_pressures,
              const Scalar cell_temperature,
              const std::vector<std::vector<Scalar>>& cell_mole_fractions,
              const std::vector<std::vector<CompConnectionData> >& well_connection_data,
              const SummaryState& sumary_state,
              const CompWellState* prev_well_state = nullptr);

    const SingleWellState& operator[](const std::string& well_name) const;

    SingleWellState& operator[](const std::string& well_name);

    data::Wells report() const;

private:
    WellContainer<SingleWellState> wells_;

    const CompositionalConfig& comp_config_;

    void base_init(const std::vector<Well>& wells_ecl,
                   const std::vector<Scalar>& cell_pressures,
                   const Scalar temperature,
                   const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                   const std::vector<std::vector<CompConnectionData> >& well_connection_data,
                   const SummaryState& summary_state);

    void initSingleWell(const Well& well,
                        const std::vector<Scalar>& cell_pressures,
                        const Scalar tempearture,
                        const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                        const std::vector<CompConnectionData >& conn_data,
                        const SummaryState& summary_state);

    void initSingleInjector(const Well& well,
                            const std::vector<Scalar>& cell_pressures,
                            const Scalar temperature,
                            const std::vector<CompConnectionData>& conn_data,
                            const SummaryState& summary_state);

    void initSingleProducer(const Well& well,
                            const std::vector<Scalar>& cell_pressures,
                            const Scalar temperature,
                            const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                            const std::vector<CompConnectionData>& conn_data,
                            const SummaryState& summary_state);
};

} // end of namespace Opm

#include "CompWellState_impl.hpp"

#endif // OPM_COMP_WELL_STATES_HPP
