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

#ifndef OPM_SINGLE_COMP_WELL_STATE_HPP
#define OPM_SINGLE_COMP_WELL_STATE_HPP

#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/simulators/wells/PerforationData.hpp>

#include <string>
#include <vector>


namespace Opm {

template <typename Scalar>
class CompConnectionData
{
public:
    CompConnectionData() = default;
    CompConnectionData(std::size_t num_connection,
                       std::size_t num_phases,
                       std::size_t num_components);

    CompConnectionData(const std::vector<PerforationData<Scalar>>& connections,
                       const std::size_t num_phases,
                       const CompositionalConfig& comp_config);

    std::vector<Scalar> pressure {};
    std::vector<Scalar> surface_phase_rates {}; // surface phase rates
    std::vector<Scalar> reservoir_phase_rates {}; // phase rates
    std::vector<Scalar> total_molar_fractions {};

    // connection_transmissibility_factor
    std::vector<Scalar> transmissibility_factor {};
    std::vector<int> satnum_id {};
    std::vector<std::size_t> ecl_index {};
};

template <typename FluidSystem, typename Scalar>
class SingleCompWellState
{
public:
    SingleCompWellState(const std::string& name,
                        const CompositionalConfig& comp_config,
                        const Scalar temperature,
                        const std::vector<PerforationData<Scalar>>& connections,
                        bool is_producer);

    std::string name;
    bool producer;

    WellStatus status{WellStatus::OPEN};
    Scalar bhp{0};
    Scalar temperature{0};

    std::vector<Scalar> surface_phase_rates;
    std::vector<Scalar> phase_fractions; // V or L
    std::vector<Scalar> reservoir_phase_rates;
    // WZMF
    std::vector<Scalar> total_molar_fractions;
    // WXMF WYMF and WAMF
    std::vector<std::vector<Scalar> > phase_molar_fractions;

    CompConnectionData<Scalar> connection_data;

    WellInjectorCMode injection_cmode{WellInjectorCMode::CMODE_UNDEFINED};
    WellProducerCMode production_cmode{WellProducerCMode::CMODE_UNDEFINED};

    // TODO: the function can be reorganized so that we do not need to have initSingleInjector
    // and initSingleProducer, but we split when update the targets
    // so we have a function update_targets() to split between injector and producer
    void update_producer_targets(const Well& well,
                                 const std::vector<std::vector<Scalar>>& cell_mole_fractions,
                                 const SummaryState& st);

    void update_injector_targets(const Well& well, const SummaryState& st);

    Scalar get_total_surface_rate() const;
};

} // namespace Opm

#include "SingleCompWellState_impl.hpp"

#endif // OPM_SINGLE_COMP_WELL_STATE_HPP
