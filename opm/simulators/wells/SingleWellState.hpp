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

#ifndef OPM_SINGLE_WELL_STATE_HEADER_INCLUDED
#define OPM_SINGLE_WELL_STATE_HEADER_INCLUDED

#include <functional>
#include <vector>

#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>

#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/ALQState.hpp>

namespace Opm {

template<class Scalar> struct PerforationData;
class SummaryState;
class Well;

template<typename Scalar, typename IndexTraits>
class SingleWellState {
public:
    static const int waterPhaseIdx = PhaseUsageInfo<IndexTraits>::waterPhaseIdx;
    static const int oilPhaseIdx = PhaseUsageInfo<IndexTraits>::oilPhaseIdx;
    static const int gasPhaseIdx = PhaseUsageInfo<IndexTraits>::gasPhaseIdx;

    SingleWellState(const std::string& name,
                    const ParallelWellInfo<Scalar>& pinfo,
                    const PhaseUsageInfo<IndexTraits>& pu,
                    bool is_producer,
                    Scalar pressure_first_connection,
                    const std::vector<PerforationData<Scalar>>& perf_input,
                    Scalar temp);

    static SingleWellState serializationTestObject(const ParallelWellInfo<Scalar>& pinfo);

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(name);
        serializer(status);
        serializer(producer);
        serializer(bhp);
        serializer(thp);
        serializer(pressure_first_connection);
        serializer(temperature);
        serializer(energy_rate);
        serializer(efficiency_scaling_factor);
        serializer(phase_mixing_rates);
        serializer(well_potentials);
        serializer(productivity_index);
        serializer(implicit_ipr_a);
        serializer(implicit_ipr_b);
        serializer(surface_rates);
        serializer(reservoir_rates);
        serializer(prev_surface_rates);
        serializer(trivial_group_target);
        serializer(segments);
        serializer(events);
        serializer(injection_cmode);
        serializer(production_cmode);
        serializer(filtrate_conc);
        serializer(perf_data);
        serializer(primaryvar);
        serializer(alq_state);
        serializer(group_target);
        serializer(was_shut_before_action_applied);
    }

    bool operator==(const SingleWellState&) const;

    std::string name;
    std::reference_wrapper<const ParallelWellInfo<Scalar>> parallel_info;

    WellStatus status{WellStatus::OPEN};
    bool producer;
    PhaseUsageInfo<IndexTraits> pu;
    Scalar bhp{0};
    Scalar thp{0};
    Scalar pressure_first_connection{0};

    // thermal related
    Scalar temperature{0};
    Scalar energy_rate{0.};

    Scalar efficiency_scaling_factor{1.0};

    // filtration injection concentration
    Scalar filtrate_conc{0};

    std::array<Scalar,4> phase_mixing_rates{};
    enum RateIndices {
      dissolved_gas = 0,
      dissolved_gas_in_water = 1,
      vaporized_oil = 2,
      vaporized_water = 3
    };

    struct GroupTarget {
        std::string group_name;
        Scalar target_value;
        Group::ProductionCMode production_cmode {Group::ProductionCMode::NONE};
        Scalar target_value_fallback;
        Group::ProductionCMode production_cmode_fallback {Group::ProductionCMode::NONE};

        bool operator==(const GroupTarget& other) const {
            return (group_name == other.group_name 
                 && target_value == other.target_value 
                 && production_cmode == other.production_cmode 
                 && target_value_fallback == other.target_value_fallback 
                 && production_cmode_fallback == other.production_cmode_fallback);
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer) {
            serializer(group_name);
            serializer(target_value);
            serializer(production_cmode);
            serializer(target_value_fallback);
            serializer(production_cmode_fallback);
        }
    };

    std::vector<Scalar> well_potentials;
    std::vector<Scalar> productivity_index;
    std::vector<Scalar> implicit_ipr_a;
    std::vector<Scalar> implicit_ipr_b;
    std::vector<Scalar> surface_rates;
    std::vector<Scalar> reservoir_rates;
    std::vector<Scalar> prev_surface_rates;
    PerfData<Scalar> perf_data;
    bool trivial_group_target;
    std::optional<GroupTarget> group_target;
    SegmentState<Scalar> segments;
    Events events;
    WellInjectorCMode injection_cmode{WellInjectorCMode::CMODE_UNDEFINED};
    WellProducerCMode production_cmode{WellProducerCMode::CMODE_UNDEFINED};
    std::vector<Scalar> primaryvar;
    ALQState<Scalar> alq_state;
    // This is used to indicate whether the well was shut before applying an action
    // if it was SHUT, even the action set the well to OPEN, the data in the well state
    // is not well-defined. We do not use it to overwrite the current well state.
    bool was_shut_before_action_applied {false};

    /// Special purpose method to support dynamically rescaling a well's
    /// CTFs through WELPI.
    ///
    /// \param[in] new_perf_data New perforation data.  Only
    ///    PerforationData::connection_transmissibility_factor actually
    ///    used (overwrites existing internal values).
    void reset_connection_factors(const std::vector<PerforationData<Scalar>>& new_perf_data);
    void update_producer_targets(const Well& ecl_well, const SummaryState& st);
    void update_injector_targets(const Well& ecl_well, const SummaryState& st);
    /// \brief update the type of the well and the targets.
    ///
    /// This called after ACTIONX is executed to update well rates. The new status is
    /// in ecl_well and st.
    /// \return whether well was switched to a producer
    void update_type_and_targets(const Well& ecl_well, const SummaryState& st);
    bool updateStatus(WellStatus status);
    void init_timestep(const SingleWellState& other);
    void shut();
    void stop();
    void open();

    // The sum_xxx_rates() functions sum over all connection rates of pertinent
    // types. In the case of distributed wells this involves an MPI
    // communication.
    Scalar sum_solvent_rates() const;
    Scalar sum_polymer_rates() const;
    Scalar sum_brine_rates() const;
    Scalar sum_microbial_rates() const;
    Scalar sum_oxygen_rates() const;
    Scalar sum_urea_rates() const;
    Scalar sum_wat_mass_rates() const;

    Scalar sum_filtrate_rate() const;
    Scalar sum_filtrate_total() const;

private:
    Scalar sum_connection_rates(const std::vector<Scalar>& connection_rates) const;
};

}

#endif
