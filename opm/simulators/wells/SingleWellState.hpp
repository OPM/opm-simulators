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

#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

namespace Opm {

struct PerforationData;
class SummaryState;
class Well;

class SingleWellState {
public:
    SingleWellState(const std::string& name,
                    const ParallelWellInfo& pinfo,
                    bool is_producer,
                    double presssure_first_connection,
                    const std::vector<PerforationData>& perf_input,
                    const PhaseUsage& pu,
                    double temp);

    static SingleWellState serializationTestObject(const ParallelWellInfo& pinfo);

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(name);
        serializer(status);
        serializer(producer);
        serializer(bhp);
        serializer(thp);
        serializer(temperature);
        serializer(phase_mixing_rates);
        serializer(well_potentials);
        serializer(productivity_index);
        serializer(implicit_ipr_a);
        serializer(implicit_ipr_b);        
        serializer(surface_rates);
        serializer(reservoir_rates);
        serializer(prev_surface_rates);
        serializer(trivial_target);
        serializer(segments);
        serializer(events);
        serializer(injection_cmode);
        serializer(production_cmode);
        serializer(filtrate_conc);
        serializer(perf_data);
    }

    bool operator==(const SingleWellState&) const;

    std::string name;
    std::reference_wrapper<const ParallelWellInfo> parallel_info;

    WellStatus status{WellStatus::OPEN};
    bool producer;
    PhaseUsage pu;
    double bhp{0};
    double thp{0};
    double temperature{0};

    // filtration injection concentration
    double filtrate_conc{0};

    std::array<double,4> phase_mixing_rates{};
    enum RateIndices {
      dissolved_gas = 0,
      dissolved_gas_in_water = 1,
      vaporized_oil = 2,
      vaporized_water = 3
    };

    std::vector<double> well_potentials;
    std::vector<double> productivity_index;
    std::vector<double> implicit_ipr_a;
    std::vector<double> implicit_ipr_b;    
    std::vector<double> surface_rates;
    std::vector<double> reservoir_rates;
    std::vector<double> prev_surface_rates;
    PerfData perf_data;
    bool trivial_target;
    SegmentState segments;
    Events events;
    WellInjectorCMode injection_cmode{WellInjectorCMode::CMODE_UNDEFINED};
    WellProducerCMode production_cmode{WellProducerCMode::CMODE_UNDEFINED};


    /// Special purpose method to support dynamically rescaling a well's
    /// CTFs through WELPI.
    ///
    /// \param[in] new_perf_data New perforation data.  Only
    ///    PerforationData::connection_transmissibility_factor actually
    ///    used (overwrites existing internal values).
    void reset_connection_factors(const std::vector<PerforationData>& new_perf_data);
    void update_producer_targets(const Well& ecl_well, const SummaryState& st);
    void update_injector_targets(const Well& ecl_well, const SummaryState& st);
    void update_targets(const Well& ecl_well, const SummaryState& st);
    void updateStatus(WellStatus status);
    void init_timestep(const SingleWellState& other);
    void shut();
    void stop();
    void open();

    // The sum_xxx_rates() functions sum over all connection rates of pertinent
    // types. In the case of distributed wells this involves an MPI
    // communication.
    double sum_solvent_rates() const;
    double sum_polymer_rates() const;
    double sum_brine_rates() const;

    double sum_filtrate_rate() const;
    double sum_filtrate_total() const;

private:
    double sum_connection_rates(const std::vector<double>& connection_rates) const;
};


}



#endif
