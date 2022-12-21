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

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Events.hpp>

#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

namespace Opm {

struct PerforationData;

class SingleWellState {
public:
    SingleWellState(const std::string& name,
                    const ParallelWellInfo& pinfo,
                    bool is_producer,
                    double presssure_first_connection,
                    const std::vector<PerforationData>& perf_input,
                    const PhaseUsage& pu,
                    double temp);

    std::string name;
    std::reference_wrapper<const ParallelWellInfo> parallel_info;

    Well::Status status{Well::Status::OPEN};
    bool producer;
    PhaseUsage pu;
    double bhp{0};
    double thp{0};
    double temperature{0};
    double dissolved_gas_rate{0};
    double dissolved_gas_rate_in_water{0};
    double vaporized_oil_rate{0};
    double vaporized_wat_rate{0};
    std::vector<double> well_potentials;
    std::vector<double> productivity_index;
    std::vector<double> surface_rates;
    std::vector<double> reservoir_rates;
    PerfData perf_data;
    bool trivial_target;
    SegmentState segments;
    Events events;
    Well::InjectorCMode injection_cmode{Well::InjectorCMode::CMODE_UNDEFINED};
    Well::ProducerCMode production_cmode{Well::ProducerCMode::CMODE_UNDEFINED};


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
    void updateStatus(Well::Status status);
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
private:
    double sum_connection_rates(const std::vector<double>& connection_rates) const;
};


}



#endif
