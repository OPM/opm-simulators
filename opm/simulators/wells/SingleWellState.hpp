/*
  Copyright 2012, 2014 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 NORCE.
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_SINGLEWELLSTATE_HEADER_INCLUDED
#define OPM_SINGLEWELLSTATE_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

#include <array>
#include <vector>

namespace Opm
{

/// Encapsulate all data needed to represent the state of a well for persistence purposes,
/// i.e. from one timestep to the next or for restarting.
template <int NumActivePhases>
struct SingleWellState
{
    // ------ Types ------

    using PhaseRates = std::array<double, NumActivePhases>;

    struct Connection {
        double pressure = -1e100;
        PhaseRates surface_rates = {0.0};
        PhaseRates reservoir_rates = {0.0};
        double solvent_rate = 0.0;
        double water_throughput = 0.0;
        double skin_pressure = 0.0;
        double water_velocity = 0.0;
    };

    struct Segment {
        int segment_number = -1;
        double pressure = -1e100;
        double pressure_drop = 0.0;
        double pressure_drop_hydrostatic = 0.0;
        double pressure_drop_acceleration = 0.0;
        double pressure_drop_friction = 0.0;
        PhaseRates surface_rates = {0.0};
    };

    // ------ Data members ------

    // Flags and statuses.
    Well::Status status = Well::Status::SHUT;
    bool is_producer = true;
    bool effective_events_occurred = true;
    Well::InjectorCMode current_injection_control = Well::InjectorCMode::CMODE_UNDEFINED;
    Well::ProducerCMode current_production_control = Well::ProducerCMode::CMODE_UNDEFINED;

    // Quantities.
    double bhp = -1e100;
    double thp = -1e100;
    double temperature = 273.15 + 20; // 20 degrees Celcius by default.
    PhaseRates surface_rates = {0.0};
    PhaseRates reservoir_rates = {0.0};
    double dissolved_gas_rate = 0.0;
    double vaporized_oil_rate = 0.0;
    PhaseRates potentials = {0.0};
    PhaseRates productivity_index = {0.0};

    // Connection and segment data.
    std::vector<Connection> connections;
    std::vector<Segment> segments;
};


} // namespace Opm


#endif // OPM_SINGLEWELLSTATE_HEADER_INCLUDED
