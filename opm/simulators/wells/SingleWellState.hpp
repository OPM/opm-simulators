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

#include <vector>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>

namespace Opm {

class SingleWellState {
public:
    SingleWellState(bool is_producer, std::size_t num_phases, double temp);

    Well::Status status{Well::Status::OPEN};
    bool producer;
    double bhp{0};
    double thp{0};
    double temperature{};
    double dissolved_gas_rate{0};
    double vaporized_oil_rate{0};
    std::vector<double> well_potentials;
    std::vector<double> productivity_index;
    Events events;
    Well::InjectorCMode injection_cmode{Well::InjectorCMode::CMODE_UNDEFINED};
    Well::ProducerCMode production_cmode{Well::ProducerCMode::CMODE_UNDEFINED};


    void init_timestep(const SingleWellState& other);
    void shut();
    void stop();
};


}



#endif
