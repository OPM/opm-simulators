/*
  Copyright 2021 Equinor ASA.


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

#include <opm/simulators/wells/PerfData.hpp>

namespace Opm
{


PerfData::PerfData(std::size_t num_perf, const PhaseUsage& pu_arg):
    pu(pu_arg),
    pressure(num_perf),
    rates(num_perf),
    phase_rates(num_perf * pu.num_phases),
    solvent_rates(num_perf),
    polymer_rates(num_perf),
    brine_rates(num_perf),
    water_throughput(num_perf),
    skin_pressure(num_perf),
    water_velocity(num_perf),
    prod_index(num_perf)
{

}

}

