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

#ifndef OPM_PERFDATA_HEADER_INCLUDED
#define OPM_PERFDATA_HEADER_INCLUDED

#include <vector>

#include <opm/core/props/BlackoilPhases.hpp>

namespace Opm
{

class PerfData
{
private:
    PhaseUsage pu;

public:
    PerfData(std::size_t num_perf, const PhaseUsage& pu);
    std::size_t size() const;

    std::vector<double> pressure;
    std::vector<double> rates;
    std::vector<double> phase_rates;
    std::vector<double> solvent_rates;
    std::vector<double> polymer_rates;
    std::vector<double> brine_rates;
    std::vector<double> water_throughput;
    std::vector<double> skin_pressure;
    std::vector<double> water_velocity;
    std::vector<double> prod_index;
};

} // namespace Opm

#endif // OPM_PERFORATIONDATA_HEADER_INCLUDED
