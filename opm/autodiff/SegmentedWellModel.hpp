/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_SEGMENTEDWELLMODEL_HEADER_INCLUDED
#define OPM_SEGMENTEDWELLMODEL_HEADER_INCLUDED

#include <vector>

struct Wells;

namespace Opm
{

    class WellState;
    class PhaseUsage;


    /// A class giving a well model, by which we mean a way to compute
    /// the pressure deltas of each perforation and the bottom-hole
    /// pressure. This class contains an explicit model, that uses a
    /// different density for each well segment, that is between each
    /// pair of perforations.
    class SegmentedWellModel
    {
    public:
        static std::vector<double> computeConnectionPressureDelta(const Wells& wells,
                                                                  const WellState& wstate,
                                                                  const PhaseUsage& phase_usage,
                                                                  const std::vector<double>& b_perf,
                                                                  const std::vector<double>& rsmax_perf,
                                                                  const std::vector<double>& rvmax_perf,
                                                                  const std::vector<double>& z_perf,
                                                                  const std::vector<double>& surf_dens,
                                                                  const double gravity);
    };

} // namespace Opm

#endif // OPM_SEGMENTEDWELLMODEL_HEADER_INCLUDED
