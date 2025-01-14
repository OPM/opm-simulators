/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef OPM_BLACKOILMODEL_CONV_MONITOR_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_CONV_MONITOR_HEADER_INCLUDED

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>

namespace Opm {

/// Implementation of penalty cards for three-phase black oil.
template <class Scalar>
class BlackoilModelConvergenceMonitor
{
public:
    using MonitorParams = typename BlackoilModelParameters<Scalar>::ConvergenceMonitorParams;
    explicit BlackoilModelConvergenceMonitor(const MonitorParams& param);

    void checkPenaltyCard(ConvergenceReport& report, int iteration);

    void reset();

private:
    const MonitorParams& param_;
    ConvergenceReport::PenaltyCard total_penaltyCard_;
    double prev_distance_;
    int prev_above_tolerance_;
};

} // namespace Opm

#endif
