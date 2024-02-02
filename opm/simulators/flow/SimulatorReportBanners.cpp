/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#include <config.h>
#include <opm/simulators/flow/SimulatorReportBanners.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <boost/date_time.hpp>

#include <sstream>

namespace Opm::details {

void outputReportStep(const SimulatorTimer& timer)
{
    std::ostringstream stepMsg;
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
    stepMsg.imbue(std::locale(std::locale::classic(), facet));
    stepMsg << "\nReport step " << std::setw(2) << timer.currentStepNum()
            << "/" << timer.numSteps()
            << " at day " << unit::convert::to(timer.simulationTimeElapsed(), unit::day)
            << "/" << unit::convert::to(timer.totalTime(), unit::day)
            << ", date = " << timer.currentDateTime();
    OpmLog::info(stepMsg.str());
}

} // namespace Opm::details
