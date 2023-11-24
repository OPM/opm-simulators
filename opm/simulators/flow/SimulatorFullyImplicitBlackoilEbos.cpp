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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <boost/date_time.hpp>

#include <algorithm>
#include <sstream>
#include <string>

namespace Opm {

void outputReportStep(const SimulatorTimer& timer)
{
    std::ostringstream stepMsg;
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
    stepMsg.imbue(std::locale(std::locale::classic(), facet));
    stepMsg << "\nReport step " << std::setw(2) <<timer.currentStepNum()
          << "/" << timer.numSteps()
          << " at day " << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day)
          << "/" << (double)unit::convert::to(timer.totalTime(), unit::day)
          << ", date = " << timer.currentDateTime();
    OpmLog::info(stepMsg.str());
}

void outputTimestampFIP(const SimulatorTimer& timer,
                        const std::string& title,
                        const std::string& version)
{
    std::ostringstream ss;
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d %b %Y");
    ss.imbue(std::locale(std::locale::classic(), facet));
    ss << "\n                              **************************************************************************\n"
    << "  Balance  at" << std::setw(10) << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day) << "  Days"
    << " *" << std::setw(30) << title << "                                          *\n"
    << "  Report " << std::setw(4) << timer.reportStepNum() << "    " << timer.currentDateTime()
    << "  *                                             Flow  version " << std::setw(11) << version << "  *\n"
    << "                              **************************************************************************\n";
    OpmLog::note(ss.str());
}

void checkSerializedCmdLine(const std::string& current,
                            const std::string& stored)
{
    auto filter_strings = [](const std::vector<std::string>& input)
    {
        std::vector<std::string> output;
        output.reserve(input.size());
        std::copy_if(input.begin(), input.end(), std::back_inserter(output),
                     [](const std::string& line)
                     {
                        return line.compare(0, 11, "EclDeckFile") != 0 &&
                               line.compare(0, 8, "LoadStep") != 0 &&
                               line.compare(0, 9, "OutputDir") != 0 &&
                               line.compare(0, 8, "SaveFile") != 0 &&
                               line.compare(0, 8, "SaveStep") != 0;
                     });
        return output;
    };

    auto curr_strings = split_string(current, '\n');
    auto stored_strings = split_string(stored, '\n');
    std::sort(curr_strings.begin(), curr_strings.end());
    std::sort(stored_strings.begin(), stored_strings.end());
    curr_strings = filter_strings(curr_strings);
    stored_strings = filter_strings(stored_strings);

    std::vector<std::string> difference;
    std::set_symmetric_difference(stored_strings.begin(), stored_strings.end(),
                                  curr_strings.begin(), curr_strings.end(),
                                  std::back_inserter(difference));

    std::vector<std::string> only_stored, only_curr;
    if (!difference.empty()) {
        for (std::size_t i = 0; i < difference.size(); ) {
            auto stored_it = std::find(stored_strings.begin(),
                                       stored_strings.end(), difference[i]);
            auto pos = difference[i].find_first_of('=');
            if (i < difference.size() - 1 &&
                difference[i].compare(0, pos, difference[i+1], 0, pos) == 0) {
                if (stored_it == stored_strings.end()) {
                    std::swap(difference[i], difference[i+1]);
                }
                i += 2;
            } else {
                if (stored_it == stored_strings.end()) {
                    only_curr.push_back(difference[i]);
                } else {
                    only_stored.push_back(difference[i]);
                }
                difference.erase(difference.begin() + i);
            }
        }
        std::stringstream str;
        str << "Differences:\n";
        for (std::size_t i = 0; i < difference.size(); ++i) {
            str << '\t' << (i % 2 == 0 ? '-' : '+') << difference[i] << '\n';
        }
        if (!only_stored.empty()) {
            str << "Only in serialized parameters:\n";
            for (const std::string& line : only_stored)
                str << '\t' << line << '\n';
        }
        if (!only_curr.empty()) {
            str << "Only in current parameters:\n";
            for (const std::string& line : only_curr)
                str << '\t' << line << '\n';
        }
        OpmLog::warning("Command line parameters mismatch:\n" + str.str());
    }
}

} // namespace Opm
