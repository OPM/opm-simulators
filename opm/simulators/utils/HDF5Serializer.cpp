/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include <config.h>
#include <opm/simulators/utils/HDF5Serializer.hpp>

#include <algorithm>
#include <cstdlib>

namespace Opm {

void HDF5Serializer::writeHeader(const std::string& simulator_name,
                                 const std::string& module_version,
                                 const std::string& time_stamp,
                                 const std::string& case_name,
                                 const std::string& params,
                                 int num_procs)
{
    try {
        this->pack(simulator_name, module_version, time_stamp,
                   case_name, params, num_procs);
    } catch (...) {
        m_packSize = std::numeric_limits<std::size_t>::max();
        throw;
    }
    m_h5file.write("/", "simulator_info", m_buffer, HDF5File::DataSetMode::ROOT_ONLY);
}

int HDF5Serializer::lastReportStep() const
{
    const auto entries = m_h5file.list("/report_step");
    int last = -1;
    for (const auto& entry : entries) {
        int num = std::atoi(entry.c_str());
        last = std::max(last, num);
    }

    return last;
}

std::vector<int> HDF5Serializer::reportSteps() const
{
    const auto entries = m_h5file.list("/report_step");
    std::vector<int> result(entries.size());
    std::transform(entries.begin(), entries.end(), result.begin(),
                   [](const std::string& input)
                   {
                      return std::atoi(input.c_str());
                   });
    std::sort(result.begin(), result.end());
    return result;
}

}
