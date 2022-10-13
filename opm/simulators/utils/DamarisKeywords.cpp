/*
  Copyright 2021 Equinor.

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

#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <string>
#include <map>

/*
    Below is the Damaris Keywords supported by Damaris to be filled
    in the built-in XML file.

    The entries in the map below will be filled by the corresponding
    Damaris Keywords. Yet, only output directory and FileMode are to
    be chosen by the user
*/

namespace Opm::DamarisOutput
{
std::map<std::string, std::string>
DamarisKeywords(std::string OutputDir, bool enableDamarisOutputCollective)
{
    if (enableDamarisOutputCollective) {
        std::map<std::string, std::string> damaris_keywords = {
            {"_SHMEM_BUFFER_BYTES_REGEX_", "536870912"},
            {"_DC_REGEX_", "1"},
            {"_DN_REGEX_", "0"},
            {"_File_Mode", "Collective"},
            {"_MORE_VARIABLES_REGEX_", ""},
            {"_PATH_REGEX_", OutputDir},
            {"_MYSTORE_OR_EMPTY_REGEX_", "MyStore"},
        };
        return damaris_keywords;
    } else {
        std::map<std::string, std::string> damaris_keywords = {
            {"_SHMEM_BUFFER_BYTES_REGEX_", "536870912"},
            {"_DC_REGEX_", "1"},
            {"_DN_REGEX_", "0"},
            {"_File_Mode", "FilePerCore"},
            {"_MORE_VARIABLES_REGEX_", ""},
            {"_PATH_REGEX_", OutputDir},
            {"_MYSTORE_OR_EMPTY_REGEX_", "MyStore"},
        };
        return damaris_keywords;
    }
}

} // namespace Opm::DamarisOutput
