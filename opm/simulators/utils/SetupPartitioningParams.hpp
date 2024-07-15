/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_SETUP_PARTITIONING_PARAMS_HPP
#define OPM_SETUP_PARTITIONING_PARAMS_HPP

#include <map>
#include <string>

namespace Opm {

std::map<std::string,std::string> setupZoltanParams(const std::string& conf);
std::map<std::string,std::string> setupMetisParams(const std::string& conf);

} // namespace Opm

#endif // OPM_SETUP_PARTITIONING_PARAMS_HPP
