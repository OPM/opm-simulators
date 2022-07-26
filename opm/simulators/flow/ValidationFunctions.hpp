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

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include <opm/simulators/flow/KeywordValidation.hpp>


namespace Opm
{
namespace KeywordValidation
{


void validateBRINE(const DeckKeyword& keyword, std::vector<ValidationError>& errors);



// This is a mapping between keyword names and small functions
// for validation of special keywords.
std::unordered_map<std::string, std::function<void(const DeckKeyword& keyword, std::vector<KeywordValidation::ValidationError>& errors)>> specialValidation() {
    return {{"BRINE", validateBRINE}};
};

}
}
