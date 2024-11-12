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
#include <config.h>
#include <opm/simulators/flow/ValidationFunctions.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

namespace {
    void validateBRINE(const Opm::DeckKeyword& keyword,
                       std::vector<Opm::KeywordValidation::ValidationError>& errors)
    {
        if (keyword.empty()) {
            return;
        }

        errors.emplace_back(Opm::KeywordValidation::ValidationError {
            false,
            keyword.location(),
            0,  // not relevant
            0,  // not relevant
            std::nullopt,
            std::string{"The BRINE keyword does not accept any salt name arguments"}}
        );
    }
}

namespace Opm::KeywordValidation {

std::unordered_map<std::string, ValidationFunction>
specialValidation()
{
    return {{"BRINE", validateBRINE}};
}

}
