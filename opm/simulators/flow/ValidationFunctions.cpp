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
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

#include <fmt/format.h>

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

    // Special case since we support the parsing of the items, which can be UDAs.
    void validateGSATPROD(const Opm::DeckKeyword& keyword,
                          std::vector<Opm::KeywordValidation::ValidationError>& errors)
    {
        if (keyword.empty()) {
            return;
        }

        using Kw = Opm::ParserKeywords::GSATPROD;

        const auto& record = keyword.getRecord(0);

        const auto& resv = record.getItem<Kw::RES_FLUID_VOL_PRODUCTION_RATE>().get<Opm::UDAValue>(0);
        if (resv.is_defined()) {
            const auto& resv_val = resv.is<double>() ? fmt::format("{}", resv.get<double>()) : resv.get<std::string>();
            errors.emplace_back(Opm::KeywordValidation::ValidationError {
                true,
                keyword.location(),
                0,  // not relevant
                5,
                resv_val,
                std::string{"Reservoir volume rate is not supported and should be defaulted (1*)"}}
            );
        }

        const auto& gaslift = record.getItem<Kw::LIFT_GAS_SUPPLY_RATE>().get<Opm::UDAValue>(0);
        if (gaslift.is_defined()) {
            const auto& gaslift_val = gaslift.is<double>() ? fmt::format("{}", gaslift.get<double>()) : gaslift.get<std::string>();
            errors.emplace_back(Opm::KeywordValidation::ValidationError {
                true,
                keyword.location(),
                0,  // not relevant
                6,
                gaslift_val,
                std::string{"Gaslift rate is not supported and should be defaulted (1*)"}}
            );
        }

        const auto& calrate = record.getItem<Kw::MEAN_CALORIFIC_VALUE>().get<Opm::UDAValue>(0);
        if (calrate.is_defined()) {
            const auto& calrate_val = calrate.is<double>() ? fmt::format("{}", calrate.get<double>()) : calrate.get<std::string>();
            errors.emplace_back(Opm::KeywordValidation::ValidationError {
                true,
                keyword.location(),
                0,  // not relevant
                7,
                calrate_val,
                std::string{"Calorific rate is not used and should be defaulted (1*)"}}
            );
        }
    }
}

namespace Opm::KeywordValidation {

std::unordered_map<std::string, ValidationFunction>
specialValidation()
{
    return {{"BRINE", validateBRINE}, {"GSATPROD", validateGSATPROD}};
}

}
