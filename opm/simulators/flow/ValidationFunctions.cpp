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
#include <opm/input/eclipse/Parser/ParserKeywords/F.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

#include <fmt/format.h>

#include <cstddef>
#include <string>

namespace {
    void validateBRINE(const Opm::Deck&,
                       const Opm::DeckKeyword& keyword,
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

    // FACTLI carries one multiplier per equilibration region in a single item,
    // and the per-item validation only ever inspects the first value.  Flow
    // labels single phase cells as if every multiplier were one, so every
    // region has to be checked here.
    void validateFACTLI(const Opm::Deck&,
                        const Opm::DeckKeyword& keyword,
                        std::vector<Opm::KeywordValidation::ValidationError>& errors)
    {
        if (keyword.empty()) {
            return;
        }

        using Kw = Opm::ParserKeywords::FACTLI;

        const auto& item = keyword.getRecord(0).getItem<Kw::DATA>();

        for (std::size_t i = 0; i < item.data_size(); ++i) {
            if (item.defaultApplied(i) || item.get<double>(i) == Kw::DATA::defaultValue) {
                continue;
            }

            errors.emplace_back(Opm::KeywordValidation::ValidationError {
                true,
                keyword.location(),
                0,  // a single record
                1,
                fmt::format("{}", item.get<double>(i)),
                fmt::format("FACTLI(DATA): the Li phase labelling is not scaled, so only "
                            "the default multiplier of 1.0 is supported (equilibration "
                            "region {})", i + 1)}
            );
        }
    }

    // PARACHOR only feeds the surface tension calculation, which is requested
    // with MISCIBLE.  Flow computes no surface tensions either way, so report
    // the keyword - but say plainly when the values could not have been used
    // anyway, so a deck that merely carries them is not read as a problem.
    void validatePARACHOR(const Opm::Deck& deck,
                          const Opm::DeckKeyword& keyword,
                          std::vector<Opm::KeywordValidation::ValidationError>& errors)
    {
        if (keyword.empty()) {
            return;
        }

        const bool miscible = deck.hasKeyword<Opm::ParserKeywords::MISCIBLE>();

        errors.emplace_back(Opm::KeywordValidation::ValidationError {
            miscible,
            keyword.location(),
            0,  // not relevant
            0,  // not relevant
            std::nullopt,
            miscible
                ? std::string{"Surface tensions are not calculated, so the parachors "
                              "MISCIBLE asks for are not used"}
                : std::string{"Surface tensions are not calculated; without MISCIBLE the "
                              "parachors are unused either way"}}
        );
    }

    // Special case since we support the parsing of the items, which can be UDAs.
    void validateGSATPROD(const Opm::Deck&,
                          const Opm::DeckKeyword& keyword,
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
    return {{"BRINE", validateBRINE},
            {"FACTLI", validateFACTLI},
            {"GSATPROD", validateGSATPROD},
            {"PARACHOR", validatePARACHOR}};
}

}
