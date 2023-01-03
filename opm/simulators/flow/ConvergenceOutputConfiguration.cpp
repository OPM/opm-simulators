/*
  Copyright 2022 Equinor ASA.

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

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>

#include <algorithm>
#include <cstddef>
#include <regex>
#include <stdexcept>
#include <string_view>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {
    std::vector<std::string> tokenizeOptionValues(std::string_view options)
    {
        const auto split = std::regex { R"(\s*,\s*)" };
        return {
            std::cregex_token_iterator {
                options.begin(), options.end(), split, -1
            },
            std::cregex_token_iterator {}
        };
    }

    void reportUnsupportedOptionValuesAndThrow(std::vector<std::string> unsupp,
                                               std::string_view         optionName)
    {
        std::sort(unsupp.begin(), unsupp.end());
        auto u = std::unique(unsupp.begin(), unsupp.end());

        const auto pl = (std::distance(unsupp.begin(), u) != 1) ? "s" : "";

        if (optionName.empty()) {
            throw std::invalid_argument {
                fmt::format("Unsupported convergence output "
                            "option value{}: {}\n"
                            "Supported values are \"none\", "
                            "\"steps\", and \"iterations\"",
                            pl, fmt::join(unsupp.begin(), u, ", "))
            };
        }

        throw std::invalid_argument {
            fmt::format("Option {}:\n - Unsupported value{}: {}\n"
                        " - Supported values are \"none\", "
                        "\"steps\", and \"iterations\"",
                        optionName, pl,
                        fmt::join(unsupp.begin(), u, ", "))
        };
    }

    std::vector<Opm::ConvergenceOutputConfiguration::Option>
    getOptions(std::string_view options, std::string_view optionName)
    {
        using Option = Opm::ConvergenceOutputConfiguration::Option;

        auto opt = std::vector<Option>{};

        const auto values = std::unordered_map<std::string, Option> {
            { "none"      , Option::None       },
            { "step"      , Option::Steps      }, // Alias for 'steps' (plural)
            { "steps"     , Option::Steps      },
            { "iteration" , Option::Iterations }, // Alias for 'iterations' (plural)
            { "iterations", Option::Iterations },
        };

        auto unsupp = std::vector<std::string>{};
        for (const auto& value : tokenizeOptionValues(options)) {
            if (auto option = values.find(value); option != values.end()) {
                opt.push_back(option->second);
            }
            else {
                unsupp.push_back(value);
            }
        }

        if (! unsupp.empty()) {
            reportUnsupportedOptionValuesAndThrow(std::move(unsupp), optionName);
        }

        return opt;
    }
} // Anonymous namespace

// ===========================================================================
// Public Interface Below Separator
// ===========================================================================

Opm::ConvergenceOutputConfiguration::
ConvergenceOutputConfiguration(std::string_view options,
                               std::string_view optionName)
{
    auto is_none = false;
    for (const auto& option : getOptions(options, optionName)) {
        if (option == Option::None) {
            is_none = true;
            break;
        }

        this->flag_ |= static_cast<std::byte>(option);
    }

    if (is_none) {
        // Recall: "none" overrides all other options.
        this->flag_ = std::byte{0};
    }
}
