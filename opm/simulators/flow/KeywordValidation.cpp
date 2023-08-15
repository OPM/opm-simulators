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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <map>
#include <string>
#include <type_traits>

#include <fmt/format.h>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

namespace Opm
{



namespace KeywordValidation
{


    std::string get_error_report(const std::vector<ValidationError>& errors, const bool include_noncritical, const bool include_critical)
    {
        const std::string keyword_format = "  {keyword}: keyword not supported\n";
        const std::string item_format1 = "  {{keyword}}: invalid value '{}' for item {}\n";
        const std::string item_format2 = "  {{keyword}}: invalid value '{}' in record {} for item {}\n";
        const std::string location_format = "  In file: {file}, line {line}\n";

        std::string report;
        for (const ValidationError& err : errors) {
            if ((err.critical && include_critical) || (!err.critical && include_noncritical)) {
                if (err.item_number && err.item_value) {
                    std::string message;
                    if (err.record_number == 0) {
                        message = fmt::format(item_format1, *(err.item_value), *(err.item_number));
                    } else {
                        message = fmt::format(item_format2, *(err.item_value), err.record_number, *(err.item_number));
                    }
                    report.append(OpmInputError::format(message, err.location));
                } else {
                    report.append(OpmInputError::format(keyword_format, err.location));
                }
                report.append(OpmInputError::format(location_format, err.location));
                if (err.user_message) {
                    report.append("  " + *(err.user_message) + "\n");
                }
                report.append("\n");
            }
        }

        if (!report.empty()) {
            // Remove the last two newlines.
            report.erase(report.length() - 2);
            // Prepend header and file name.
            report.insert(0, "Unsupported keywords or keyword items:\n\n");
        }
        return report;
    }

    void
    KeywordValidator::validateDeck(const Deck& deck,
                                   const ParseContext& parse_context,
                                   const bool treat_critical_as_noncritical,
                                   ErrorGuard& error_guard) const
    {
        // Make a vector with all problems encountered in the deck.
        std::vector<ValidationError> errors;
        for (const auto& keyword : deck) {
            validateDeckKeyword(keyword, errors);

            const auto& special_it = this->m_special_validation.find(keyword.name());
            if (special_it != this->m_special_validation.end()) {
                const auto& validator = special_it->second;
                validator(keyword, errors);
            }
        }

        if (treat_critical_as_noncritical) {
            // Get both critical and noncritical errors.
            auto warning_report = get_error_report(errors, true, true);
            if (!warning_report.empty()) {
                // Report all as warnings.
                parse_context.handleError(ParseContext::SIMULATOR_KEYWORD_NOT_SUPPORTED, warning_report, std::nullopt, error_guard);
            }
        } else {
            // First report non-critical problems as a warning.
            auto warning_report = get_error_report(errors, true, false);
            if (!warning_report.empty()) {
                parse_context.handleError(
                                          ParseContext::SIMULATOR_KEYWORD_NOT_SUPPORTED, warning_report, std::nullopt, error_guard);
            }

            // Then report critical problems as an error.
            auto error_report = get_error_report(errors, false, true);
            if (!error_report.empty()) {
                OpmLog::info("\nOPM Flow will terminate due to unsupported critical keywords.\n"
                             "These are keywords that would change the simulator results if supported.\n"
                             "If you need to override this behaviour, you can use the command line argument --parsing-strictness=low,\n"
                             "this will reduce the severity of this to a warning instead of an error.");
                parse_context.handleError(
                                          ParseContext::SIMULATOR_KEYWORD_NOT_SUPPORTED_CRITICAL, error_report, std::nullopt, error_guard);
            }
        }
    }

    void KeywordValidator::validateDeckKeyword(const DeckKeyword& keyword, std::vector<ValidationError>& errors) const
    {
        const auto& it = m_keywords.find(keyword.name());
        if (it != m_keywords.end()) {
            // If the keyword is not supported, add an error for that.
            const auto& properties = it->second;
            errors.push_back(ValidationError {
                properties.critical, keyword.location(), 1, std::nullopt, std::nullopt, properties.message});
        } else {
            // Otherwise, check all its items.
            validateKeywordItems(keyword, m_string_items, errors);
            validateKeywordItems(keyword, m_int_items, errors);
            validateKeywordItems(keyword, m_double_items, errors);
        }
    }


    template <typename T>
    void KeywordValidator::validateKeywordItems(const DeckKeyword& keyword,
                                                const PartiallySupportedKeywords<T>& partially_supported_items,
                                                std::vector<ValidationError>& errors) const
    {
        const auto& keyword_properties = partially_supported_items.find(keyword.name());
        if (keyword_properties != partially_supported_items.end()) {
            // If this keyworcs has partially supported items, iterate over all of them.
            for (std::size_t record_index = 0; record_index < keyword.size(); record_index++) {
                const auto& record = keyword.getRecord(record_index);
                for (std::size_t item_index = 0; item_index < record.size(); item_index++) {
                    const auto& item = record.getItem(item_index);
                    // Find the index number, which starts counting at one, so item_index + 1
                    const auto& item_properties = keyword_properties->second.find(item_index + 1);
                    if (item_properties != keyword_properties->second.end()) {
                        if (item.hasValue(0)) {
                            // Validate the item, if it is partially supported.
                            validateKeywordItem<T>(keyword,
                                                   item_properties->second,
                                                   keyword.size() > 1,
                                                   record_index,
                                                   item_index,
                                                   item.get<T>(0),
                                                   errors);
                        }
                    }
                }
            }
        }
    }


    template <typename T>
    void KeywordValidator::validateKeywordItem(const DeckKeyword& keyword,
                                               const PartiallySupportedKeywordProperties<T>& properties,
                                               const bool multiple_records,
                                               const std::size_t record_index,
                                               const std::size_t item_index,
                                               const T& item_value,
                                               std::vector<ValidationError>& errors) const
    {
        if (!properties.validator(item_value)) {
            // If the value is not permitted, format the value to report it.
            std::string formatted_value;
            if constexpr (std::is_arithmetic<T>::value)
                formatted_value = std::to_string(item_value);
            else
                formatted_value = item_value;
            // Add the relevant information to the vector of errors. Record and
            // index numbers start at 1, so add 1. Pass zero for the record
            // index if there is only a single record.
            errors.push_back(ValidationError {properties.critical,
                                              keyword.location(),
                                              multiple_records ? record_index + 1 : 0,
                                              item_index + 1,
                                              formatted_value,
                                              properties.message});
        }
    }


} // namespace KeywordValidation

} // namespace Opm
