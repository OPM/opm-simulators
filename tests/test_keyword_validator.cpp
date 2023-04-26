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

#include <string>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/simulators/flow/KeywordValidation.hpp>

#define BOOST_TEST_MODULE KeywordValidatorTest

#include <boost/test/unit_test.hpp>

using namespace Opm;
using namespace Opm::KeywordValidation;

const UnsupportedKeywords test_unsupported_keywords = {
    {"ECHO", {false, "This is not a critical error"}},
    {"NOECHO", {true, std::nullopt}},
};


const PartiallySupportedKeywords<std::string> test_string_items = {
    {
        "PINCH",
        {
            {2, {false, allow_values<std::string> {"GAP", "BAR"}, std::nullopt}}, // GAP
            {4, {true, allow_values<std::string> {"TOPBOT"}, "This is a critical error"}}, // PINCHOUT_OPTION
        },
    },
    {
        "EQLOPTS",
        {
            {1, {false, allow_values<std::string> {}, std::nullopt}},
            {3, {false, allow_values<std::string> {"THPRESS"}, std::nullopt}},
        },
    },
};


const PartiallySupportedKeywords<int> test_int_items = {
    {
        "ENDSCALE",
        {
            {3, {false, allow_values<int> {1}, std::nullopt}}, // NTENDP
            {4, {true, allow_values<int> {20, 30, 40}, std::nullopt}}, // NSENDP
        },
    },
    {
        "COMPDAT",
        {
            {2, {false, allow_values<int> {1}, std::nullopt}}, // I
            {3, {false, allow_values<int> {1}, std::nullopt}}, // J
        },
    },
};


const PartiallySupportedKeywords<double> test_double_items = {
    {
        "EHYSTR",
        {
            {1, {false, [](double x) { return x > 0; }, std::nullopt}},
            {3, {false, allow_values<double> {1.0, 2.0}, std::nullopt}},
        },
    },
};



BOOST_AUTO_TEST_CASE(non_critical_keyword)
{
    const auto keywords_string = std::string {R"(
ECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(!errors[0].critical);
    BOOST_CHECK(!errors[0].item_number);
    BOOST_CHECK(!errors[0].item_value);
    BOOST_CHECK(*(errors[0].user_message) == "This is not a critical error");
}


BOOST_AUTO_TEST_CASE(critical_keyword)
{
    const auto keywords_string = std::string {R"(
NOECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["NOECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors[0].critical);
    BOOST_CHECK(!errors[0].item_number);
    BOOST_CHECK(!errors[0].item_value);
    BOOST_CHECK(!errors[0].user_message);
}


BOOST_AUTO_TEST_CASE(non_critical_keyword_item_string)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   FOO   1*   1* /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["PINCH"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(!errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 2);
    BOOST_CHECK(errors[0].item_value == "FOO");
    BOOST_CHECK(!errors[0].user_message);
}


BOOST_AUTO_TEST_CASE(critical_keyword_item_string)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   GAP   1*   FOO /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["PINCH"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 4);
    BOOST_CHECK(errors[0].item_value == "FOO");
    BOOST_CHECK(*(errors[0].user_message) == "This is a critical error");
}


BOOST_AUTO_TEST_CASE(non_critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(!errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 3);
    BOOST_CHECK(errors[0].item_value == "0");
    BOOST_CHECK(!errors[0].user_message);
}


BOOST_AUTO_TEST_CASE(critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  1 0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 4);
    BOOST_CHECK(errors[0].item_value == "0");
    BOOST_CHECK(!errors[0].user_message);
}

BOOST_AUTO_TEST_CASE(non_critical_keyword_item_double_ok)
{
    const auto keywords_string = std::string {R"(
EHYSTR
   1.0 0 1.0/
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EHYSTR"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors.size() == 0);
}

BOOST_AUTO_TEST_CASE(non_critical_keyword_item_double)
{
    const auto keywords_string = std::string {R"(
EHYSTR
   1.0 0 10.0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EHYSTR"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(!errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 3);
    BOOST_CHECK(errors[0].item_value == "10.000000");
    BOOST_CHECK(!errors[0].user_message);
}

BOOST_AUTO_TEST_CASE(non_critical_keyword_item_double_lambda)
{
    const auto keywords_string = std::string {R"(
EHYSTR
   -1.0 0 1.0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EHYSTR"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(!errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 1);
    BOOST_CHECK(errors[0].item_value == "-1.000000");
    BOOST_CHECK(!errors[0].user_message);
}

BOOST_AUTO_TEST_CASE(two_keyword_errors)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   GAP   1*   FOO /
ENDSCALE
   NODIR   REVERS  1 0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword1 = deck["PINCH"].back();
    const auto& test_keyword2 = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword1, errors);
    validator.validateDeckKeyword(test_keyword2, errors);
    BOOST_CHECK(errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 4);
    BOOST_CHECK(errors[0].item_value == "FOO");
    BOOST_CHECK(*(errors[0].user_message) == "This is a critical error");
    BOOST_CHECK(errors[1].critical);
    BOOST_CHECK(errors[1].item_number == 4);
    BOOST_CHECK(errors[1].item_value == "0");
    BOOST_CHECK(!errors[1].user_message);
}


BOOST_AUTO_TEST_CASE(report_not_critical)
{
    const auto keywords_string = std::string {R"(
ECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  ECHO: keyword not supported\n"
                   "  In file: <memory string>, line 2\n"
                   "  This is not a critical error");
}


BOOST_AUTO_TEST_CASE(report_critical_missing)
{
    const auto keywords_string = std::string {R"(
ECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(report_critical)
{
    const auto keywords_string = std::string {R"(
NOECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["NOECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  NOECHO: keyword not supported\n"
                   "  In file: <memory string>, line 2");
}


BOOST_AUTO_TEST_CASE(report_not_critical_missing)
{
    const auto keywords_string = std::string {R"(
NOECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["NOECHO"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_string)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   FOO   1*   1* /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["PINCH"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  PINCH: invalid value 'FOO' for item 2\n"
                   "  In file: <memory string>, line 2");
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_string_two_records)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   GAP   1*   1* /
COMPDAT
   C-4H   0   1   1    1    OPEN     1*    4.088536E+1   0.21600   1*   0.00000   1*   'Z' /
   C-4H   1   0   2    2    OPEN     1*    7.072475E+1   0.21600   1*   0.00000   1*   'Z' /
/
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["COMPDAT"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  COMPDAT: invalid value '0' in record 1 for item 2\n"
                   "  In file: <memory string>, line 4\n\n"
                   "  COMPDAT: invalid value '0' in record 2 for item 3\n"
                   "  In file: <memory string>, line 4");
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_string_missing)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   FOO   1*   1* /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["PINCH"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  ENDSCALE: invalid value '0' for item 3\n"
                   "  In file: <memory string>, line 2");
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_int_missing)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  1 0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  ENDSCALE: invalid value '0' for item 4\n"
                   "  In file: <memory string>, line 2");
}


BOOST_AUTO_TEST_CASE(error_report_critical_keyword_item_int_missing)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  1 0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_all_not_critical)
{
    const auto keywords_string = std::string {R"(
ECHO
NOECHO
PINCH
   0.41   GAP   1*   FOO /
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword1 = deck["ECHO"].back();
    const auto& test_keyword2 = deck["NOECHO"].back();
    const auto& test_keyword3 = deck["PINCH"].back();
    const auto& test_keyword4 = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword1, errors);
    validator.validateDeckKeyword(test_keyword2, errors);
    validator.validateDeckKeyword(test_keyword3, errors);
    validator.validateDeckKeyword(test_keyword4, errors);
    const auto report = get_error_report(errors, true, false);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  ECHO: keyword not supported\n"
                   "  In file: <memory string>, line 2\n"
                   "  This is not a critical error\n\n"
                   "  ENDSCALE: invalid value '0' for item 3\n"
                   "  In file: <memory string>, line 6");
}

BOOST_AUTO_TEST_CASE(error_report_not_critical)
{
    const auto keywords_string = std::string {R"(
ECHO
NOECHO
PINCH
   0.41   GAP   1*   FOO /
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword1 = deck["ECHO"].back();
    const auto& test_keyword2 = deck["NOECHO"].back();
    const auto& test_keyword3 = deck["PINCH"].back();
    const auto& test_keyword4 = deck["ENDSCALE"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword1, errors);
    validator.validateDeckKeyword(test_keyword2, errors);
    validator.validateDeckKeyword(test_keyword3, errors);
    validator.validateDeckKeyword(test_keyword4, errors);
    const auto report = get_error_report(errors, false, true);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  NOECHO: keyword not supported\n"
                   "  In file: <memory string>, line 3\n\n"
                   "  PINCH: invalid value 'FOO' for item 4\n"
                   "  In file: <memory string>, line 4\n"
                   "  This is a critical error");
}


BOOST_AUTO_TEST_CASE(keyword_without_default)
{
    const auto keywords_string = std::string {R"(
EQLOPTS
   1* /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EQLOPTS"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors.size() == 0);
}



BOOST_AUTO_TEST_CASE(keyword_without_default_with_wrong_value)
{
    const auto keywords_string = std::string {R"(
EQLOPTS
   FOO /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EQLOPTS"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors.size() == 1);
    BOOST_CHECK(errors[0].item_number == 1);
    BOOST_CHECK(errors[0].item_value == "FOO");
    BOOST_CHECK(!errors[0].user_message);
}



BOOST_AUTO_TEST_CASE(keyword_without_default_with_correct_value)
{
    const auto keywords_string = std::string {R"(
EQLOPTS
   1* 1* THPRESS /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck["EQLOPTS"].back();
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items, test_double_items, {});
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors.size() == 0);
}
