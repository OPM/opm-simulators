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

#include <string>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
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
            {2, {false, {"GAP"}, std::nullopt}}, // GAP
            {4, {true, {"TOPBOT"}, "This is a critical error"}}, // PINCHOUT_OPTION
        },
    },
};


const PartiallySupportedKeywords<int> test_int_items = {
    {
        "ENDSCALE",
        {
            {3, {false, {1}, std::nullopt}}, // NTENDP
            {4, {true, {20}, std::nullopt}}, // NSENDP
        },
    },
    {
        "COMPDAT",
        {
            {2, {false, {1}, std::nullopt}}, // I
            {3, {false, {1}, std::nullopt}}, // J
        },
    },
};


BOOST_AUTO_TEST_CASE(non_critical_keyword)
{
    const auto keywords_string = std::string {R"(
ECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck.getKeyword("ECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("NOECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("PINCH");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("PINCH");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    BOOST_CHECK(errors[0].critical);
    BOOST_CHECK(errors[0].item_number == 4);
    BOOST_CHECK(errors[0].item_value == "0");
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
    const auto& test_keyword1 = deck.getKeyword("PINCH");
    const auto& test_keyword2 = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
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
    const auto& test_keyword = deck.getKeyword("ECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword = deck.getKeyword("ECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(report_critical)
{
    const auto keywords_string = std::string {R"(
NOECHO
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck.getKeyword("NOECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true);
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
    const auto& test_keyword = deck.getKeyword("NOECHO");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_string)
{
    const auto keywords_string = std::string {R"(
PINCH
   0.41   FOO   1*   1* /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck.getKeyword("PINCH");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword = deck.getKeyword("COMPDAT");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword = deck.getKeyword("PINCH");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_non_critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  0 20 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true);
    BOOST_CHECK(report.empty());
}


BOOST_AUTO_TEST_CASE(error_report_critical_keyword_item_int)
{
    const auto keywords_string = std::string {R"(
ENDSCALE
   NODIR   REVERS  1 0 /
)"};
    const auto deck = Parser {}.parseString(keywords_string);
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, true);
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
    const auto& test_keyword = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword1 = deck.getKeyword("ECHO");
    const auto& test_keyword2 = deck.getKeyword("NOECHO");
    const auto& test_keyword3 = deck.getKeyword("PINCH");
    const auto& test_keyword4 = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword1, errors);
    validator.validateDeckKeyword(test_keyword2, errors);
    validator.validateDeckKeyword(test_keyword3, errors);
    validator.validateDeckKeyword(test_keyword4, errors);
    const auto report = get_error_report(errors, false);
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
    const auto& test_keyword1 = deck.getKeyword("ECHO");
    const auto& test_keyword2 = deck.getKeyword("NOECHO");
    const auto& test_keyword3 = deck.getKeyword("PINCH");
    const auto& test_keyword4 = deck.getKeyword("ENDSCALE");
    KeywordValidator validator(test_unsupported_keywords, test_string_items, test_int_items);
    std::vector<ValidationError> errors;
    validator.validateDeckKeyword(test_keyword1, errors);
    validator.validateDeckKeyword(test_keyword2, errors);
    validator.validateDeckKeyword(test_keyword3, errors);
    validator.validateDeckKeyword(test_keyword4, errors);
    const auto report = get_error_report(errors, true);
    BOOST_CHECK(report
                == "Unsupported keywords or keyword items:\n\n"
                   "  NOECHO: keyword not supported\n"
                   "  In file: <memory string>, line 3\n\n"
                   "  PINCH: invalid value 'FOO' for item 4\n"
                   "  In file: <memory string>, line 4\n"
                   "  This is a critical error");
}
