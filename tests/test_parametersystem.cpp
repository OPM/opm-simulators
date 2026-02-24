/*
  Copyright 2024 Equinor.

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

#include <opm/models/utils/parametersystem.hpp>

#define BOOST_TEST_MODULE ParameterSystemTest
#include <boost/test/unit_test.hpp>

namespace Opm::Parameters {

struct SimpleParamBool { static constexpr bool value = false; };
struct SimpleParamDouble { static constexpr double value = 1.0; };
struct SimpleParamFloat { static constexpr float value = 2.0; };
struct SimpleParamInt { static constexpr int value = 1; };
struct SimpleParamString  { static constexpr auto value = "foo"; };
struct SimpleParamBoolN2
{
  static constexpr auto name = "SimpleB2";
  static constexpr bool value = true;
};

}

namespace {

struct Fixture
{
    Fixture()
    {
        Opm::Parameters::reset();
        Opm::Parameters::Register<Opm::Parameters::SimpleParamBool>("Simple bool parameter");
        Opm::Parameters::Register<Opm::Parameters::SimpleParamBoolN2>("Simpler bool parameter");
        Opm::Parameters::Register<Opm::Parameters::SimpleParamDouble>("Simple double parameter");
        Opm::Parameters::Register<Opm::Parameters::SimpleParamFloat>("Simple float parameter");
        Opm::Parameters::Register<Opm::Parameters::SimpleParamInt>("Simple int parameter");
        Opm::Parameters::Register<Opm::Parameters::SimpleParamString>("Simple string parameter");
        Opm::Parameters::SetDefault<Opm::Parameters::SimpleParamInt>(10);
        Opm::Parameters::Hide<Opm::Parameters::SimpleParamInt>();
        Opm::Parameters::endRegistration();
    }
};

std::string trimString(const std::string& input)
{
    std::string result(input);
    std::erase(result, ' ');
    std::erase(result, '\n');
    return result;
}

}

BOOST_FIXTURE_TEST_CASE(GetLists, Fixture)
{
  const char* argv[] = {
      "test_parametersystem",
      "--simple-param-bool=true",
      "--simple-param-float=3.0",
      "--simple-param-string=bar",
      "--unused-param=foo",
  };

    auto noPositional = [](std::function<void(const std::string&,
                                              const std::string&)>,
                           std::set<std::string>&,
                           std::string&,
                           int,
                           const char**,
                           int,
                           int) -> int
                        {
                            assert("Should not be here!");
                            return 0;
                        };


  Opm::Parameters::parseCommandLineOptions(5, argv, noPositional);

  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamBool>(), true);
  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamFloat>(), true);
  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamString>(), true);
  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamBoolN2>(), false);
  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamDouble>(), false);
  BOOST_CHECK_EQUAL(Opm::Parameters::IsSet<Opm::Parameters::SimpleParamInt>(), false);

  using SettingMap = std::vector<Opm::Parameters::Parameter>;

  const SettingMap set_ref = {
      {"SimpleParamBool", "true"},
      {"SimpleParamFloat", "3.0"},
      {"SimpleParamString", "bar"},
  };

  const SettingMap unused_ref = {
      {"UnusedParam", "foo"},
  };

  SettingMap  set, unused;
  Opm::Parameters::getLists(set, unused);

  BOOST_CHECK_EQUAL_COLLECTIONS(set.begin(), set.end(),
                                set_ref.begin(), set_ref.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(unused.begin(), unused.end(),
                                unused_ref.begin(), unused_ref.end());
}

BOOST_FIXTURE_TEST_CASE(ParseParameterFile, Fixture)
{
  Opm::Parameters::parseParameterFile("parametersystem.ini", true);

  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamBool>(), true);
  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamFloat>(), 3.f);
  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamString>(), "bar");
  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamBoolN2>(), true);
  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamDouble>(), 1.0);
  BOOST_CHECK_EQUAL(Opm::Parameters::Get<Opm::Parameters::SimpleParamInt>(), 10);
}

BOOST_FIXTURE_TEST_CASE(PrintUsage, Fixture)
{
  std::stringstream usage;
  Opm::Parameters::printUsage("", usage);
  BOOST_CHECK_EQUAL(trimString(usage.str()),
trimString(R"(
Recognized options:
    --simple-b2=BOOLEAN                           Simpler bool parameter. Default: true
    --simple-param-bool=BOOLEAN                   Simple bool parameter. Default: false
    --simple-param-double=SCALAR                  Simple double parameter. Default: 1
    --simple-param-float=SCALAR                   Simple float parameter. Default: 2
    --simple-param-string=STRING                  Simple string parameter. Default: "foo"
)"));
}

BOOST_FIXTURE_TEST_CASE(PrintUsageAll, Fixture)
{
  std::stringstream usage;
  Opm::Parameters::printUsage("===foobar===", usage, "", true);
  BOOST_CHECK_EQUAL(trimString(usage.str()),
trimString(R"(===foobar===
Recognized options:
    -h,--help                                     Print this help message and exit
    --help-all                                    Print all parameters, including obsolete, hidden and deprecated ones.
    --simple-b2=BOOLEAN                           Simpler bool parameter. Default: true
    --simple-param-bool=BOOLEAN                   Simple bool parameter. Default: false
    --simple-param-double=SCALAR                  Simple double parameter. Default: 1
    --simple-param-float=SCALAR                   Simple float parameter. Default: 2
    --simple-param-int=INTEGER                    Simple int parameter. Default: 10
    --simple-param-string=STRING                  Simple string parameter. Default: "foo"
)"));
}

BOOST_FIXTURE_TEST_CASE(PrintValues, Fixture)
{
  std::stringstream values;
  Opm::Parameters::printValues(values);
  BOOST_CHECK_EQUAL(trimString(values.str()),
trimString(R"(# [parameters which were specified at compile-time]
SimpleB2="1"
SimpleParamBool="0"
SimpleParamDouble="1"
SimpleParamFloat="2"
SimpleParamInt="10"
SimpleParamString="foo"
)"));
}
