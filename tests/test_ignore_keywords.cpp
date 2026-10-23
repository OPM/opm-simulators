/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE TestIgnoreKeywords
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#include <opm/common/utility/OpmInputError.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Fixture {
    Fixture()
    {
        // UNKNOWNA (RUNSPEC) and UNKNOWNB (SCHEDULE) are not known to the parser.
        const std::string deck = R"(RUNSPEC
DIMENS
  2 2 1 /
START
  1 JAN 2020 /
UNKNOWNA
  1 2 3 /
GRID
DXV
  2*100.0 /
DYV
  2*100.0 /
DZV
  10.0 /
DEPTHZ
  9*2000.0 /
PERMX
  4*100.0 /
PERMY
  4*100.0 /
PERMZ
  4*10.0 /
PORO
  4*0.3 /
PROPS
SOLUTION
SCHEDULE
UNKNOWNB
TSTEP
  10 /
END
)";

        input_path = std::filesystem::temp_directory_path() / "ignore_keywords_test";
        std::filesystem::remove_all(input_path);
        std::filesystem::create_directories(input_path);

        deck_file = (input_path / "INPUT.DATA").string();
        std::ofstream(deck_file) << deck;
    }

    ~Fixture()
    {
        std::filesystem::remove_all(input_path);
    }

    int initialize(const std::string& option) const
    {
        std::vector<const char*> argv {"test_ignore_keywords", deck_file.c_str()};
        if (!option.empty()) {
            argv.push_back(option.c_str());
        }
        const int argc = static_cast<int>(argv.size());
        argv.push_back(nullptr);

        Opm::Parameters::reset();
        Opm::ThreadManager::registerParameters();

        Opm::Main main(argc, const_cast<char**>(argv.data()), false);
        return main.justInitialize();
    }

    std::filesystem::path input_path;
    std::string deck_file;
};

} // namespace

BOOST_FIXTURE_TEST_CASE(UnknownKeywordsFail, Fixture)
{
    BOOST_CHECK_EQUAL(initialize(""), EXIT_FAILURE);
}

BOOST_FIXTURE_TEST_CASE(OnlyListedKeywordsAreIgnored, Fixture)
{
    BOOST_CHECK_EQUAL(initialize("--ignore-keywords=UNKNOWNA"), EXIT_FAILURE);
}

BOOST_FIXTURE_TEST_CASE(IgnoredKeywordsParse, Fixture)
{
    BOOST_CHECK_EQUAL(initialize("--ignore-keywords=UNKNOWNA:UNKNOWNB"), EXIT_SUCCESS);
}

BOOST_FIXTURE_TEST_CASE(SpacesAroundIgnoredKeywords, Fixture)
{
    BOOST_CHECK_EQUAL(initialize("--ignore-keywords= UNKNOWNA : UNKNOWNB "), EXIT_SUCCESS);
}

BOOST_FIXTURE_TEST_CASE(RunspecParseIgnoresListedKeywords, Fixture)
{
    // flow_comp and flowexp_comp parse RUNSPEC alone before reading the deck.
    const auto parseRunspec = [this](const std::string& ignored) {
        Opm::ParseContext parseContext;
        Opm::ignoreKeywords(parseContext, ignored);
        Opm::Parser {}.parseFile(deck_file, parseContext,
                                 std::vector {Opm::Ecl::SectionType::RUNSPEC});
    };
    BOOST_CHECK_THROW(parseRunspec(""), Opm::OpmInputError);
    BOOST_CHECK_NO_THROW(parseRunspec("UNKNOWNA"));
}

BOOST_FIXTURE_TEST_CASE(VanguardReadDeckIgnoresListedKeywords, Fixture)
{
    // The flowexp variants read the deck through the vanguard.
    const auto readDeck = [this](const std::string& ignored) {
        Opm::FlowGenericVanguard::modelParams_ = {};
        Opm::FlowGenericVanguard::readDeck(deck_file, ignored);
    };
    BOOST_CHECK_THROW(readDeck(""), std::runtime_error);
    BOOST_CHECK_NO_THROW(readDeck("UNKNOWNA:UNKNOWNB"));
}

bool
init_unit_test_func()
{
    return true;
}

int
main(int argc, char** argv)
{
    int argcDummy = 1;
    const char* tmp[] = {"test_ignore_keywords", nullptr};
    char** argvDummy = const_cast<char**>(tmp);
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
    Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

    Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
