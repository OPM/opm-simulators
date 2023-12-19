/*
  Copyright 2023 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE TestOutputDir
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/Main.hpp>

#include <filesystem>
#include <fstream>
#include <string>

namespace {

struct Fixture {
    Fixture()
    {
        const std::string deck = R"(RUNSPEC
DIMENS
  10 10 3 /
START
 8 OCT 2020 /
GRID
DXV
  10*100.0 /
DYV
  10*100.0 /
DZV
  3*10.0 /
DEPTHZ
  121*2000.0 /
PERMX
  300*100.0 /
PERMY
  300*100.0 /
PERMZ
  300*10.0 /
PORO
  300*0.3 /
PROPS
SOLUTION
SCHEDULE

TSTEP
  10
/
END
)";

        input_path = std::filesystem::temp_directory_path() / "outputdir_test/";

        std::filesystem::remove_all(input_path);
        std::filesystem::create_directories(input_path / "subdir" / "subdir");

        for (const auto& file_path : {input_path / "INPUT.DATA",
                                      input_path / "subdir" / "INPUT.DATA",
                                      input_path / "subdir" / "subdir" / "INPUT.DATA"}) {
            std::ofstream of(file_path);
            of << deck;
        }
    }

    ~Fixture()
    {
        std::filesystem::remove_all(input_path);
    }

    std::filesystem::path input_path;
};

}

BOOST_FIXTURE_TEST_CASE(WithOutputDir, Fixture)
{
    std::filesystem::current_path(input_path);

    using PathPair = std::pair<std::filesystem::path, std::filesystem::path>;

    for (const auto& Case : {PathPair{input_path, input_path / "output1"},
                             PathPair{input_path / "subdir", input_path / "output2"},
                             PathPair{input_path / "subdir" / "subdir", input_path / "output3"}}) {
        const std::string output_path = "--output-dir=" + Case.second.string();
        const std::string input_file_path = (Case.first / "INPUT.DATA");
        const char* no_param[] = {"test_outputdir", input_file_path.c_str(),
                                  output_path.c_str(), nullptr};

        using ParamsMeta = Opm::GetProp<Opm::Properties::TTag::FlowEarlyBird,
                                        Opm::Properties::ParameterMetaData>;
        ParamsMeta::clear();

        Opm::Main main(3, const_cast<char**>(no_param), false);

        BOOST_CHECK_EQUAL(main.justInitialize(), EXIT_SUCCESS);
        BOOST_CHECK(!std::filesystem::exists(input_path / "INPUT.PRT"));
        BOOST_CHECK(!std::filesystem::exists(input_path / "INPUT.DBG"));
        BOOST_CHECK(!std::filesystem::exists(Case.first / "INPUT.PRT"));
        BOOST_CHECK(!std::filesystem::exists(Case.first / "INPUT.DBG"));
        BOOST_CHECK(std::filesystem::exists(Case.second / "INPUT.PRT"));
        BOOST_CHECK(std::filesystem::exists(Case.second / "INPUT.DBG"));
    }
}

BOOST_FIXTURE_TEST_CASE(NoOutputDir, Fixture)
{
    std::filesystem::current_path(input_path);

    for (const auto& Case : {input_path / "subdir" / "subdir",
                             input_path / "subdir"}) {
        const std::string input_file_path = (Case / "INPUT.DATA");
        const char* no_param[] = {"test_outputdir", input_file_path.c_str(), nullptr};

        using ParamsMeta = Opm::GetProp<Opm::Properties::TTag::FlowEarlyBird,
                                        Opm::Properties::ParameterMetaData>;
        ParamsMeta::clear();

        Opm::Main main(2, const_cast<char**>(no_param), false);

        BOOST_CHECK_EQUAL(main.justInitialize(), EXIT_SUCCESS);
        BOOST_CHECK(!std::filesystem::exists(input_path / "INPUT.PRT"));
        BOOST_CHECK(!std::filesystem::exists(input_path / "INPUT.DBG"));
        BOOST_CHECK(std::filesystem::exists(Case/ "INPUT.PRT"));
        BOOST_CHECK(std::filesystem::exists(Case/ "INPUT.DBG"));
    }
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    // MPI setup.
    int argcDummy = 1;
    const char *tmp[] = {"test_outputdir"};
    char **argvDummy = const_cast<char**>(tmp);
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
    Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

    Opm::EclGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
