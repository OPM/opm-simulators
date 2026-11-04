/*
  Copyright 2026, SINTEF Digital

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
/*!
 * \file
 * \brief Helpers for tests that run a deck through the simulator's main entry point
 */
#ifndef OPM_RUN_DECK_HELPER_HPP
#define OPM_RUN_DECK_HELPER_HPP

#include "config.h"

#include "TemporaryDirectory.hpp"

#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/io/eclipse/ESmry.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \brief Global fixture for tests that use runDeck()
 *
 * Initializes MPI and sets up communication. Unlike SimulatorFixture, it adds
 * no log backend, since each run sets up its own logging and removes all
 * backends when it ends.
 */
struct RunDeckFixture
{
    RunDeckFixture()
    {
        auto& suite = boost::unit_test::framework::master_test_suite();
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(suite.argc, suite.argv);
#else
        Dune::MPIHelper::instance(suite.argc, suite.argv);
#endif
        FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());
    }
};

/*!
 * \brief Run a deck through the simulator's main entry point
 *
 * Writes the deck to <dir>/<name>.DATA, runs it with its output in dir and
 * returns its summary with all vectors loaded. A failed run ends the test case.
 *
 * \tparam TypeTag The type tag of the simulator to run
 */
template <class TypeTag>
EclIO::ESmry runDeck(const std::filesystem::path& dir,
                     const std::string& name,
                     const std::string& deck_text)
{
    const auto deck_file = dir / (name + ".DATA");
    std::ofstream(deck_file) << deck_text;

    std::vector<std::string> args = {boost::unit_test::framework::master_test_suite().argv[0],
                                     deck_file.string(),
                                     "--output-dir=" + dir.string()};
    std::vector<char*> argv;
    std::ranges::transform(args, std::back_inserter(argv),
                           [](std::string& arg) { return arg.data(); });
    argv.push_back(nullptr);
    Parameters::reset();
    const int status = Main(static_cast<int>(args.size()), argv.data(), /*ownMPI=*/false)
                           .runStatic<TypeTag>();
    BOOST_REQUIRE_EQUAL(status, EXIT_SUCCESS);

    EclIO::ESmry smry(dir / (name + ".SMSPEC"));
    smry.loadData();
    return smry;
}

} // namespace Opm

#endif // OPM_RUN_DECK_HELPER_HPP
