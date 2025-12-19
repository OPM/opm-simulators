// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Common fixtures and utilities for simulator-based tests
 */
#ifndef OPM_SIMULATOR_FIXTURE_HPP
#define OPM_SIMULATOR_FIXTURE_HPP

#include "config.h"
#include "TestTypeTag.hpp"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <memory>
#include <string>

namespace Opm {

/*!
 * \brief Global fixture for simulator-based tests
 *
 * Initializes MPI and sets up communication. Use with BOOST_GLOBAL_FIXTURE.
 */
struct SimulatorFixture
{
    SimulatorFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
#else
        Dune::MPIHelper::instance(argc, argv);
#endif
        FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());
    }
};

/*!
 * \brief Initialize a simulator from an Eclipse deck file
 *
 * \tparam TypeTag The type tag for the simulator
 * \param filename Path to the Eclipse deck file
 * \param test_name Name of the test (used in argv[0])
 * \param threads_per_process Number of threads per process (default: 1)
 * \return Unique pointer to the initialized simulator
 */
template <class TypeTag>
std::unique_ptr<GetPropType<TypeTag, Properties::Simulator>>
initSimulator(const char* filename,
              const char* test_name = "test_simulator",
              int threads_per_process = 1)
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    std::string filename_arg = "--ecl-deck-file-name=";
    filename_arg += filename;

    const char* argv[] = {
        test_name,
        filename_arg.c_str()
    };

    Parameters::reset();
    registerAllParameters_<TypeTag>(false);
    registerEclTimeSteppingParameters<double>();
    BlackoilModelParameters<double>::registerParameters();
    Parameters::Register<Parameters::EnableTerminalOutput>("Do *NOT* use!");
    Parameters::SetDefault<Parameters::ThreadsPerProcess>(threads_per_process);
    Parameters::endRegistration();
    setupParameters_<TypeTag>(/*argc=*/sizeof(argv) / sizeof(argv[0]),
                              argv,
                              /*registerParams=*/false,
                              /*allowUnused=*/false,
                              /*handleHelp=*/true,
                              /*myRank=*/0);

    FlowGenericVanguard::readDeck(filename);
    return std::make_unique<Simulator>();
}

} // namespace Opm

#endif // OPM_SIMULATOR_FIXTURE_HPP
