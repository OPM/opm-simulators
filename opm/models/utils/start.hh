// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
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
 * \brief Provides convenience routines to bring up the simulation at runtime.
 */
#ifndef EWOMS_START_HH
#define EWOMS_START_HH

#include <dune/common/parallel/mpihelper.hh>

#include <opm/material/common/ResetLocale.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/terminal.hpp>
#include <opm/models/utils/simulator.hh>

#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#if HAVE_MPI
#include <mpi.h>
#endif

//! \cond SKIP_THIS

namespace Opm {
/*!
 * \brief Announce all runtime parameters to the registry but do not specify them yet.
 */
template <class TypeTag>
static inline void registerAllParameters_(bool finalizeRegistration)
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using TM = GetPropType<TypeTag, Properties::ThreadManager>;

    Parameters::Register<Parameters::ParameterFile>
        ("An .ini file which contains a set of run-time parameters");
    Parameters::Register<Parameters::PrintParameters>
        ("Print the values of the run-time parameters at the "
         "start of the simulation");

    TM::registerParameters();
    Simulator::registerParameters();

    if (finalizeRegistration) {
        Parameters::endRegistration();
    }
}

/*!
 * \brief Register all runtime parameters, parse the command line
 *        arguments and the parameter file.
 *
 * \param argc The number of command line arguments
 * \param argv Array with the command line argument strings
 * \return A negative value if --help or --print-properties was provided,
 *         a positive value for errors or 0 for success.
 */
template <class TypeTag>
static inline int setupParameters_(int argc,
                                   const char **argv,
                                   bool registerParams,
                                   bool allowUnused,
                                   bool handleHelp,
                                   const int myRank)
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    ////////////////////////////////////////////////////////////
    // Register all parameters
    ////////////////////////////////////////////////////////////
    if (registerParams) {
        registerAllParameters_<TypeTag>(true);
    }

    ////////////////////////////////////////////////////////////
    // set the parameter values
    ////////////////////////////////////////////////////////////

    // fill the parameter tree with the options from the command line
    const auto& positionalParamCallback = Problem::handlePositionalParameter;
    std::string helpPreamble; // print help if non-empty!
    if (myRank == 0 && handleHelp) {
        helpPreamble = Problem::helpPreamble(argc, argv);
    }
    const std::string s =
        Parameters::parseCommandLineOptions(argc,
                                            argv,
                                            positionalParamCallback,
                                            helpPreamble);
    if (!s.empty()) {
        int status = 1;
        if (s == "Help called") { // only on master process
            status = -1; // Use negative values to indicate --help argument
        }
#if HAVE_MPI
        // Force -1 if the master process has that.
        int globalStatus;
        MPI_Allreduce(&status, &globalStatus, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        return globalStatus;
#endif
        return status;
    }

    const std::string paramFileName = Parameters::Get<Parameters::ParameterFile>(false);
    if (!paramFileName.empty()) {
        ////////////////////////////////////////////////////////////
        // add the parameters specified using an .ini file
        ////////////////////////////////////////////////////////////

        // check whether the parameter file is readable.
        if (!Parameters::parseParameterFile(paramFileName, /*overwrite=*/false)) {
            std::ostringstream oss;
            if (myRank == 0) {
                oss << "Parameter file \"" << paramFileName
                    << "\" does not exist or is not readable.";
                Parameters::printUsage(argv[0], std::cerr, oss.str());
            }
            return /*status=*/1;
        }
    }

    // make sure that no unknown parameters are encountered
    using ParamList = std::vector<Parameters::Parameter>;

    ParamList usedParams;
    ParamList unusedParams;

    Parameters::getLists(usedParams, unusedParams);
    if (!allowUnused && !unusedParams.empty()) {
        if (myRank == 0) {
            if (unusedParams.size() == 1) {
                std::cerr << "The following explicitly specified parameter is unknown:\n";
            }
            else {
                std::cerr << "The following " << unusedParams.size()
                          << " explicitly specified parameters are unknown:\n";
            }

            std::cerr << "\n";
            for (const auto& keyValue : unusedParams) {
                std::cerr << "   " << keyValue << "\n";
            }
            std::cerr << "\n";

            std::cerr << "Use\n"
                      << "\n"
                      << "  " << argv[0] << " --help\n"
                      << "\n"
                      <<"to obtain the list of recognized command line parameters.\n\n";
        }
        return /*status=*/1;
    }

    return /*status=*/0;
}

/*!
 * \ingroup Common
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file and runs the simulation
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments
 * \param argv The array of the command line arguments
 */
template <class TypeTag>
static inline int start(int argc, char **argv,  bool registerParams)
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using TM = GetPropType<TypeTag, Properties::ThreadManager>;

    assignResetTerminalSignalHandlers();

    resetLocale();

    int myRank = 0;
    try {
        // initialize MPI, finalize is done automatically on exit
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
        myRank = Dune::Fem::MPIManager::rank();
#else
        myRank = Dune::MPIHelper::instance(argc, argv).rank();
#endif
        const int paramStatus =
            setupParameters_<TypeTag>(argc,
                                      const_cast<const char**>(argv),
                                      registerParams,
                                      false,
                                      true,
                                      myRank);
        if (paramStatus == 1) {
            return 1;
        }
        if (paramStatus == 2) {
            return 0;
        }

        TM::init();

        // setting up backend for STDCOUT logger
        if (myRank == 0) {
            setupStreamLogging("STDOUT_LOGGER");
        }

        // read the initial time step and the end time
        const Scalar endTime = Parameters::Get<Parameters::EndTime<Scalar>>();
        if (endTime < -1e50) {
            if (myRank == 0) {
                Parameters::printUsage(argv[0], std::cerr,
                                       "Mandatory parameter '--end-time' not specified!");
            }
            return 1;
        }

        const Scalar initialTimeStepSize = Parameters::Get<Parameters::InitialTimeStepSize<Scalar>>();
        if (initialTimeStepSize < -1e50) {
            if (myRank == 0) {
                Parameters::printUsage(argv[0], std::cerr,
                                       "Mandatory parameter '--initial-time-step-size' "
                                       "not specified!");
            }
            return 1;
        }

        if (myRank == 0) {
#ifdef EWOMS_VERSION
            const std::string versionString = EWOMS_VERSION;
#else
            const std::string versionString;
#endif
            const std::string briefDescription = Problem::briefDescription();
            if (!briefDescription.empty()) {
                const std::string tmp = breakLines(briefDescription,
                                                   /*indentWidth=*/0,
                                                   getTtyWidth());
                std::cout << tmp << std::endl << std::endl;
            }
            else {
                std::cout << "opm models " << versionString
                          << " will now start the simulation. " << std::endl;
            }
        }

        // print the parameters if requested
        const int printParams = Parameters::Get<Parameters::PrintParameters>();
        if (myRank == 0) {
            const std::string endParametersSeparator("# [end of parameters]\n");
            if (printParams) {
                bool printSeparator = false;
                if (printParams == 1 || !isatty(fileno(stdout))) {
                    Parameters::printValues(std::cout);
                    printSeparator = true;
                }
                else {
                    // always print the list of specified but unused parameters
                    printSeparator = printSeparator || Parameters::printUnused(std::cout);
                }
                if (printSeparator) {
                    std::cout << endParametersSeparator;
                }
            }
            else {
                // always print the list of specified but unused parameters
                if (Parameters::printUnused(std::cout)) {
                    std::cout << endParametersSeparator;
                }
            }
        }

        // instantiate and run the concrete problem. make sure to
        // deallocate the problem and before the time manager and the
        // grid
        Simulator simulator;
        simulator.run();

        if (myRank == 0) {
            std::cout << "Simulation completed" << std::endl;                                 
        }
        return 0;
    }
    catch (std::exception& e) {
        if (myRank == 0) {
            std::cout << e.what() << ". Abort!\n" << std::flush;

            std::cout << "Trying to reset TTY.\n";
            resetTerminal();
        }

        return 1;
    }
    catch (...) {
        if (myRank == 0) {
            std::cout << "Unknown exception thrown!\n" << std::flush;

            std::cout << "Trying to reset TTY.\n";
            resetTerminal();
        }

        return 3;
    }
}

} // namespace Opm

#endif
