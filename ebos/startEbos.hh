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
#ifndef EWOMS_STARTEBOS_HH
#define EWOMS_STARTEBOS_HH

#include "ebos.hh"
#include <opm/models/utils/start.hh>

#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_ECL_INPUT
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#endif

#include <opm/parser/eclipse/Utility/String.hpp>

BEGIN_PROPERTIES

// forward declaration of property tags
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(ThreadManager);
NEW_PROP_TAG(PrintProperties);
NEW_PROP_TAG(PrintParameters);
NEW_PROP_TAG(ParameterFile);
NEW_PROP_TAG(Problem);
END_PROPERTIES

//! \cond SKIP_THIS

namespace Opm {

enum class FileOutputMode {
    //! \brief No output to files.
    OUTPUT_NONE = 0,
    //! \brief Output only to log files, no eclipse output.
    OUTPUT_LOG_ONLY = 1,
    //! \brief Output to all files.
    OUTPUT_ALL = 3
};


static void ensureOutputDirExists(const std::string& cmdline_output_dir)
{
    if (!Opm::filesystem::is_directory(cmdline_output_dir)) {
        try {
            Opm::filesystem::create_directories(cmdline_output_dir);
        }
        catch (...) {
            throw std::runtime_error("Creation of output directory '" + cmdline_output_dir + "' failed\n");
        }
    }
}

// Setup the OpmLog backends
static FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

    if (!cmdline_output_dir.empty()) {
        ensureOutputDirExists(cmdline_output_dir);
    }

    // create logFile
    using Opm::filesystem::path;
    path fpath(deck_filename);
    std::string baseName;
    std::ostringstream debugFileStream;
    std::ostringstream logFileStream;

    // Strip extension "." or ".DATA"
    std::string extension = uppercase(fpath.extension().string());
    if (extension == ".DATA" || extension == ".") {
        baseName = uppercase(fpath.stem().string());
    } else {
        baseName = uppercase(fpath.filename().string());
    }

    std::string output_dir = cmdline_output_dir;
    if (output_dir.empty()) {
        output_dir = absolute(path(baseName).parent_path()).string();
    }

    logFileStream << output_dir << "/" << baseName;
    debugFileStream << output_dir << "/" << baseName;

    if (mpi_rank_ != 0) {
        // Added rank to log file for non-zero ranks.
        // This prevents message loss.
        debugFileStream << "." << mpi_rank_;
        // If the following file appears then there is a bug.
        logFileStream << "." << mpi_rank_;
    }
    logFileStream << ".PRT";
    debugFileStream << ".DBG";

    FileOutputMode output;
    {
        static std::map<std::string, FileOutputMode> stringToOutputMode =
            { {"none", FileOutputMode::OUTPUT_NONE },
              {"false", FileOutputMode::OUTPUT_LOG_ONLY },
              {"log", FileOutputMode::OUTPUT_LOG_ONLY },
              {"all" , FileOutputMode::OUTPUT_ALL },
              {"true" , FileOutputMode::OUTPUT_ALL }};
        auto outputModeIt = stringToOutputMode.find(cmdline_output);
        if (outputModeIt != stringToOutputMode.end()) {
            output = outputModeIt->second;
        }
        else {
            output = FileOutputMode::OUTPUT_ALL;
            std::cerr << "Value " << cmdline_output <<
                " is not a recognized output mode. Using \"all\" instead."
                      << std::endl;
        }
    }

    if (output > FileOutputMode::OUTPUT_NONE) {
        std::shared_ptr<Opm::EclipsePRTLog> prtLog = std::make_shared<Opm::EclipsePRTLog>(logFileStream.str(), Opm::Log::NoDebugMessageTypes, false, output_cout_);
        Opm::OpmLog::addBackend("ECLIPSEPRTLOG", prtLog);
        prtLog->setMessageLimiter(std::make_shared<Opm::MessageLimiter>());
        prtLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(false));
    }

    if (output >= FileOutputMode::OUTPUT_LOG_ONLY) {
        std::string debugFile = debugFileStream.str();
        std::shared_ptr<Opm::StreamLog> debugLog = std::make_shared<Opm::EclipsePRTLog>(debugFileStream.str(), Opm::Log::DefaultMessageTypes, false, output_cout_);
        Opm::OpmLog::addBackend("DEBUGLOG", debugLog);
    }

    std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(std::cout, Opm::Log::StdoutMessageTypes);
    Opm::OpmLog::addBackend(stdout_log_id, streamLog);
    streamLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(true));

    return output;
}

//! \endcond

/*!
 * \ingroup Common
 *
 * \brief Wrapper around the main function that set up the OPM
 *        logging (.PRT, .DBG) for ebos.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments
 * \param argv The array of the command line arguments
 */
template <class TypeTag>
static inline int startEbos(int argc, char **argv)
{

    int myRank = 0;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
        myRank = Dune::Fem::MPIManager::rank();
#else
        myRank = Dune::MPIHelper::instance(argc, argv).rank();
#endif

        int paramStatus = setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*registerParams=*/true);
        if (paramStatus == 1)
            return 1;
        if (paramStatus == 2)
            return 0;

    bool outputCout = false;
    if (myRank == 0)
        outputCout = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);

    std::string deckFilename = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
    setupLogging(myRank,
                 deckFilename,
                 EWOMS_GET_PARAM(TypeTag, std::string, OutputDir),
                 EWOMS_GET_PARAM(TypeTag, std::string, OutputMode),
                 outputCout, "STDOUT_LOGGER");

    // Call the main function. Parameters are already registered
    // They should not be registered again
    return Opm::start<TypeTag>(argc, argv, /*registerParams=*/false);

}

} // namespace Opm

#endif
