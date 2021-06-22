/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#include "config.h"

#if HAVE_MPI
#include "mpi.h"
#endif
#include "readDeck.hpp"

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/utility/String.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/FileSystem.hpp>

#include <opm/io/eclipse/EclIOdata.hpp>

#include <opm/output/eclipse/RestartIO.hpp>
#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/RestartFileView.hpp>
#include <opm/io/eclipse/rst/aquifer.hpp>
#include <opm/io/eclipse/rst/state.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>

#include "UnsupportedFlowKeywords.hpp"
#include "PartiallySupportedFlowKeywords.hpp"
#include <opm/simulators/flow/KeywordValidation.hpp>

#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <fmt/format.h>

#include <cstdlib>
#include <memory>
#include <utility>

namespace Opm
{

void ensureOutputDirExists_(const std::string& cmdline_output_dir)
{
    if (!Opm::filesystem::is_directory(cmdline_output_dir)) {
        try {
            Opm::filesystem::create_directories(cmdline_output_dir);
        }
        catch (...) {
            throw std::runtime_error(fmt::format("Creation of output directory '{}' failed\n", cmdline_output_dir));
        }
    }
}

// Setup the OpmLog backends
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

    if (!cmdline_output_dir.empty()) {
        ensureOutputDirExists_(cmdline_output_dir);
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
        output_dir = fpath.has_parent_path()
            ? absolute(fpath.parent_path()).generic_string()
            : Opm::filesystem::current_path().generic_string();
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

    if (mpi_rank_ == 0) {
        std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(std::cout, Opm::Log::StdoutMessageTypes);
        Opm::OpmLog::addBackend(stdout_log_id, streamLog);

        bool use_color_coding = OpmLog::stdoutIsTerminal();
        streamLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(use_color_coding));
    }
    return output;
}


namespace {
void setupMessageLimiter(const Opm::MessageLimits msgLimits,  const std::string& stdout_log_id) {
    std::shared_ptr<Opm::StreamLog> stream_log = Opm::OpmLog::getBackend<Opm::StreamLog>(stdout_log_id);

    const std::map<int64_t, int> limits = {{Opm::Log::MessageType::Note,
                                            msgLimits.getCommentPrintLimit()},
                                           {Opm::Log::MessageType::Info,
                                            msgLimits.getMessagePrintLimit()},
                                           {Opm::Log::MessageType::Warning,
                                            msgLimits.getWarningPrintLimit()},
                                           {Opm::Log::MessageType::Error,
                                            msgLimits.getErrorPrintLimit()},
                                           {Opm::Log::MessageType::Problem,
                                            msgLimits.getProblemPrintLimit()},
                                           {Opm::Log::MessageType::Bug,
                                            msgLimits.getBugPrintLimit()}};
    stream_log->setMessageLimiter(std::make_shared<Opm::MessageLimiter>(10, limits));
}
}


void readDeck(int rank, std::string& deckFilename, std::unique_ptr<Opm::Deck>& deck, std::unique_ptr<Opm::EclipseState>& eclipseState,
              std::unique_ptr<Opm::Schedule>& schedule, std::unique_ptr<UDQState>& udqState, std::unique_ptr<Opm::SummaryConfig>& summaryConfig,
              std::unique_ptr<ErrorGuard> errorGuard, std::shared_ptr<Opm::Python>& python, std::unique_ptr<ParseContext> parseContext,
              bool initFromRestart, bool checkDeck, const std::optional<int>& outputInterval)
{
    if (!errorGuard)
    {
        errorGuard = std::make_unique<ErrorGuard>();
    }

    int parseSuccess = 1; // > 0 is success
    std::string failureMessage;

    if (rank==0) {
        try
        {
            if ( (!deck || !schedule || !summaryConfig ) && !parseContext)
            {
                OPM_THROW(std::logic_error, "We need a parse context if deck, schedule, or summaryConfig are not initialized");
            }

            if (!deck)
            {
                Opm::Parser parser;
                deck = std::make_unique<Opm::Deck>( parser.parseFile(deckFilename , *parseContext, *errorGuard));

                Opm::KeywordValidation::KeywordValidator keyword_validator(
                    Opm::FlowKeywordValidation::unsupportedKeywords(),
                    Opm::FlowKeywordValidation::partiallySupported<std::string>(),
                    Opm::FlowKeywordValidation::partiallySupported<int>(),
                    Opm::FlowKeywordValidation::partiallySupported<double>());
                keyword_validator.validateDeck(*deck, *parseContext, *errorGuard);

                if ( checkDeck )
                    Opm::checkDeck(*deck, parser, *parseContext, *errorGuard);
            }

            if (!eclipseState) {
#if HAVE_MPI
                eclipseState = std::make_unique<Opm::ParallelEclipseState>(*deck);
#else
                eclipseState = std::make_unique<Opm::EclipseState>(*deck);
#endif
            }

            const auto& init_config = eclipseState->getInitConfig();
            if (init_config.restartRequested()) {
                // Analytic aquifers must always be loaded from the restart
                // file in restarted runs and the corresponding keywords
                // (e.g., AQUANCON and AQUCT) do not exist in the input deck
                // in this case.  In other words, there's no way to check if
                // there really are analytic aquifers in the run until we
                // attempt to read the specifications from the restart file.
                // If the loader determines that there are no analytic
                // aquifers, then 'EclipseState::loadRestartAquifers()' does
                // nothing.
                const int report_step = init_config.getRestartStep();
                const auto rst_filename = eclipseState->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                auto rst_file = std::make_shared<EclIO::ERst>(rst_filename);
                auto rst_view = std::make_shared<EclIO::RestartFileView>(std::move(rst_file), report_step);

                // Note: RstState::load() will just *read* from the grid
                // structure, and only do so if the case actually includes
                // analytic aquifers.  The pointer to the input grid is just
                // to allow 'nullptr' to signify "don't load aquifers" in
                // certain unit tests.  Passing an optional<EclipseGrid> is
                // too expensive however since doing so will create a copy
                // of the grid inside the optional<>.
                const auto rst_state = RestartIO::RstState::
                    load(std::move(rst_view), &eclipseState->getInputGrid());

                eclipseState->loadRestartAquifers(rst_state.aquifers);

                // For the time being initializing wells and groups from the
                // restart file is not possible.  Work is underway and the
                // ability is included here contingent on user-level switch
                // 'initFromRestart' (i.e., setting "--sched-restart=false"
                // as a command line invocation parameter).
                const auto* init_state = initFromRestart ? &rst_state : nullptr;
                if (!schedule) {
                    schedule = std::make_unique<Schedule>(*deck, *eclipseState,
                                                          *parseContext, *errorGuard,
                                                          python, outputInterval, init_state);
                }

                udqState = std::make_unique<UDQState>((*schedule)[0].udq().params().undefinedValue());
                udqState->load_rst(rst_state);
            }
            else {
                if (!schedule) {
                    schedule = std::make_unique<Schedule>(*deck, *eclipseState,
                                                          *parseContext, *errorGuard,
                                                          python);
                }

                udqState = std::make_unique<UDQState>((*schedule)[0].udq().params().undefinedValue());
            }


            if (Opm::OpmLog::hasBackend("STDOUT_LOGGER")) {
                // loggers might not be set up!
                setupMessageLimiter((*schedule)[0].message_limits(), "STDOUT_LOGGER");
            }

            if (!summaryConfig) {
                summaryConfig = std::make_unique<Opm::SummaryConfig>(*deck, *schedule, eclipseState->fieldProps(),
                                                                     eclipseState->aquifer(), *parseContext, *errorGuard);
            }

            Opm::checkConsistentArrayDimensions(*eclipseState, *schedule, *parseContext, *errorGuard);
        }
        catch(const OpmInputError& input_error) {
            failureMessage = input_error.what();
            parseSuccess = 0;
        }
        catch(const std::exception& std_error)
        {
            failureMessage = std_error.what();
            parseSuccess = 0;
        }
    }
#if HAVE_MPI
    else {
        if (!summaryConfig)
            summaryConfig = std::make_unique<Opm::SummaryConfig>();
        if (!schedule)
            schedule = std::make_unique<Opm::Schedule>(python);
        if (!eclipseState)
            eclipseState = std::make_unique<Opm::ParallelEclipseState>();
    }

    try
    {
        Opm::eclStateBroadcast(*eclipseState, *schedule, *summaryConfig);
    }
    catch(const std::exception& broadcast_error)
    {
        failureMessage = broadcast_error.what();
        OpmLog::error(fmt::format("Distributing properties to all processes failed\n"
                                  "Internal error message: {}", broadcast_error.what()));
        parseSuccess = 0;
    }

#endif

    if (*errorGuard) { // errors encountered
        parseSuccess = 0;
        errorGuard->dump();
        errorGuard->clear();
    }


    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    parseSuccess = comm.min(parseSuccess);

    if (!parseSuccess)
    {
        if (rank == 0)
        {
            OpmLog::error("Unrecoverable errors were encountered while loading input");
        }
#if HAVE_MPI
        MPI_Finalize();
#endif
        std::exit(EXIT_FAILURE);
    }
}
} // end namespace Opm
