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

#include <opm/simulators/utils/readDeck.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/io/eclipse/EclIOdata.hpp>
#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/RestartFileView.hpp>
#include <opm/io/eclipse/rst/aquifer.hpp>
#include <opm/io/eclipse/rst/state.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/EclipseState/checkDeck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldData.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/ArrayDimChecker.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/flow/KeywordValidation.hpp>
#include <opm/simulators/flow/ValidationFunctions.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/utils/PartiallySupportedFlowKeywords.hpp>
#include <opm/simulators/utils/UnsupportedFlowKeywords.hpp>

#include <fmt/format.h>

#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {

    void setupMessageLimiter(const Opm::MessageLimits& msgLimits,
                             const std::string& stdout_log_id)
    {
        const auto limits = std::map<std::int64_t, int> {
            {Opm::Log::MessageType::Note,     msgLimits.getCommentPrintLimit()},
            {Opm::Log::MessageType::Info,     msgLimits.getMessagePrintLimit()},
            {Opm::Log::MessageType::Warning,  msgLimits.getWarningPrintLimit()},
            {Opm::Log::MessageType::Error,    msgLimits.getErrorPrintLimit()},
            {Opm::Log::MessageType::Problem,  msgLimits.getProblemPrintLimit()},
            {Opm::Log::MessageType::Bug,      msgLimits.getBugPrintLimit()},
        };

        Opm::OpmLog::getBackend<Opm::StreamLog>(stdout_log_id)
            ->setMessageLimiter(std::make_shared<Opm::MessageLimiter>(10, limits));
    }

    void loadObjectsFromRestart(const Opm::Deck&                     deck,
                                const Opm::Parser&                   parser,
                                const Opm::ParseContext&             parseContext,
                                const bool                           initFromRestart,
                                const std::optional<int>&            outputInterval,
                                Opm::EclipseState&                   eclipseState,
                                std::shared_ptr<Opm::Python>         python,
                                std::shared_ptr<Opm::Schedule>&      schedule,
                                std::unique_ptr<Opm::UDQState>&      udqState,
                                std::unique_ptr<Opm::Action::State>& actionState,
                                std::unique_ptr<Opm::WellTestState>& wtestState,
                                Opm::ErrorGuard&                     errorGuard)
    {
        // Analytic aquifers must always be loaded from the restart file in
        // restarted runs and the corresponding keywords (e.g., AQUANCON and
        // AQUCT) do not exist in the input deck in this case.  In other
        // words, there's no way to check if there really are analytic
        // aquifers in the run until we attempt to read the specifications
        // from the restart file.  If the loader determines that there are
        // no analytic aquifers, then 'EclipseState::loadRestartAquifers()'
        // does nothing.
        const auto& init_config = eclipseState.getInitConfig();
        const int report_step = init_config.getRestartStep();
        const auto rst_filename = eclipseState.getIOConfig()
            .getRestartFileName(init_config.getRestartRootName(), report_step, false);

        auto rst_file = std::make_shared<Opm::EclIO::ERst>(rst_filename);
        auto rst_view = std::make_shared<Opm::EclIO::RestartFileView>
            (std::move(rst_file), report_step);

        // Note: RstState::load() will just *read* from the grid structure,
        // and only do so if the case actually includes analytic aquifers.
        // The pointer to the input grid is just to allow 'nullptr' to
        // signify "don't load aquifers" in certain unit tests.  Passing an
        // optional<EclipseGrid> is too expensive however since doing so
        // will create a copy of the grid inside the optional<>.
        const auto rst_state = Opm::RestartIO::RstState::
            load(std::move(rst_view),
                 eclipseState.runspec(), parser,
                 &eclipseState.getInputGrid());

        eclipseState.loadRestartAquifers(rst_state.aquifers);

        // For the time being initializing wells and groups from the restart
        // file is not possible.  Work is underway and the ability is
        // included here contingent on user-level switch 'initFromRestart'
        // (i.e., setting "--sched-restart=false" as a command line
        // invocation parameter).
        const auto* init_state = initFromRestart ? &rst_state : nullptr;
        if (schedule == nullptr) {
            schedule = std::make_shared<Opm::Schedule>
                (deck, eclipseState, parseContext, errorGuard,
                 std::move(python), outputInterval, init_state);
        }

        // Read network pressures from restart
        if (rst_state.network.isActive()) {
            eclipseState.loadRestartNetworkPressures(rst_state.network);
        }

        udqState = std::make_unique<Opm::UDQState>
            ((*schedule)[0].udq().params().undefinedValue());
        udqState->load_rst(rst_state);

        actionState = std::make_unique<Opm::Action::State>();
        actionState->load_rst((*schedule)[report_step].actions(), rst_state);

        wtestState = std::make_unique<Opm::WellTestState>(schedule->runspec().start_time(), rst_state);
    }

    void createNonRestartDynamicObjects(const Opm::Deck&                     deck,
                                        const Opm::EclipseState&             eclipseState,
                                        const Opm::ParseContext&             parseContext,
                                        std::shared_ptr<Opm::Python>         python,
                                        std::shared_ptr<Opm::Schedule>&      schedule,
                                        std::unique_ptr<Opm::UDQState>&      udqState,
                                        std::unique_ptr<Opm::Action::State>& actionState,
                                        std::unique_ptr<Opm::WellTestState>& wtestState,
                                        Opm::ErrorGuard&                     errorGuard)
    {
        if (schedule == nullptr) {
            schedule = std::make_shared<Opm::Schedule>
                (deck, eclipseState, parseContext,
                 errorGuard, std::move(python));
        }

        udqState = std::make_unique<Opm::UDQState>
            ((*schedule)[0].udq().params().undefinedValue());

        actionState = std::make_unique<Opm::Action::State>();
        wtestState = std::make_unique<Opm::WellTestState>();
    }

    Opm::Deck
    readDeckFile(const std::string&       deckFilename,
                 const bool               checkDeck,
                 const Opm::Parser&       parser,
                 const Opm::ParseContext& parseContext,
                 const bool               treatCriticalAsNonCritical,
                 Opm::ErrorGuard&         errorGuard)
    {
        Opm::Deck deck(parser.parseFile(deckFilename, parseContext, errorGuard));

        auto keyword_validator = Opm::KeywordValidation::KeywordValidator {
            Opm::FlowKeywordValidation::unsupportedKeywords(),
            Opm::FlowKeywordValidation::partiallySupported<std::string>(),
            Opm::FlowKeywordValidation::partiallySupported<int>(),
            Opm::FlowKeywordValidation::partiallySupported<double>(),
            Opm::KeywordValidation::specialValidation()
        };

        keyword_validator.validateDeck(deck, parseContext, treatCriticalAsNonCritical, errorGuard);

        if (checkDeck) {
            Opm::checkDeck(deck, parser, parseContext, errorGuard);
        }

        return deck;
    }

    std::shared_ptr<Opm::EclipseState>
    createEclipseState([[maybe_unused]] Opm::Parallel::Communication comm,
                       const Opm::Deck&                              deck)
    {
#if HAVE_MPI
        return std::make_shared<Opm::ParallelEclipseState>(deck, comm);
#else
        return std::make_shared<Opm::EclipseState>(deck);
#endif
    }

    void readOnIORank(Opm::Parallel::Communication         comm,
                      const std::string&                   deckFilename,
                      const Opm::ParseContext*             parseContext,
                      std::shared_ptr<Opm::EclipseState>&  eclipseState,
                      std::shared_ptr<Opm::Schedule>&      schedule,
                      std::unique_ptr<Opm::UDQState>&      udqState,
                      std::unique_ptr<Opm::Action::State>& actionState,
                      std::unique_ptr<Opm::WellTestState>& wtestState,
                      std::shared_ptr<Opm::SummaryConfig>& summaryConfig,
                      std::shared_ptr<Opm::Python>         python,
                      const bool                           initFromRestart,
                      const bool                           checkDeck,
                      const bool                           treatCriticalAsNonCritical,
                      const std::optional<int>&            outputInterval,
                      Opm::ErrorGuard&                     errorGuard)
    {
        OPM_TIMEBLOCK(readDeck);
        if (((schedule == nullptr) || (summaryConfig == nullptr)) &&
            (parseContext == nullptr))
        {
            OPM_THROW(std::logic_error,
                      "We need a parse context if schedule "
                      "or summaryConfig are not initialized");
        }

        auto parser = Opm::Parser{};
        const auto deck = readDeckFile(deckFilename, checkDeck, parser,
                                       *parseContext, treatCriticalAsNonCritical, errorGuard);

        if (eclipseState == nullptr) {
            OPM_TIMEBLOCK(createEclState);
            eclipseState = createEclipseState(comm, deck);
        }

        if (eclipseState->getInitConfig().restartRequested()) {
            loadObjectsFromRestart(deck, parser, *parseContext,
                                   initFromRestart, outputInterval,
                                   *eclipseState, std::move(python),
                                   schedule, udqState, actionState, wtestState,
                                   errorGuard);
        }
        else {
            createNonRestartDynamicObjects(deck, *eclipseState,
                                           *parseContext, std::move(python),
                                           schedule, udqState, actionState, wtestState,
                                           errorGuard);
        }

        eclipseState->appendAqufluxSchedule(schedule->getAquiferFluxSchedule());

        if (Opm::OpmLog::hasBackend("STDOUT_LOGGER")) {
            // loggers might not be set up!
            setupMessageLimiter((*schedule)[0].message_limits(), "STDOUT_LOGGER");
        }

        if (summaryConfig == nullptr) {
            summaryConfig = std::make_shared<Opm::SummaryConfig>
                (deck, *schedule, eclipseState->fieldProps(),
                 eclipseState->aquifer(), *parseContext, errorGuard);
        }

        Opm::checkConsistentArrayDimensions(*eclipseState, *schedule,
                                            *parseContext, errorGuard);
    }

#if HAVE_MPI
    void defineStateObjectsOnNonIORank(Opm::Parallel::Communication         comm,
                                       std::shared_ptr<Opm::Python>         python,
                                       std::shared_ptr<Opm::EclipseState>&  eclipseState,
                                       std::shared_ptr<Opm::Schedule>&      schedule,
                                       std::unique_ptr<Opm::UDQState>&      udqState,
                                       std::unique_ptr<Opm::Action::State>& actionState,
                                       std::unique_ptr<Opm::WellTestState>& wtestState,
                                       std::shared_ptr<Opm::SummaryConfig>& summaryConfig)
    {
        if (eclipseState == nullptr) {
            eclipseState = std::make_shared<Opm::ParallelEclipseState>(comm);
        }

        if (schedule == nullptr) {
            schedule = std::make_shared<Opm::Schedule>(std::move(python));
        }

        if (udqState == nullptr) {
            udqState = std::make_unique<Opm::UDQState>(0);
        }

        if (actionState == nullptr) {
            actionState = std::make_unique<Opm::Action::State>();
        }

        if (wtestState == nullptr) {
            wtestState = std::make_unique<Opm::WellTestState>();
        }

        if (summaryConfig == nullptr) {
            summaryConfig = std::make_shared<Opm::SummaryConfig>();
        }
    }
#endif

    std::pair<bool, std::string>
    gridHasValidCellGeometry(const Opm::EclipseGrid& inputGrid,
                             const Opm::UnitSystem&  usys)
    {
        const auto numActive = inputGrid.getNumActive();

        if (numActive == 0)
        {
            return {false, R"(Input grid has no active cells at all.
Please check geometry keywords and modifications to e.g. pore volume,
especially if grid is imported through GDFILE.)"};
        }

        for (auto activeCell = 0*numActive; activeCell < numActive; ++activeCell) {
            if (inputGrid.isValidCellGeomtry(inputGrid.getGlobalIndex(activeCell), usys)) {
                return {true, ""};
            }
        }

        return {false, R"(No active cell in input grid has valid/finite cell geometry
Please check geometry keywords and modifications to e.g. pore volume,
especially if the grid is imported through GDFILE.)"};
    }

    std::pair<bool, std::string>
    gridHasValidCellGeometry(Opm::Parallel::Communication comm,
                                  const Opm::EclipseState&     eclipseState)
    {
        bool hasValidCells = false;
        std::string message;

        if (comm.rank() == 0) {
            std::tie(hasValidCells, message) =
                gridHasValidCellGeometry(eclipseState.getInputGrid(),
                                         eclipseState.getDeckUnitSystem());
        }

#if HAVE_MPI
        const auto status = comm.broadcast(&hasValidCells, 1, 0);

        if (status != MPI_SUCCESS) {
            throw std::invalid_argument {
                "Unable to establish cell geometry validity across MPI ranks"
            };
        }
#endif  // HAVE_MPI

        return {hasValidCells, message};
    }
}

// ---------------------------------------------------------------------------


void Opm::ensureOutputDirExists(const std::string& cmdline_output_dir)
{
    namespace fs = std::filesystem;

    if (! fs::is_directory(cmdline_output_dir)) {
        try {
            fs::create_directories(cmdline_output_dir);
        }
        catch (...) {
            throw std::runtime_error {
                fmt::format("Creation of output directory '{}' failed",
                            cmdline_output_dir)
                    };
        }
    }
}

// Setup the OpmLog backends
Opm::FileOutputMode
Opm::setupLogging(const int          mpi_rank_,
                  const std::string& deck_filename,
                  const std::string& cmdline_output_dir,
                  const std::string& cmdline_output,
                  const bool         output_cout_,
                  const std::string& stdout_log_id,
                  const bool         allRanksDbgLog)
{
    if (!cmdline_output_dir.empty()) {
        ensureOutputDirExists(cmdline_output_dir);
    }

    // create logFile
    using std::filesystem::path;
    path fpath(deck_filename);
    std::string baseName;
    std::ostringstream debugFileStream;
    std::ostringstream logFileStream;

    // Strip extension "." or ".DATA"
    std::string extension = uppercase(fpath.extension().string());
    if (extension == ".DATA" || extension == ".") {
        baseName = uppercase(fpath.stem().string());
    }
    else {
        baseName = uppercase(fpath.filename().string());
    }

    std::string output_dir = cmdline_output_dir;
    if (output_dir.empty()) {
        output_dir = fpath.has_parent_path()
            ? absolute(fpath.parent_path()).generic_string()
            : std::filesystem::current_path().generic_string();
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
            std::cerr << "Value " << cmdline_output
                      << " is not a recognized output mode. Using \"all\" instead.\n";
        }
        if (!allRanksDbgLog && mpi_rank_ != 0)
        {
            output = FileOutputMode::OUTPUT_NONE;
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
        // Set a tag limit of 10 (no category limit). Will later in
        // the run be replaced by calling setupMessageLimiter(), after
        // the deck is read and the (possibly user-set) category
        // limits are known.
        streamLog->setMessageLimiter(std::make_shared<Opm::MessageLimiter>(10));
        bool use_color_coding = OpmLog::stdoutIsTerminal();
        streamLog->setMessageFormatter(std::make_shared<Opm::SimpleMessageFormatter>(use_color_coding));
    }

    return output;
}

void Opm::readDeck(Opm::Parallel::Communication    comm,
                   const std::string&              deckFilename,
                   std::shared_ptr<EclipseState>&  eclipseState,
                   std::shared_ptr<Schedule>&      schedule,
                   std::unique_ptr<UDQState>&      udqState,
                   std::unique_ptr<Action::State>& actionState,
                   std::unique_ptr<WellTestState>& wtestState,
                   std::shared_ptr<SummaryConfig>& summaryConfig,
                   std::shared_ptr<Python>         python,
                   const std::string&              parsingStrictness,
                   const bool                      initFromRestart,
                   const bool                      checkDeck,
                   const std::optional<int>&       outputInterval)
{
    auto errorGuard = std::make_unique<ErrorGuard>();

    int parseSuccess = 1; // > 0 is success
    std::string failureMessage;

    if (parsingStrictness != "high" && parsingStrictness != "normal" && parsingStrictness != "low") {
        OPM_THROW(std::runtime_error,
                  fmt::format("Incorrect value {} for parameter ParsingStrictness, must be 'high', 'normal', or 'low'", parsingStrictness));
    }

    if (comm.rank() == 0) { // Always true when !HAVE_MPI
        const bool exitOnAllErrors = (parsingStrictness == "high");
        const bool treatCriticalAsNonCritical = (parsingStrictness == "low");
        try {
            auto parseContext = setupParseContext(exitOnAllErrors);
            if (treatCriticalAsNonCritical) { // Continue with invalid names if parsing strictness is set to low
                parseContext->update(ParseContext::SCHEDULE_INVALID_NAME, InputErrorAction::WARN);
            }
            readOnIORank(comm, deckFilename, parseContext.get(),
                         eclipseState, schedule, udqState, actionState, wtestState,
                         summaryConfig, std::move(python), initFromRestart,
                         checkDeck, treatCriticalAsNonCritical, outputInterval, *errorGuard);
        }
        catch (const OpmInputError& input_error) {
            failureMessage = input_error.what();
            parseSuccess = 0;
        }
        catch (const std::exception& std_error) {
            failureMessage = std_error.what();
            parseSuccess = 0;
        }
    }

#if HAVE_MPI
    else {
        defineStateObjectsOnNonIORank(comm, std::move(python), eclipseState,
                                      schedule, udqState, actionState, wtestState,
                                      summaryConfig);
    }

    // In case of parse errors eclipseState/schedule might be null and
    // trigger segmentation faults in parallel during broadcast (e.g. when
    // serializing the non-existent TableManager)
    parseSuccess = comm.min(parseSuccess);
    try {
        if (parseSuccess) {
            OPM_TIMEBLOCK(eclBcast);
            eclStateBroadcast(comm, *eclipseState, *schedule,
                              *summaryConfig, *udqState, *actionState, *wtestState);
        }
    }
    catch (const std::exception& broadcast_error) {
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

    parseSuccess = comm.min(parseSuccess);

    if (! parseSuccess) {
        if (comm.rank() == 0) {
            OpmLog::error(fmt::format("Unrecoverable errors while loading input: {}", failureMessage));
        }

#if HAVE_MPI
        MPI_Finalize();
#endif

        std::exit(EXIT_FAILURE);
    }
}

void Opm::verifyValidCellGeometry(Parallel::Communication comm,
                                  const EclipseState&     eclipseState)
{
    auto [hasValidCells, message] = gridHasValidCellGeometry(comm, eclipseState);
    if (hasValidCells) {
        return;
    }

    throw std::invalid_argument { message };
}

std::unique_ptr<Opm::ParseContext> Opm::setupParseContext(const bool strictParsing)
{
    auto parseContext =
        std::make_unique<ParseContext>(std::vector<std::pair<std::string , InputErrorAction>>
                                       {{ParseContext::PARSE_RANDOM_SLASH, InputErrorAction::IGNORE},
                                       {ParseContext::PARSE_MISSING_DIMS_KEYWORD, InputErrorAction::WARN},
                                       {ParseContext::SUMMARY_UNKNOWN_WELL, InputErrorAction::WARN},
                                       {ParseContext::SUMMARY_UNKNOWN_GROUP, InputErrorAction::WARN}});
    if (strictParsing)
        parseContext->update(InputErrorAction::DELAYED_EXIT1);

    return parseContext;
}
