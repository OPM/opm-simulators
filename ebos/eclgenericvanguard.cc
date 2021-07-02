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

#include <config.h>
#include <ebos/eclgenericvanguard.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/FileSystem.hpp>
#include <opm/common/utility/TimeService.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquiferCell.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/OilVaporizationProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Python/Python.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <stdexcept>

namespace Opm {

double EclGenericVanguard::externalSetupTime_ = 0.0;
std::unique_ptr<ParseContext> EclGenericVanguard::externalParseContext_;
std::unique_ptr<ErrorGuard> EclGenericVanguard::externalErrorGuard_;
std::unique_ptr<Deck> EclGenericVanguard::externalDeck_;
bool EclGenericVanguard::externalDeckSet_ = false;
std::unique_ptr<EclipseState> EclGenericVanguard::externalEclState_;
std::unique_ptr<Schedule> EclGenericVanguard::externalEclSchedule_;
std::unique_ptr<SummaryConfig> EclGenericVanguard::externalEclSummaryConfig_;
std::unique_ptr<UDQState> EclGenericVanguard::externalUDQState_;

EclGenericVanguard::EclGenericVanguard()
    : python(std::make_shared<Python>())
{
}

EclGenericVanguard::~EclGenericVanguard() = default;

void EclGenericVanguard::setExternalParseContext(std::unique_ptr<ParseContext> parseContext)
{
    externalParseContext_ = std::move(parseContext);
}

void EclGenericVanguard::setExternalErrorGuard(std::unique_ptr<ErrorGuard> errorGuard)
{
    externalErrorGuard_ = std::move(errorGuard);
}

void EclGenericVanguard::setExternalSchedule(std::unique_ptr<Schedule> schedule)
{
    externalEclSchedule_ = std::move(schedule);
}

void EclGenericVanguard::setExternalSummaryConfig(std::unique_ptr<SummaryConfig> summaryConfig)
{
    externalEclSummaryConfig_ = std::move(summaryConfig);
}

void EclGenericVanguard::setExternalDeck(std::unique_ptr<Deck> deck)
{
    externalDeck_ = std::move(deck);
    externalDeckSet_ = true;
}

void EclGenericVanguard::setExternalEclState(std::unique_ptr<EclipseState> eclState)
{
    externalEclState_ = std::move(eclState);
}

void EclGenericVanguard::setExternalUDQState(std::unique_ptr<UDQState> udqState)
{
    externalUDQState_ = std::move(udqState);
}

std::string EclGenericVanguard::canonicalDeckPath(const std::string& caseName)
{
    const auto fileExists = [](const filesystem::path& f) -> bool
    {
        if (!filesystem::exists(f))
            return false;

        if (filesystem::is_regular_file(f))
            return true;

        return filesystem::is_symlink(f) && filesystem::is_regular_file(filesystem::read_symlink(f));
    };

    auto simcase = filesystem::path(caseName);
    if (fileExists(simcase))
        return simcase.string();

    for (const auto& ext : { std::string("data"), std::string("DATA") }) {
        if (fileExists(simcase.replace_extension(ext)))
            return simcase.string();
    }

    throw std::invalid_argument("Cannot find input case '"+caseName+"'");
}

std::unique_ptr<ParseContext> EclGenericVanguard::createParseContext(const std::string& ignoredKeywords,
                                                                     bool eclStrictParsing)
{
    typedef std::pair<std::string, InputError::Action> ParseModePair;
    typedef std::vector<ParseModePair> ParseModePairs;
    ParseModePairs tmp;
    tmp.emplace_back(ParseContext::PARSE_RANDOM_SLASH, InputError::IGNORE);
    tmp.emplace_back(ParseContext::PARSE_MISSING_DIMS_KEYWORD, InputError::WARN);
    tmp.emplace_back(ParseContext::SUMMARY_UNKNOWN_WELL, InputError::WARN);
    tmp.emplace_back(ParseContext::SUMMARY_UNKNOWN_GROUP, InputError::WARN);
    tmp.emplace_back(ParseContext::PARSE_EXTRA_RECORDS, InputError::WARN);

    auto parseContext = std::make_unique<ParseContext>(tmp);

    if (!ignoredKeywords.empty()) {
        size_t pos;
        size_t offset = 0;
        while (true) {
            pos = ignoredKeywords.find(':', offset);
            if (pos == std::string::npos) {
                parseContext->ignoreKeyword(ignoredKeywords.substr(offset));
                break;
            }
            parseContext->ignoreKeyword(ignoredKeywords.substr(offset, pos - offset));
            offset = pos + 1;
        }
    }

    if (eclStrictParsing)
        parseContext->update(InputError::DELAYED_EXIT1);

    return parseContext;
}

void EclGenericVanguard::updateOutputDir_(std::string outputDir,
                                          bool enableEclCompatFile)
{
    // update the location for output
    auto& ioConfig = eclState_->getIOConfig();
    if (outputDir.empty())
        // If no output directory parameter is specified, use the output directory
        // which Opm::IOConfig thinks that should be used. Normally this is the
        // directory in which the input files are located.
        outputDir = ioConfig.getOutputDir();

    // ensure that the output directory exists and that it is a directory
    if (!filesystem::is_directory(outputDir)) {
        try {
            filesystem::create_directories(outputDir);
        }
        catch (...) {
            throw std::runtime_error("Creation of output directory '"+outputDir+"' failed\n");
        }
    }

    // specify the directory output. This is not a very nice mechanism because
    // the eclState is supposed to be immutable here, IMO.
    ioConfig.setOutputDir(outputDir);

    ioConfig.setEclCompatibleRST(enableEclCompatFile);
}

void EclGenericVanguard::init()
{
    int myRank = 0;

#if HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

    // Make proper case name.
    {
        if (fileName_.empty())
            throw std::runtime_error("No input deck file has been specified as a command line argument,"
                                        " or via '--ecl-deck-file-name=CASE.DATA'");

        fileName_ = canonicalDeckPath(fileName_);

        // compute the base name of the input file name
        const char directorySeparator = '/';
        long int i;
        for (i = fileName_.size(); i >= 0; -- i)
            if (fileName_[i] == directorySeparator)
                break;
        std::string baseName = fileName_.substr(i + 1, fileName_.size());

        // remove the extension from the input file
        for (i = baseName.size(); i >= 0; -- i)
            if (baseName[i] == '.')
                break;
        std::string rawCaseName;
        if (i < 0)
            rawCaseName = baseName;
        else
            rawCaseName = baseName.substr(0, i);

        // transform the result to ALL_UPPERCASE
        caseName_ = rawCaseName;
        std::transform(caseName_.begin(), caseName_.end(), caseName_.begin(), ::toupper);
    }

    std::unique_ptr<ErrorGuard> errorGuard = nullptr;

    // Check that we are in one of the known configurations for external variables
    // and move them to internal
    if (externalDeck_)
    {
        deck_ = std::move(externalDeck_);

        if (externalParseContext_ && externalErrorGuard_ )
        {
            parseContext_ = std::move(externalParseContext_);
            errorGuard = std::move(externalErrorGuard_);
        }
        else if(externalEclState_ && externalEclSchedule_ && externalEclSummaryConfig_)
        {
            eclState_ = std::move(externalEclState_);
            eclSchedule_ = std::move(externalEclSchedule_);
            eclSummaryConfig_ = std::move(externalEclSummaryConfig_);
        }
        else
        {
            OPM_THROW(std::logic_error, "Either parse context and error guard or ECL state, schedule, and summary config need to be"
                          << " set externally.");
        }
    }
    else if (externalParseContext_)
    {
        parseContext_ = std::move(externalParseContext_);
    }
    else
    {
        parseContext_ = createParseContext(ignoredKeywords_, eclStrictParsing_);
    }

    readDeck(myRank, fileName_, deck_, eclState_, eclSchedule_, udqState_,
             eclSummaryConfig_, std::move(errorGuard), python,
             std::move(parseContext_), /* initFromRestart = */ false,
             /* checkDeck = */ enableExperiments_, outputInterval_);

    if (EclGenericVanguard::externalUDQState_)
        this->udqState_ = std::move(EclGenericVanguard::externalUDQState_);
    else
        this->udqState_ = std::make_unique<UDQState>( this->eclSchedule_->getUDQConfig(0).params().undefinedValue() );
    this->summaryState_ = std::make_unique<SummaryState>( TimeService::from_time_t(this->eclSchedule_->getStartTime() ));
    this->actionState_ = std::make_unique<Action::State>();

    // Initialize parallelWells with all local wells
    const auto& schedule_wells = schedule().getWellsatEnd();
    parallelWells_.reserve(schedule_wells.size());

    for (const auto& well: schedule_wells)
    {
        parallelWells_.emplace_back(well.name(), true);
    }
    std::sort(parallelWells_.begin(), parallelWells_.end());

    // Check whether allowing distribute wells makes sense
    if (enableDistributedWells() )
    {
        int hasMsWell = false;

        if (useMultisegmentWell_)
        {
            if (myRank == 0)
            {
                const auto& wells = this->schedule().getWellsatEnd();
                for ( const auto& well: wells)
                {
                    hasMsWell = hasMsWell || well.isMultiSegment();
                }
            }
        }
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
        const auto& comm = Dune::MPIHelper::getCommunication();
#else
        const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
#endif
        hasMsWell = comm.max(hasMsWell);

        if (hasMsWell)
        {
            if (myRank == 0)
            {
                std::string message =
                        std::string("Option --allow-distributed-wells=true is only allowed if model\n")
                        + "only has only standard wells. You need to provide option \n"
                        + " with --enable-multisegement-wells=false to treat existing \n"
                        + "multisegment wells as standard wells.";
                OpmLog::error(message);
            }
            comm.barrier();
            OPM_THROW(std::invalid_argument, "All wells need to be standard wells!");
        }
    }
}

bool EclGenericVanguard::drsdtconEnabled() const
{
  for (const auto& schIt : this->schedule()) {
      const auto& oilVaporizationControl = schIt.oilvap();
      if (oilVaporizationControl.getType() == OilVaporizationProperties::OilVaporization::DRSDTCON) {
          return true;
      }
  }

  return false;
}

std::unordered_map<size_t, const NumericalAquiferCell*> EclGenericVanguard::allAquiferCells() const
{
  return this->eclState_->aquifer().numericalAquifers().allAquiferCells();
}

} // namespace Opm
