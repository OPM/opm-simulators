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
#include <opm/common/utility/TimeService.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquiferCell.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <filesystem>
#include <stdexcept>

namespace Opm {

double EclGenericVanguard::setupTime_ = 0.0;
std::shared_ptr<EclipseState> EclGenericVanguard::eclState_;
std::shared_ptr<Schedule> EclGenericVanguard::eclSchedule_;
std::shared_ptr<SummaryConfig> EclGenericVanguard::eclSummaryConfig_;
std::unique_ptr<UDQState> EclGenericVanguard::udqState_;
std::unique_ptr<Action::State> EclGenericVanguard::actionState_;
std::unique_ptr<WellTestState> EclGenericVanguard::wtestState_;
std::unique_ptr<Parallel::Communication> EclGenericVanguard::comm_;

EclGenericVanguard::EclGenericVanguard()
    : python(std::make_shared<Python>())
{
}

EclGenericVanguard::~EclGenericVanguard() = default;

void EclGenericVanguard::setParams(double setupTime,
                                   std::shared_ptr<EclipseState> eclState,
                                   std::shared_ptr<Schedule> schedule,
                                   std::unique_ptr<UDQState> udqState,
                                   std::unique_ptr<Action::State> actionState,
                                   std::unique_ptr<WellTestState> wtestState,
                                   std::shared_ptr<SummaryConfig> summaryConfig)
{
    EclGenericVanguard::setupTime_ = setupTime;
    EclGenericVanguard::eclState_ = std::move(eclState);
    EclGenericVanguard::eclSchedule_ = std::move(schedule);
    EclGenericVanguard::udqState_ = std::move(udqState);
    EclGenericVanguard::actionState_ = std::move(actionState);
    EclGenericVanguard::wtestState_ = std::move(wtestState);
    EclGenericVanguard::eclSummaryConfig_ = std::move(summaryConfig);
}

void EclGenericVanguard::readDeck(const std::string& filename)
{
    Dune::Timer setupTimer;
    setupTimer.start();

    std::shared_ptr<Opm::EclipseState> eclipseState;
    std::shared_ptr<Opm::Schedule> schedule;
    std::unique_ptr<Opm::UDQState> udqState;
    std::unique_ptr<Opm::Action::State> actionState;
    std::unique_ptr<Opm::WellTestState> wtestState;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig;

    auto parseContext =
        std::make_unique<ParseContext>(std::vector<std::pair<std::string , InputError::Action>>
                                            {{ParseContext::PARSE_RANDOM_SLASH, InputError::IGNORE},
                                             {ParseContext::PARSE_MISSING_DIMS_KEYWORD, InputError::WARN},
                                             {ParseContext::SUMMARY_UNKNOWN_WELL, InputError::WARN},
                                             {ParseContext::SUMMARY_UNKNOWN_GROUP, InputError::WARN}});

    Opm::readDeck(EclGenericVanguard::comm(),
                  filename, eclipseState, schedule, udqState,
                  actionState, wtestState,
                  summaryConfig, nullptr, nullptr, std::move(parseContext),
                  false, false, {});

    EclGenericVanguard::setParams(setupTimer.elapsed(),
                                  eclipseState, schedule,
                                  std::move(udqState),
                                  std::move(actionState),
                                  std::move(wtestState), summaryConfig);
}

std::string EclGenericVanguard::canonicalDeckPath(const std::string& caseName)
{
    const auto fileExists = [](const std::filesystem::path& f) -> bool
    {
        if (!std::filesystem::exists(f))
            return false;

        if (std::filesystem::is_regular_file(f))
            return true;

        return std::filesystem::is_symlink(f) && std::filesystem::is_regular_file(std::filesystem::read_symlink(f));
    };

    auto simcase = std::filesystem::path(caseName);
    if (fileExists(simcase))
        return simcase.string();

    for (const auto& ext : { std::string("data"), std::string("DATA") }) {
        if (fileExists(simcase.replace_extension(ext)))
            return simcase.string();
    }

    throw std::invalid_argument("Cannot find input case '"+caseName+"'");
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
    if (!std::filesystem::is_directory(outputDir)) {
        try {
            std::filesystem::create_directories(outputDir);
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

    this->summaryState_ = std::make_unique<SummaryState>( TimeService::from_time_t(this->eclSchedule_->getStartTime() ));

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
        const auto& comm = EclGenericVanguard::comm();

        if (useMultisegmentWell_)
        {
            if (comm.rank() == 0)
            {
                const auto& wells = this->schedule().getWellsatEnd();
                for ( const auto& well: wells)
                {
                    hasMsWell = hasMsWell || well.isMultiSegment();
                }
            }
        }

        hasMsWell = comm.max(hasMsWell);

        if (hasMsWell)
        {
            if (comm.rank() == 0)
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
