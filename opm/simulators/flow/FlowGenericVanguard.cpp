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
#include <opm/simulators/flow/FlowGenericVanguard.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <opm/common/utility/MemPacker.hpp>
#include <opm/common/utility/Serializer.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquiferCell.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQParams.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/input/eclipse/Parser/InputErrorAction.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/utils/readDeck.hpp>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <filesystem>
#include <stdexcept>

namespace Opm {

std::unique_ptr<Parallel::Communication> FlowGenericVanguard::comm_;
FlowGenericVanguard::SimulationModelParams FlowGenericVanguard::modelParams_;

FlowGenericVanguard::FlowGenericVanguard()
    : FlowGenericVanguard(std::move(modelParams_))
{}

FlowGenericVanguard::FlowGenericVanguard(SimulationModelParams&& params)
    : python(std::make_shared<Python>())
{
    defineSimulationModel(std::move(params));
}

FlowGenericVanguard::SimulationModelParams
FlowGenericVanguard::serializationTestParams()
{
    SimulationModelParams result;
    result.setupTime_ = 1.234;
    result.actionState_ = std::make_unique<Action::State>(Action::State::serializationTestObject());
    result.eclSchedule_ = std::make_unique<Schedule>(Schedule::serializationTestObject());
    result.summaryState_ = std::make_unique<SummaryState>(SummaryState::serializationTestObject());
    result.udqState_ = std::make_unique<UDQState>(UDQState::serializationTestObject());
    // Remaining members left as null pointers: wtestState_, eclState_ and eclSummaryConfig_.

    return result;
}

FlowGenericVanguard::~FlowGenericVanguard() = default;

void FlowGenericVanguard::defineSimulationModel(SimulationModelParams&& params)
{
    actionState_ = std::move(params.actionState_);
    eclSchedule_ = std::move(params.eclSchedule_);
    eclState_ = std::move(params.eclState_);
    eclSummaryConfig_ = std::move(params.eclSummaryConfig_);
    setupTime_ = params.setupTime_;
    udqState_ = std::move(params.udqState_);
    wtestState_ = std::move(params.wtestState_);
    summaryState_ = std::move(params.summaryState_);
}

void FlowGenericVanguard::readDeck(const std::string& filename)
{
    Dune::Timer setupTimer;
    setupTimer.start();

    Opm::readDeck(comm(),
                  filename,
                  modelParams_.eclState_,
                  modelParams_.eclSchedule_,
                  modelParams_.udqState_,
                  modelParams_.actionState_,
                  modelParams_.wtestState_,
                  modelParams_.eclSummaryConfig_,
                  nullptr, "normal", "normal", "100", false, false, false, {}, /*slaveMode=*/false);
    modelParams_.setupTime_ = setupTimer.stop();
}

std::string FlowGenericVanguard::canonicalDeckPath(const std::string& caseName)
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

void FlowGenericVanguard::updateNOSIM_(std::string_view dryRunString)
{
    try {
        // Possible to force initialization only behavior (NOSIM).
        if (dryRunString != "" && dryRunString != "auto") {
            bool enableDryRun;
            if (dryRunString == "true"
                || dryRunString == "t"
                || dryRunString == "1")
                enableDryRun = true;
            else if (dryRunString == "false"
                        || dryRunString == "f"
                        || dryRunString == "0")
                enableDryRun = false;
            else
                throw std::invalid_argument("Invalid value for parameter EnableDryRun: '"
                                            + std::string(dryRunString) + "'");
            auto& ioConfig = eclState().getIOConfig();
            ioConfig.overrideNOSIM(enableDryRun);
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Failed to create valid EclipseState object" << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
        throw;
    }
}

void FlowGenericVanguard::updateOutputDir_(std::string outputDir,
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

void FlowGenericVanguard::init()
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

    // set communicator if not set as in opm flow
    if(!comm_){
        FlowGenericVanguard::setCommunication(std::make_unique<Parallel::Communication>());
    }
    
    // set eclState if not already set as in opm flow
    // it means that setParams is called
    if(!eclState_){
        this->readDeck(fileName_);
        this->defineSimulationModel(std::move(this->modelParams_));
    }

    
    if (!this->summaryState_) {
        this->summaryState_ = std::make_unique<SummaryState>
            (TimeService::from_time_t(this->eclSchedule_->getStartTime()),
             this->eclState_->runspec().udqParams().undefinedValue());
    }

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
        const auto& comm = FlowGenericVanguard::comm();

        if (useMultisegmentWell_)
        {
            if (comm.rank() == 0)
            {
                const auto& wells = this->schedule().getWellsatEnd();
                for (const auto& well : wells)
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

bool FlowGenericVanguard::drsdtconEnabled() const
{
  for (const auto& schIt : this->schedule()) {
      const auto& oilVaporizationControl = schIt.oilvap();
      if (oilVaporizationControl.getType() == OilVaporizationProperties::OilVaporization::DRSDTCON) {
          return true;
      }
  }

  return false;
}

std::unordered_map<size_t, const NumericalAquiferCell*>
FlowGenericVanguard::allAquiferCells() const
{
  return this->eclState_->aquifer().numericalAquifers().allAquiferCells();
}

template<>
void FlowGenericVanguard::
serializeOp<Serializer<Serialization::MemPacker>>(Serializer<Serialization::MemPacker>& serializer)
{
    serializer(*summaryState_);
    serializer(*udqState_);
    serializer(*actionState_);
    serializer(*eclSchedule_);
}

bool FlowGenericVanguard::operator==(const FlowGenericVanguard& rhs) const
{
    auto cmp_ptr = [](const auto& a, const auto& b)
    {
        if (!a && !b) {
            return true;
        }

        if (a && b) {
            return *a == *b;
        }

        return false;
    };
    return cmp_ptr(this->summaryState_, rhs.summaryState_);
           cmp_ptr(this->udqState_, rhs.udqState_) &&
           cmp_ptr(this->actionState_, rhs.actionState_) &&
           cmp_ptr(this->eclSchedule_, rhs.eclSchedule_);
}

template<class Scalar>
void FlowGenericVanguard::registerParameters_()
{
    Parameters::Register<Parameters::EclDeckFileName>
        ("The name of the file which contains the ECL deck to be simulated");
    Parameters::Register<Parameters::EclOutputInterval>
        ("The number of report steps that ought to be skipped between two writes of ECL results");
    Parameters::Register<Parameters::EnableDryRun>
        ("Specify if the simulation ought to be actually run, or just pretended to be");
    Parameters::Register<Parameters::EnableOpmRstFile>
        ("Include OPM-specific keywords in the ECL restart file to "
         "enable restart of OPM simulators from these files");
    Parameters::Register<Parameters::IgnoreKeywords>
        ("List of Eclipse keywords which should be ignored. As a ':' separated string.");
    Parameters::Register<Parameters::ParsingStrictness>
        ("Set strictness of parsing process. Available options are "
         "normal (stop for critical errors), "
         "high (stop for all errors) and "
         "low (as normal, except do not stop due to unsupported "
         "keywords even if marked critical");
    Parameters::Register<Parameters::ActionParsingStrictness>
        ("Set strictness of parsing process for ActionX and PyAction. Available options are "
         "normal (do not apply keywords that have not been tested for ActionX or PyAction) and "
         "low (try to apply all keywords, beware: the simulation outcome might be incorrect).");
    Parameters::Register<Parameters::InputSkipMode>
        ("Set compatibility mode for the SKIP100/SKIP300 keywords. Options are "
         "100 (skip SKIP100..ENDSKIP, keep SKIP300..ENDSKIP) [default], "
         "300 (skip SKIP300..ENDSKIP, keep SKIP100..ENDSKIP) and "
         "all (skip both SKIP100..ENDSKIP and SKIP300..ENDSKIP) ");
    Parameters::Register<Parameters::SchedRestart>
        ("When restarting: should we try to initialize wells and "
         "groups from historical SCHEDULE section.");
    Parameters::Register<Parameters::EdgeWeightsMethod>
        ("Choose edge-weighing strategy: 0=uniform, 1=trans, 2=log(trans).");

#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
    Parameters::Register<Parameters::NumJacobiBlocks>
        ("Number of blocks to be created for the Block-Jacobi preconditioner.");
#endif

    Parameters::Register<Parameters::OwnerCellsFirst>
        ("Order cells owned by rank before ghost/overlap cells.");
#if HAVE_MPI
    Parameters::Register<Parameters::PartitionMethod>
        ("Choose partitioning strategy: 0=simple, 1=Zoltan, 2=METIS.");
    Parameters::Register<Parameters::SerialPartitioning>
        ("Perform partitioning for parallel runs on a single process.");
    Parameters::Register<Parameters::ZoltanImbalanceTol<Scalar>>
        ("Tolerable imbalance of the loadbalancing provided by Zoltan. DEPRECATED: Use --imbalance-tol instead");
    Parameters::Register<Parameters::ZoltanParams>
        ("Configuration of Zoltan partitioner. "
         "Valid options are: graph, hypergraph or scotch. "
         "Alternatively, you can request a configuration to be read "
         "from a JSON file by giving the filename here, ending with '.json.' "
         "See https://sandialabs.github.io/Zoltan/ug_html/ug.html "
         "for available Zoltan options.");
    Parameters::Register<Parameters::ImbalanceTol<Scalar>>
        ("Tolerable imbalance of the loadbalancing (default: 1.1).");
    Parameters::Register<Parameters::MetisParams>
        ("Configuration of Metis partitioner. "
         "You can request a configuration to be read "
         "from a JSON file by giving the filename here, ending with '.json.' "
         "See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf"
         "for available METIS options.");
    Parameters::Register<Parameters::ExternalPartition>
        ("Name of file from which to load an externally generated "
         "partitioning of the model's active cells for MPI "
         "distribution purposes. If empty, the built-in partitioning "
         "method will be employed.");
    Parameters::Hide<Parameters::ExternalPartition>();

    Parameters::Hide<Parameters::ZoltanImbalanceTol<Scalar>>();
    Parameters::Hide<Parameters::ZoltanParams>();
#endif
    Parameters::Register<Parameters::AllowDistributedWells>
        ("Allow the perforations of a well to be distributed to interior of multiple processes");
    // register here for the use in the tests without BlackoilModelParameters
    Parameters::Register<Parameters::UseMultisegmentWell>
        ("Use the well model for multi-segment wells instead of the one for single-segment wells");
}

template void FlowGenericVanguard::registerParameters_<double>();

#if FLOW_INSTANTIATE_FLOAT
template void FlowGenericVanguard::registerParameters_<float>();
#endif

} // namespace Opm
