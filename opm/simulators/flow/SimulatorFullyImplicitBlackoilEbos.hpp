/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILEBOS_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILEBOS_HEADER_INCLUDED

#include <dune/common/hash.hh>

#include <opm/simulators/flow/BlackoilModelEbos.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>
#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>
#include <opm/simulators/flow/NonlinearSolverEbos.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#if HAVE_HDF5
#include <ebos/hdf5serializer.hh>
#endif

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableAdaptiveTimeStepping {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableTuning {
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct OutputExtraConvergenceInfo
{
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct SaveStep
{
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct LoadStep
{
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct SaveFile
{
    using type = UndefinedProperty;
};

template<class TypeTag>
struct EnableTerminalOutput<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableAdaptiveTimeStepping<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableTuning<TypeTag, TTag::EclFlowProblem> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct OutputExtraConvergenceInfo<TypeTag, TTag::EclFlowProblem>
{
    static constexpr auto* value = "none";
};

template <class TypeTag>
struct SaveStep<TypeTag, TTag::EclFlowProblem>
{
    static constexpr auto* value = "";
};

template <class TypeTag>
struct SaveFile<TypeTag, TTag::EclFlowProblem>
{
    static constexpr auto* value = "";
};

template <class TypeTag>
struct LoadStep<TypeTag, TTag::EclFlowProblem>
{
    static constexpr int value = -1;
};

} // namespace Opm::Properties

namespace Opm {

void outputReportStep(const SimulatorTimer& timer);
void outputTimestampFIP(const SimulatorTimer& timer,
                        const std::string& title,
                        const std::string& version);

/// a simulator for the blackoil model
template<class TypeTag>
class SimulatorFullyImplicitBlackoilEbos
{
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using AquiferModel = GetPropType<TypeTag, Properties::EclAquiferModel>;

    typedef AdaptiveTimeSteppingEbos<TypeTag> TimeStepper;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;
    typedef BlackOilMICPModule<TypeTag> MICPModule;

    typedef BlackoilModelEbos<TypeTag> Model;
    typedef NonlinearSolverEbos<TypeTag, Model> Solver;
    typedef typename Model::ModelParameters ModelParameters;
    typedef typename Solver::SolverParameters SolverParameters;
    typedef BlackoilWellModel<TypeTag> WellModel;


    /// Initialise from parameters and objects to observe.
    /// \param[in] param       parameters, this class accepts the following:
    ///     parameter (default)            effect
    ///     -----------------------------------------------------------
    ///     output (true)                  write output to files?
    ///     output_dir ("output")          output directoty
    ///     output_interval (1)            output every nth step
    ///     nl_pressure_residual_tolerance (0.0) pressure solver residual tolerance (in Pascal)
    ///     nl_pressure_change_tolerance (1.0)   pressure solver change tolerance (in Pascal)
    ///     nl_pressure_maxiter (10)       max nonlinear iterations in pressure
    ///     nl_maxiter (30)                max nonlinear iterations in transport
    ///     nl_tolerance (1e-9)            transport solver absolute residual tolerance
    ///     num_transport_substeps (1)     number of transport steps per pressure step
    ///     use_segregation_split (false)  solve for gravity segregation (if false,
    ///                                    segregation is ignored).
    ///
    /// \param[in] props         fluid and rock properties
    /// \param[in] linsolver     linear solver
    /// \param[in] eclipse_state the object which represents an internalized ECL deck
    /// \param[in] output_writer
    /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
    SimulatorFullyImplicitBlackoilEbos(Simulator& ebosSimulator)
        : ebosSimulator_(ebosSimulator)
    {
        phaseUsage_ = phaseUsageFromDeck(eclState());

        // Only rank 0 does print to std::cout, and only if specifically requested.
        this->terminalOutput_ = false;
        if (this->grid().comm().rank() == 0) {
            this->terminalOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);

            this->startConvergenceOutputThread(EWOMS_GET_PARAM(TypeTag, std::string,
                                                               OutputExtraConvergenceInfo),
                                               R"(OutputExtraConvergenceInfo (--output-extra-convergence-info))");
        }

        const std::string saveSpec = EWOMS_GET_PARAM(TypeTag, std::string, SaveStep);
        if (saveSpec == "all") {
            saveStride_ = 1;
        } else if (!saveSpec.empty() && saveSpec[0] == ':') {
            saveStride_ = std::atoi(saveSpec.c_str()+1);
        } else if (!saveSpec.empty()) {
            saveStep_ = std::atoi(saveSpec.c_str());
        }

        loadStep_ = EWOMS_GET_PARAM(TypeTag, int, LoadStep);

        saveFile_ = EWOMS_GET_PARAM(TypeTag, std::string, SaveFile);
        if (saveFile_.empty()) {
            const auto& ioconfig = ebosSimulator_.vanguard().eclState().getIOConfig();
            saveFile_ = ioconfig.fullBasePath() + ".OPMRST";
            if (loadStep_ != -1 && !std::filesystem::exists(saveFile_)) {
                std::filesystem::path path(ioconfig.getInputDir() + "/");
                path.replace_filename(ioconfig.getBaseName() + ".OPMRST");
                saveFile_ = path;
                if (!std::filesystem::exists(saveFile_)) {
                    OPM_THROW(std::runtime_error, "Error locating serialized restart file");
                }
            }
        }
    }

    ~SimulatorFullyImplicitBlackoilEbos()
    {
        // Safe to call on all ranks, not just the I/O rank.
        this->endConvergenceOutputThread();
    }

    static void registerParameters()
    {
        ModelParameters::registerParameters();
        SolverParameters::registerParameters();
        TimeStepper::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTerminalOutput,
                             "Print high-level information about the simulation's progress to the terminal");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableAdaptiveTimeStepping,
                             "Use adaptive time stepping between report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTuning,
                             "Honor some aspects of the TUNING keyword.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, OutputExtraConvergenceInfo,
                             "Provide additional convergence output "
                             "files for diagnostic purposes. "
                             "\"none\" gives no extra output and "
                             "overrides all other options, "
                             "\"steps\" generates an INFOSTEP file, "
                             "\"iterations\" generates an INFOITER file. "
                             "Combine options with commas, e.g., "
                             "\"steps,iterations\" for multiple outputs.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SaveStep,
                             "Save serialized state to .OPMRST file. "
                             "Either a specific report step, \"all\" to save "
                             "all report steps or \":x\" to save every x'th step.");
        EWOMS_REGISTER_PARAM(TypeTag, int, LoadStep,
                             "Load serialized state from .OPMRST file. "
                             "Either a specific report step, or 0 to load last "
                             "stored report step.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SaveFile,
                             "FileName for .OPMRST file used for serialized state. "
                             "If empty, CASENAME.OPMRST is used.");
        EWOMS_HIDE_PARAM(TypeTag, SaveFile);
    }

    /// Run the simulation.
    /// This will run succesive timesteps until timer.done() is true. It will
    /// modify the reservoir and well states.
    /// \param[in,out] timer       governs the requested reporting timesteps
    /// \param[in,out] state       state of reservoir: pressure, fluxes
    /// \return                    simulation report, with timing data
    SimulatorReport run(SimulatorTimer& timer)
    {
        init(timer);
        // Make cache up to date. No need for updating it in elementCtx.
        ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        // Main simulation loop.
        while (!timer.done()) {
            bool continue_looping = runStep(timer);
            if (!continue_looping) break;
        }
        return finalize();
    }

    void init(SimulatorTimer &timer)
    {
        ebosSimulator_.setEpisodeIndex(-1);

        // Create timers and file for writing timing info.
        solverTimer_ = std::make_unique<time::StopWatch>();
        totalTimer_ = std::make_unique<time::StopWatch>();
        totalTimer_->start();

        // adaptive time stepping
        bool enableAdaptive = EWOMS_GET_PARAM(TypeTag, bool, EnableAdaptiveTimeStepping);
        bool enableTUNING = EWOMS_GET_PARAM(TypeTag, bool, EnableTuning);
        if (enableAdaptive) {
            const UnitSystem& unitSystem = this->ebosSimulator_.vanguard().eclState().getUnits();
            if (enableTUNING) {
                const auto& sched_state = schedule()[timer.currentStepNum()];
                auto max_next_tstep = sched_state.max_next_tstep();
                adaptiveTimeStepping_ = std::make_unique<TimeStepper>(max_next_tstep,
                                                                      sched_state.tuning(),
                                                                      unitSystem, terminalOutput_);
            }
            else {
                adaptiveTimeStepping_ = std::make_unique<TimeStepper>(unitSystem, terminalOutput_);
            }

            if (isRestart()) {
                // For restarts the ebosSimulator may have gotten some information
                // about the next timestep size from the OPMEXTRA field
                adaptiveTimeStepping_->setSuggestedNextStep(ebosSimulator_.timeStepSize());
            }
        }
    }

    bool runStep(SimulatorTimer& timer)
    {
        if (schedule().exitStatus().has_value()) {
            if (terminalOutput_) {
                OpmLog::info("Stopping simulation since EXIT was triggered by an action keyword.");
            }
            report_.success.exit_status = schedule().exitStatus().value();
            return false;
        }

        if (loadStep_ > -1) {
            loadTimerInfo(timer);
        }

        // Report timestep.
        if (terminalOutput_) {
            std::ostringstream ss;
            timer.report(ss);
            OpmLog::debug(ss.str());
        }

        if (terminalOutput_) {
            outputReportStep(timer);
        }

        // write the inital state at the report stage
        if (timer.initialStep()) {
            Dune::Timer perfTimer;
            perfTimer.start();

            ebosSimulator_.setEpisodeIndex(-1);
            ebosSimulator_.setEpisodeLength(0.0);
            ebosSimulator_.setTimeStepSize(0.0);
            wellModel_().beginReportStep(timer.currentStepNum());
            ebosSimulator_.problem().writeOutput();

            report_.success.output_write_time += perfTimer.stop();
        }

        // Run a multiple steps of the solver depending on the time step control.
        solverTimer_->start();

        auto solver = createSolver(wellModel_());

        ebosSimulator_.startNextEpisode(
            ebosSimulator_.startTime()
               + schedule().seconds(timer.currentStepNum()),
            timer.currentStepLength());
        ebosSimulator_.setEpisodeIndex(timer.currentStepNum());
        if (loadStep_> -1) {
            wellModel_().prepareDeserialize(loadStep_ - 1);
            loadSimulatorState();
            loadStep_ = -1;
            ebosSimulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }
        solver->model().beginReportStep();

        bool enableTUNING = EWOMS_GET_PARAM(TypeTag, bool, EnableTuning);

        // If sub stepping is enabled allow the solver to sub cycle
        // in case the report steps are too large for the solver to converge
        //
        // \Note: The report steps are met in any case
        // \Note: The sub stepping will require a copy of the state variables
        if (adaptiveTimeStepping_) {
            const auto& events = schedule()[timer.currentStepNum()].events();
            if (enableTUNING) {
                if (events.hasEvent(ScheduleEvents::TUNING_CHANGE)) {
                    const auto& sched_state = schedule()[timer.currentStepNum()];
                    const auto& tuning = sched_state.tuning();
                    const auto& max_next_tstep = sched_state.max_next_tstep();
                    adaptiveTimeStepping_->updateTUNING(max_next_tstep, tuning);
                }
            }
            bool event = events.hasEvent(ScheduleEvents::NEW_WELL) ||
                events.hasEvent(ScheduleEvents::INJECTION_TYPE_CHANGED) ||
                events.hasEvent(ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER) ||
                events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE);
            auto stepReport = adaptiveTimeStepping_->step(timer, *solver, event, nullptr);
            report_ += stepReport;
            //Pass simulation report to eclwriter for summary output
            ebosSimulator_.problem().setSimulationReport(report_);
        } else {
            // solve for complete report step
            auto stepReport = solver->step(timer);
            report_ += stepReport;
            if (terminalOutput_) {
                std::ostringstream ss;
                stepReport.reportStep(ss);
                OpmLog::info(ss.str());
            }
        }

        // write simulation state at the report stage
        Dune::Timer perfTimer;
        perfTimer.start();
        const double nextstep = adaptiveTimeStepping_ ? adaptiveTimeStepping_->suggestedNextStep() : -1.0;
        ebosSimulator_.problem().setNextTimeStepSize(nextstep);
        ebosSimulator_.problem().writeOutput();
        report_.success.output_write_time += perfTimer.stop();

        solver->model().endReportStep();

        // take time that was used to solve system for this reportStep
        solverTimer_->stop();

        // update timing.
        report_.success.solver_time += solverTimer_->secsSinceStart();

        if (this->grid().comm().rank() == 0) {
            // Destructively grab the step convergence reports.  The solver
            // object and the model object contained therein are about to go
            // out of scope.
            this->writeConvergenceOutput(solver->model().getStepReportsDestructively());
        }

        // Increment timer, remember well state.
        ++timer;

        if (terminalOutput_) {
            if (!timer.initialStep()) {
                const std::string version = moduleVersionName();
                outputTimestampFIP(timer, eclState().getTitle(), version);
            }

            std::string msg =
                "Time step took " + std::to_string(solverTimer_->secsSinceStart()) + " seconds; "
                "total solver time " + std::to_string(report_.success.solver_time) + " seconds.";
            OpmLog::debug(msg);
        }

        handleSave(timer);

        return true;
    }

    SimulatorReport finalize()
    {
        // make sure all output is written to disk before run is finished
        {
            Dune::Timer finalOutputTimer;
            finalOutputTimer.start();

            ebosSimulator_.problem().finalizeOutput();
            report_.success.output_write_time += finalOutputTimer.stop();
        }

        // Stop timer and create timing report
        totalTimer_->stop();
        report_.success.total_time = totalTimer_->secsSinceStart();
        report_.success.converged = true;

        return report_;
    }

    const Grid& grid() const
    { return ebosSimulator_.vanguard().grid(); }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(ebosSimulator_);
        serializer(report_);
        serializer(adaptiveTimeStepping_);
    }

protected:

    std::unique_ptr<Solver> createSolver(WellModel& wellModel)
    {
        auto model = std::make_unique<Model>(ebosSimulator_,
                                             modelParam_,
                                             wellModel,
                                             terminalOutput_);

        return std::make_unique<Solver>(solverParam_, std::move(model));
    }

    const EclipseState& eclState() const
    { return ebosSimulator_.vanguard().eclState(); }


    const Schedule& schedule() const
    { return ebosSimulator_.vanguard().schedule(); }

    bool isRestart() const
    {
        const auto& initconfig = eclState().getInitConfig();
        return initconfig.restartRequested();
    }

    WellModel& wellModel_()
    { return ebosSimulator_.problem().wellModel(); }

    const WellModel& wellModel_() const
    { return ebosSimulator_.problem().wellModel(); }

    void startConvergenceOutputThread(std::string_view convOutputOptions,
                                      std::string_view optionName)
    {
        const auto config = ConvergenceOutputConfiguration {
            convOutputOptions, optionName
        };
        if (! config.want(ConvergenceOutputConfiguration::Option::Iterations)) {
            return;
        }

        auto getPhaseName = ConvergenceOutputThread::ComponentToPhaseName {
            [compNames = typename Model::ComponentName{}](const int compIdx)
            { return std::string_view { compNames.name(compIdx) }; }
        };

        auto convertTime = ConvergenceOutputThread::ConvertToTimeUnits {
            [usys = this->eclState().getUnits()](const double time)
            { return usys.from_si(UnitSystem::measure::time, time); }
        };

        this->convergenceOutputQueue_.emplace();
        this->convergenceOutputObject_.emplace
            (this->eclState().getIOConfig().getOutputDir(),
             this->eclState().getIOConfig().getBaseName(),
             std::move(getPhaseName),
             std::move(convertTime),
             config, *this->convergenceOutputQueue_);

        this->convergenceOutputThread_
            .emplace(&ConvergenceOutputThread::writeASynchronous,
                     &this->convergenceOutputObject_.value());
    }

    void writeConvergenceOutput(std::vector<StepReport>&& reports)
    {
        if (! this->convergenceOutputThread_.has_value()) {
            return;
        }

        auto requests = std::vector<ConvergenceReportQueue::OutputRequest>{};
        requests.reserve(reports.size());

        for (auto&& report : reports) {
            requests.push_back({ report.report_step, report.current_step, std::move(report.report) });
        }

        this->convergenceOutputQueue_->enqueue(std::move(requests));
    }

    void endConvergenceOutputThread()
    {
        if (! this->convergenceOutputThread_.has_value()) {
            return;
        }

        this->convergenceOutputQueue_->signalLastOutputRequest();
        this->convergenceOutputThread_->join();
    }

    //! \brief Serialization of simulator data to .OPMRST files at end of report steps.
    void handleSave(SimulatorTimer& timer)
    {
        if (saveStride_ == -1 && saveStep_ == -1) {
            return;
        }

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        int nextStep = timer.currentStepNum();
        if ((saveStep_ != -1 && nextStep == saveStep_)  ||
            (saveStride_ != -1 && (nextStep % saveStride_) == 0)) {
#if !HAVE_HDF5
            OpmLog::error("Saving of serialized state requested, but no HDF5 support available.");
#else
            const std::string groupName = "/report_step/" + std::to_string(nextStep);
            if (nextStep == saveStride_ || nextStep == saveStep_) {
                std::filesystem::remove(saveFile_);
            }
            HDF5Serializer writer(saveFile_,
                                  HDF5File::OpenMode::APPEND,
                                  EclGenericVanguard::comm());
            if (nextStep == saveStride_ || nextStep == saveStep_) {
                std::ostringstream str;
                Parameters::printValues<TypeTag>(str);
                writer.writeHeader("OPM Flow",
                                   moduleVersion(),
                                   compileTimestamp(),
                                   ebosSimulator_.vanguard().caseName(),
                                   str.str(),
                                   EclGenericVanguard::comm().size());

                if (EclGenericVanguard::comm().size() > 1) {
                    const auto& cellMapping = ebosSimulator_.vanguard().globalCell();
                    std::size_t hash = Dune::hash_range(cellMapping.begin(), cellMapping.end());
                    writer.write(hash, "/", "grid_checksum");
                }
            }
            writer.write(*this, groupName, "simulator_data");
            writer.write(timer, groupName, "simulator_timer",
                         HDF5File::DataSetMode::ROOT_ONLY);
            OpmLog::info("Serialized state written for report step " + std::to_string(nextStep));
#endif
        }

        OPM_END_PARALLEL_TRY_CATCH("Error saving serialized state: ",
                                   EclGenericVanguard::comm());
    }

    //! \brief Load timer info from serialized state.
    void loadTimerInfo([[maybe_unused]] SimulatorTimer& timer)
    {
#if !HAVE_HDF5
        OpmLog::error("Loading of serialized state requested, but no HDF5 support available.");
        loadStep_ = -1;
#else
        OPM_BEGIN_PARALLEL_TRY_CATCH();

        HDF5Serializer reader(saveFile_,
                              HDF5File::OpenMode::READ,
                              EclGenericVanguard::comm());

        if (loadStep_ == 0) {
            loadStep_ = reader.lastReportStep();
        }

        OpmLog::info("Loading serialized state for report step " + std::to_string(loadStep_));
        const std::string groupName = "/report_step/" + std::to_string(loadStep_);
        reader.read(timer, groupName, "simulator_timer", HDF5File::DataSetMode::ROOT_ONLY);

        if (EclGenericVanguard::comm().size() > 1) {
            std::size_t stored_hash;
            reader.read(stored_hash, "/", "grid_checksum");
            const auto& cellMapping = ebosSimulator_.vanguard().globalCell();
            std::size_t hash = Dune::hash_range(cellMapping.begin(), cellMapping.end());
            if (hash != stored_hash) {
                throw std::runtime_error("Grid hash mismatch, .SAVE file cannot be used");
            }
        }

        OPM_END_PARALLEL_TRY_CATCH("Error loading serialized state: ",
                                   EclGenericVanguard::comm());
#endif
    }

    //! \brief Load simulator state from serialized state.
    void loadSimulatorState()
    {
#if HAVE_HDF5
        OPM_BEGIN_PARALLEL_TRY_CATCH();

        HDF5Serializer reader(saveFile_,
                              HDF5File::OpenMode::READ,
                              EclGenericVanguard::comm());
        const std::string groupName = "/report_step/" + std::to_string(loadStep_);
        reader.read(*this, groupName, "simulator_data");

        OPM_END_PARALLEL_TRY_CATCH("Error loading serialized state: ",
                                   EclGenericVanguard::comm());
#endif
    }

    // Data.
    Simulator& ebosSimulator_;
    std::unique_ptr<WellConnectionAuxiliaryModule<TypeTag>> wellAuxMod_;

    ModelParameters modelParam_;
    SolverParameters solverParam_;

    // Observed objects.
    PhaseUsage phaseUsage_;
    // Misc. data
    bool terminalOutput_;

    SimulatorReport report_;
    std::unique_ptr<time::StopWatch> solverTimer_;
    std::unique_ptr<time::StopWatch> totalTimer_;
    std::unique_ptr<TimeStepper> adaptiveTimeStepping_;

    std::optional<ConvergenceReportQueue> convergenceOutputQueue_{};
    std::optional<ConvergenceOutputThread> convergenceOutputObject_{};
    std::optional<std::thread> convergenceOutputThread_{};

    int saveStride_ = -1; //!< Stride to save serialized state at
    int saveStep_ = -1; //!< Specific step to save serialized state at
    int loadStep_ = -1; //!< Step to load serialized state from
    std::string saveFile_; //!< File to load/save serialized state from/to.
};

} // namespace Opm

#endif // OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_EBOS_HPP
