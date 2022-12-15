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

#include <opm/simulators/flow/NonlinearSolverEbos.hpp>
#include <opm/simulators/flow/BlackoilModelEbos.hpp>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/grid/utility/StopWatch.hpp>

#include <opm/common/ErrorMacros.hpp>

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
struct ExtraConvergenceOutput
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
struct ExtraConvergenceOutput<TypeTag, TTag::EclFlowProblem>
{
    static constexpr auto* value = "none";
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

        // Only rank 0 does print to std::cout
        const auto& comm = grid().comm();
        terminalOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);
        terminalOutput_ = terminalOutput_ && (comm.rank() == 0);
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
        EWOMS_REGISTER_PARAM(TypeTag, std::string, ExtraConvergenceOutput,
                             "Provide additional convergence output "
                             "files for diagnostic purposes. "
                             "\"none\" gives no extra output and "
                             "overrides all other options, "
                             "\"steps\" generates an INFOSTEP file, "
                             "\"iterations\" generates an INFOITER file. "
                             "Combine options with commas, e.g., "
                             "\"steps,iterations\" for multiple outputs.");
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

        // Increment timer, remember well state.
        ++timer;

        if (terminalOutput_) {
            if (!timer.initialStep()) {
                const std::string version = moduleVersionName();
                outputTimestampFIP(timer, eclState().getTitle(), version);
            }
        }

        if (terminalOutput_) {
            std::string msg =
                "Time step took " + std::to_string(solverTimer_->secsSinceStart()) + " seconds; "
                "total solver time " + std::to_string(report_.success.solver_time) + " seconds.";
            OpmLog::debug(msg);
        }
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
};

} // namespace Opm

#endif // OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_EBOS_HPP
