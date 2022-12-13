/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP

#include <iostream>
#include <utility>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/input/eclipse/Schedule/ScheduleState.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>

namespace Opm::Properties {

namespace TTag {
struct FlowTimeSteppingParameters {};
}

template<class TypeTag, class MyTypeTag>
struct SolverRestartFactor {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverGrowthFactor {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverMaxGrowth {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverMaxTimeStepInDays {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverMinTimeStep {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverContinueOnConvergenceFailure {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverMaxRestarts {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SolverVerbosity {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepVerbosity {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct InitialTimeStepInDays {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct FullTimeStepInitially {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepAfterEventInDays {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControl {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlTolerance {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlTargetIterations {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlTargetNewtonIterations {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlDecayRate {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlGrowthRate {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlDecayDampingFactor {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlGrowthDampingFactor {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct TimeStepControlFileName {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MinTimeStepBeforeShuttingProblematicWellsInDays {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct MinTimeStepBasedOnNewtonIterations {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct SolverRestartFactor<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.33;
};
template<class TypeTag>
struct SolverGrowthFactor<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2.0;
};
template<class TypeTag>
struct SolverMaxGrowth<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3.0;
};
template<class TypeTag>
struct SolverMaxTimeStepInDays<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 365.0;
};
template<class TypeTag>
struct SolverMinTimeStep<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0e-12;
};
template<class TypeTag>
struct SolverContinueOnConvergenceFailure<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct SolverMaxRestarts<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr int value = 10;
};
template<class TypeTag>
struct SolverVerbosity<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr int value = 1;
};
template<class TypeTag>
struct TimeStepVerbosity<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr int value = 1;
};
template<class TypeTag>
struct InitialTimeStepInDays<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct FullTimeStepInitially<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct TimeStepAfterEventInDays<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = -1.0;
};
template<class TypeTag>
struct TimeStepControl<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr auto value = "pid+newtoniteration";
};
template<class TypeTag>
struct TimeStepControlTolerance<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-1;
};
template<class TypeTag>
struct TimeStepControlTargetIterations<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr int value = 30;
};
template<class TypeTag>
struct TimeStepControlTargetNewtonIterations<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr int value = 8;
};
template<class TypeTag>
struct TimeStepControlDecayRate<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.75;
};
template<class TypeTag>
struct TimeStepControlGrowthRate<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.25;
};
template<class TypeTag>
struct TimeStepControlDecayDampingFactor<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct TimeStepControlGrowthDampingFactor<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3.2;
};
template<class TypeTag>
struct TimeStepControlFileName<TypeTag, TTag::FlowTimeSteppingParameters> {
    static constexpr auto value = "timesteps";
};
template<class TypeTag>
struct MinTimeStepBeforeShuttingProblematicWellsInDays<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.01;
};

template<class TypeTag>
struct MinTimeStepBasedOnNewtonIterations<TypeTag, TTag::FlowTimeSteppingParameters> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

} // namespace Opm::Properties

namespace Opm {
    // AdaptiveTimeStepping
    //---------------------
    void logTimer(const AdaptiveSimulatorTimer& substepTimer);

    template<class TypeTag>
    class AdaptiveTimeSteppingEbos
    {
        template <class Solver>
        class SolutionTimeErrorSolverWrapperEbos : public RelativeChangeInterface
        {
            const Solver& solver_;
        public:
            SolutionTimeErrorSolverWrapperEbos(const Solver& solver)
              : solver_(solver)
            {}

            /// return || u^n+1 - u^n || / || u^n+1 ||
            double relativeChange() const
            { return solver_.model().relativeChange(); }
        };

        template<class E>
        void logException_(const E& exception, bool verbose)
        {
            if (verbose) {
                std::string message;
                message = "Caught Exception: ";
                message += exception.what();
                OpmLog::debug(message);
            }
        }

    public:
        //! \brief contructor taking parameter object
        AdaptiveTimeSteppingEbos(const UnitSystem& unitSystem,
                                 const bool terminalOutput = true)
            : timeStepControl_()
            , restartFactor_(EWOMS_GET_PARAM(TypeTag, double, SolverRestartFactor)) // 0.33
            , growthFactor_(EWOMS_GET_PARAM(TypeTag, double, SolverGrowthFactor)) // 2.0
            , maxGrowth_(EWOMS_GET_PARAM(TypeTag, double, SolverMaxGrowth)) // 3.0
            , maxTimeStep_(EWOMS_GET_PARAM(TypeTag, double, SolverMaxTimeStepInDays)*24*60*60) // 365.25
            , minTimeStep_(unitSystem.to_si(UnitSystem::measure::time, EWOMS_GET_PARAM(TypeTag, double, SolverMinTimeStep))) // 1e-12;
            , ignoreConvergenceFailure_(EWOMS_GET_PARAM(TypeTag, bool, SolverContinueOnConvergenceFailure)) // false;
            , solverRestartMax_(EWOMS_GET_PARAM(TypeTag, int, SolverMaxRestarts)) // 10
            , solverVerbose_(EWOMS_GET_PARAM(TypeTag, int, SolverVerbosity) > 0 && terminalOutput) // 2
            , timestepVerbose_(EWOMS_GET_PARAM(TypeTag, int, TimeStepVerbosity) > 0 && terminalOutput) // 2
            , suggestedNextTimestep_(EWOMS_GET_PARAM(TypeTag, double, InitialTimeStepInDays)*24*60*60) // 1.0
            , fullTimestepInitially_(EWOMS_GET_PARAM(TypeTag, bool, FullTimeStepInitially)) // false
            , timestepAfterEvent_(EWOMS_GET_PARAM(TypeTag, double, TimeStepAfterEventInDays)*24*60*60) // 1e30
            , useNewtonIteration_(false)
            , minTimeStepBeforeShuttingProblematicWells_(EWOMS_GET_PARAM(TypeTag, double, MinTimeStepBeforeShuttingProblematicWellsInDays)*unit::day)

        {
            init_(unitSystem);
        }



        //! \brief contructor taking parameter object
        //! \param tuning Pointer to ecl TUNING keyword
        //! \param timeStep current report step
        AdaptiveTimeSteppingEbos(double max_next_tstep,
                                 const Tuning& tuning,
                                 const UnitSystem& unitSystem,
                                 const bool terminalOutput = true)
            : timeStepControl_()
            , restartFactor_(tuning.TSFCNV)
            , growthFactor_(tuning.TFDIFF)
            , maxGrowth_(tuning.TSFMAX)
            , maxTimeStep_(tuning.TSMAXZ) // 365.25
            , minTimeStep_(tuning.TSFMIN) // 0.1;
            , ignoreConvergenceFailure_(true)
            , solverRestartMax_(EWOMS_GET_PARAM(TypeTag, int, SolverMaxRestarts)) // 10
            , solverVerbose_(EWOMS_GET_PARAM(TypeTag, int, SolverVerbosity) > 0 && terminalOutput) // 2
            , timestepVerbose_(EWOMS_GET_PARAM(TypeTag, int, TimeStepVerbosity) > 0 && terminalOutput) // 2
            , suggestedNextTimestep_(max_next_tstep) // 1.0
            , fullTimestepInitially_(EWOMS_GET_PARAM(TypeTag, bool, FullTimeStepInitially)) // false
            , timestepAfterEvent_(tuning.TMAXWC) // 1e30
            , useNewtonIteration_(false)
            , minTimeStepBeforeShuttingProblematicWells_(EWOMS_GET_PARAM(TypeTag, double, MinTimeStepBeforeShuttingProblematicWellsInDays)*unit::day)
        {
            init_(unitSystem);
        }

        static void registerParameters()
        {
            // TODO: make sure the help messages are correct (and useful)
            EWOMS_REGISTER_PARAM(TypeTag, double, SolverRestartFactor,
                                 "The factor time steps are elongated after restarts");
            EWOMS_REGISTER_PARAM(TypeTag, double, SolverGrowthFactor,
                                 "The factor time steps are elongated after a successful substep");
            EWOMS_REGISTER_PARAM(TypeTag, double, SolverMaxGrowth,
                                 "The maximum factor time steps are elongated after a report step");
            EWOMS_REGISTER_PARAM(TypeTag, double, SolverMaxTimeStepInDays,
                                 "The maximum size of a time step in days");
            EWOMS_REGISTER_PARAM(TypeTag, double, SolverMinTimeStep,
                                 "The minimum size of a time step in days for field and metric and hours for lab. If a step cannot converge without getting cut below this step size the simulator will stop");
            EWOMS_REGISTER_PARAM(TypeTag, bool, SolverContinueOnConvergenceFailure,
                                 "Continue instead of stop when minimum solver time step is reached");
            EWOMS_REGISTER_PARAM(TypeTag, int, SolverMaxRestarts,
                                 "The maximum number of breakdowns before a substep is given up and the simulator is terminated");
            EWOMS_REGISTER_PARAM(TypeTag, int, SolverVerbosity,
                                 "Specify the \"chattiness\" of the non-linear solver itself");
            EWOMS_REGISTER_PARAM(TypeTag, int, TimeStepVerbosity,
                                 "Specify the \"chattiness\" during the time integration");
            EWOMS_REGISTER_PARAM(TypeTag, double, InitialTimeStepInDays,
                                 "The size of the initial time step in days");
            EWOMS_REGISTER_PARAM(TypeTag, bool, FullTimeStepInitially,
                                 "Always attempt to finish a report step using a single substep");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepAfterEventInDays,
                                 "Time step size of the first time step after an event occurs during the simulation in days");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, TimeStepControl,
                                 "The algorithm used to determine time-step sizes. valid options are: 'pid' (default), 'pid+iteration', 'pid+newtoniteration', 'iterationcount', 'newtoniterationcount' and 'hardcoded'");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepControlTolerance,
                                 "The tolerance used by the time step size control algorithm");
            EWOMS_REGISTER_PARAM(TypeTag, int, TimeStepControlTargetIterations,
                                 "The number of linear iterations which the time step control scheme should aim for (if applicable)");
            EWOMS_REGISTER_PARAM(TypeTag, int, TimeStepControlTargetNewtonIterations,
                                 "The number of Newton iterations which the time step control scheme should aim for (if applicable)");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepControlDecayRate,
                                 "The decay rate of the time step size of the number of target iterations is exceeded");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepControlGrowthRate,
                                 "The growth rate of the time step size of the number of target iterations is undercut");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepControlDecayDampingFactor,
                                 "The decay rate of the time step decrease when the target iterations is exceeded");
            EWOMS_REGISTER_PARAM(TypeTag, double, TimeStepControlGrowthDampingFactor,
                                 "The growth rate of the time step increase when the target iterations is undercut");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, TimeStepControlFileName,
                                 "The name of the file which contains the hardcoded time steps sizes");
            EWOMS_REGISTER_PARAM(TypeTag, double, MinTimeStepBeforeShuttingProblematicWellsInDays,
                                 "The minimum time step size in days for which problematic wells are not shut");
            EWOMS_REGISTER_PARAM(TypeTag, double, MinTimeStepBasedOnNewtonIterations,
                                 "The minimum time step size (in days for field and metric unit and hours for lab unit) can be reduced to based on newton iteration counts");
        }

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps
        */
        template <class Solver>
        SimulatorReport step(const SimulatorTimer& simulatorTimer,
                             Solver& solver,
                             const bool isEvent,
                             const std::vector<int>* fipnum = nullptr)
        {
            SimulatorReport report;
            const double timestep = simulatorTimer.currentStepLength();

            // init last time step as a fraction of the given time step
            if (suggestedNextTimestep_ < 0) {
                suggestedNextTimestep_ = restartFactor_ * timestep;
            }

            if (fullTimestepInitially_) {
                suggestedNextTimestep_ = timestep;
            }

            // use seperate time step after event
            if (isEvent && timestepAfterEvent_ > 0) {
                suggestedNextTimestep_ = timestepAfterEvent_;
            }

            auto& ebosSimulator = solver.model().ebosSimulator();
            auto& ebosProblem = ebosSimulator.problem();

            // create adaptive step timer with previously used sub step size
            AdaptiveSimulatorTimer substepTimer(simulatorTimer, suggestedNextTimestep_, maxTimeStep_);

            // counter for solver restarts
            int restarts = 0;

            // sub step time loop
            while (!substepTimer.done()) {
                // get current delta t
                const double dt = substepTimer.currentStepLength() ;
                if (timestepVerbose_) {
                    logTimer(substepTimer);
                }

                SimulatorReportSingle substepReport;
                std::string causeOfFailure = "";
                try {
                    substepReport = solver.step(substepTimer);
                    if (solverVerbose_) {
                        // report number of linear iterations
                        OpmLog::debug("Overall linear iterations used: " + std::to_string(substepReport.total_linear_iterations));
                    }
                }
                catch (const TooManyIterations& e) {
                    substepReport = solver.failureReport();
                    causeOfFailure = "Solver convergence failure - Iteration limit reached";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const LinearSolverProblem& e) {
                    substepReport = solver.failureReport();
                    causeOfFailure = "Linear solver convergence failure";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const NumericalProblem& e) {
                    substepReport = solver.failureReport();
                    causeOfFailure = "Solver convergence failure - Numerical problem encountered";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const std::runtime_error& e) {
                    substepReport = solver.failureReport();

                    logException_(e, solverVerbose_);
                    // also catch linear solver not converged
                }
                catch (const Dune::ISTLError& e) {
                    substepReport = solver.failureReport();

                    logException_(e, solverVerbose_);
                    // also catch errors in ISTL AMG that occur when time step is too large
                }
                catch (const Dune::MatrixBlockError& e) {
                    substepReport = solver.failureReport();

                    logException_(e, solverVerbose_);
                    // this can be thrown by ISTL's ILU0 in block mode, yet is not an ISTLError
                }

                //Pass substep to eclwriter for summary output
                ebosSimulator.problem().setSubStepReport(substepReport);

                report += substepReport;

                bool continue_on_uncoverged_solution = ignoreConvergenceFailure_ && !substepReport.converged && dt <= minTimeStep_;

                if (continue_on_uncoverged_solution) {
                    const auto msg = std::string("Solver failed to converge but timestep ")
                            + std::to_string(dt) + " is smaller or equal to "
                            + std::to_string(minTimeStep_) + "\n which is the minimum threshold given"
                            +  "by option --solver-min-time-step= \n";
                    if (solverVerbose_) {
                        OpmLog::error(msg);
                    }
                }

                if (substepReport.converged || continue_on_uncoverged_solution) {

                    // advance by current dt
                    ++substepTimer;

                    // create object to compute the time error, simply forwards the call to the model
                    SolutionTimeErrorSolverWrapperEbos<Solver> relativeChange(solver);

                    // compute new time step estimate
                    const int iterations = useNewtonIteration_ ? substepReport.total_newton_iterations
                        : substepReport.total_linear_iterations;
                    double dtEstimate = timeStepControl_->computeTimeStepSize(dt, iterations, relativeChange,
                                                                               substepTimer.simulationTimeElapsed());

                    assert(dtEstimate > 0);
                    // limit the growth of the timestep size by the growth factor
                    dtEstimate = std::min(dtEstimate, double(maxGrowth_ * dt));
                    assert(dtEstimate > 0);
                    // further restrict time step size growth after convergence problems
                    if (restarts > 0) {
                        dtEstimate = std::min(growthFactor_ * dt, dtEstimate);
                        // solver converged, reset restarts counter
                        restarts = 0;
                    }

                    // Further restrict time step size if we are in
                    // prediction mode with THP constraints.
                    if (solver.model().wellModel().hasTHPConstraints()) {
                        const double maxPredictionTHPTimestep = 16.0 * unit::day;
                        dtEstimate = std::min(dtEstimate, maxPredictionTHPTimestep);
                    }
                    assert(dtEstimate > 0);
                    if (timestepVerbose_) {
                        std::ostringstream ss;
                        substepReport.reportStep(ss);
                        OpmLog::info(ss.str());
                    }

                    // write data if outputWriter was provided
                    // if the time step is done we do not need
                    // to write it as this will be done by the simulator
                    // anyway.
                    if (!substepTimer.done()) {
                        if (fipnum) {
                            solver.computeFluidInPlace(*fipnum);
                        }
                        time::StopWatch perfTimer;
                        perfTimer.start();

                        ebosProblem.writeOutput();

                        report.success.output_write_time += perfTimer.secsSinceStart();
                    }

                    // set new time step length
                    substepTimer.provideTimeStepEstimate(dtEstimate);

                    report.success.converged = substepTimer.done();
                    substepTimer.setLastStepFailed(false);

                }
                else { // in case of no convergence
                    substepTimer.setLastStepFailed(true);

                    // If we have restarted (i.e. cut the timestep) too
                    // many times, we have failed and throw an exception.
                    if (restarts >= solverRestartMax_) {
                        const auto msg = std::string("Solver failed to converge after cutting timestep ")
                            + std::to_string(restarts) + " times.";
                        if (solverVerbose_) {
                            OpmLog::error(msg);
                        }
                        OPM_THROW_NOLOG(NumericalProblem, msg);
                    }

                    // The new, chopped timestep.
                    const double newTimeStep = restartFactor_ * dt;


                    // If we have restarted (i.e. cut the timestep) too
                    // much, we have failed and throw an exception.
                    if (newTimeStep < minTimeStep_) {
                        const auto msg = std::string("Solver failed to converge after cutting timestep to ")
                                + std::to_string(minTimeStep_) + "\n which is the minimum threshold given"
                                +  "by option --solver-min-time-step= \n";
                        if (solverVerbose_) {
                            OpmLog::error(msg);
                        }
                        OPM_THROW_NOLOG(NumericalProblem, msg);
                    }

                    // Define utility function for chopping timestep.
                    auto chopTimestep = [&]() {
                        substepTimer.provideTimeStepEstimate(newTimeStep);
                        if (solverVerbose_) {
                            std::string msg;
                            msg = causeOfFailure + "\nTimestep chopped to "
                                + std::to_string(unit::convert::to(substepTimer.currentStepLength(), unit::day)) + " days\n";
                            OpmLog::problem(msg);
                        }
                        ++restarts;
                    };

                    const double minimumChoppedTimestep = minTimeStepBeforeShuttingProblematicWells_;
                    if (newTimeStep > minimumChoppedTimestep) {
                        chopTimestep();
                    } else {
                        // We are below the threshold, and will check if there are any
                        // wells we should close rather than chopping again.
                        std::set<std::string> failing_wells = consistentlyFailingWells(solver.model().stepReports());
                        if (failing_wells.empty()) {
                            // Found no wells to close, chop the timestep as above.
                            chopTimestep();
                        } else {
                            // Close all consistently failing wells.
                            int num_shut_wells = 0;
                            for (const auto& well : failing_wells) {
                                bool was_shut = solver.model().wellModel().forceShutWellByName(well, substepTimer.simulationTimeElapsed());
                                if (was_shut) {
                                    ++num_shut_wells;
                                }
                            }
                            if (num_shut_wells == 0) {
                                // None of the problematic wells were shut.
                                // We must fall back to chopping again.
                                chopTimestep();
                            } else {
                                substepTimer.provideTimeStepEstimate(dt);
                                if (solverVerbose_) {
                                    std::string msg;
                                    msg = "\nProblematic well(s) were shut: ";
                                    for (const auto& well : failing_wells) {
                                        msg += well;
                                        msg += " ";
                                    }
                                    msg += "(retrying timestep)\n";
                                    OpmLog::problem(msg);
                                }
                            }
                        }
                    }
                }
                ebosProblem.setNextTimeStepSize(substepTimer.currentStepLength());
            }

            // store estimated time step for next reportStep
            suggestedNextTimestep_ = substepTimer.currentStepLength();
            if (timestepVerbose_) {
                std::ostringstream ss;
                substepTimer.report(ss);
                ss << "Suggested next step size = " << unit::convert::to(suggestedNextTimestep_, unit::day) << " (days)" << std::endl;
                OpmLog::debug(ss.str());
            }

            if (! std::isfinite(suggestedNextTimestep_)) { // check for NaN
                suggestedNextTimestep_ = timestep;
            }
            return report;
        }

        /** \brief Returns the simulator report for the failed substeps of the last
         *         report step.
         */
        double suggestedNextStep() const
        { return suggestedNextTimestep_; }

        void setSuggestedNextStep(const double x)
        { suggestedNextTimestep_ = x; }

        void updateTUNING(double max_next_tstep, const Tuning& tuning)
        {
            restartFactor_ = tuning.TSFCNV;
            growthFactor_ = tuning.TFDIFF;
            maxGrowth_ = tuning.TSFMAX;
            maxTimeStep_ = tuning.TSMAXZ;
            suggestedNextTimestep_ = max_next_tstep;
            timestepAfterEvent_ = tuning.TMAXWC;
        }


    protected:
        void init_(const UnitSystem& unitSystem)
        {
            // valid are "pid" and "pid+iteration"
            std::string control = EWOMS_GET_PARAM(TypeTag, std::string, TimeStepControl); // "pid"

            const double tol =  EWOMS_GET_PARAM(TypeTag, double, TimeStepControlTolerance); // 1e-1
            if (control == "pid") {
                timeStepControl_ = TimeStepControlType(new PIDTimeStepControl(tol));
            }
            else if (control == "pid+iteration") {
                const int iterations =  EWOMS_GET_PARAM(TypeTag, int, TimeStepControlTargetIterations); // 30
                const double decayDampingFactor = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlDecayDampingFactor); // 1.0
                const double growthDampingFactor = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlGrowthDampingFactor); // 3.2
                timeStepControl_ = TimeStepControlType(new PIDAndIterationCountTimeStepControl(iterations, decayDampingFactor, growthDampingFactor, tol));
            }
            else if (control == "pid+newtoniteration") {
                const int iterations =  EWOMS_GET_PARAM(TypeTag, int, TimeStepControlTargetNewtonIterations); // 8
                const double decayDampingFactor = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlDecayDampingFactor); // 1.0
                const double growthDampingFactor = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlGrowthDampingFactor); // 3.2
                const double nonDimensionalMinTimeStepIterations = EWOMS_GET_PARAM(TypeTag, double, MinTimeStepBasedOnNewtonIterations); // 0.0 by default
                // the min time step can be reduced by the newton iteration numbers
                double minTimeStepReducedByIterations = unitSystem.to_si(UnitSystem::measure::time, nonDimensionalMinTimeStepIterations);
                timeStepControl_ = TimeStepControlType(new PIDAndIterationCountTimeStepControl(iterations, decayDampingFactor,
                                                                      growthDampingFactor, tol, minTimeStepReducedByIterations));
                useNewtonIteration_ = true;
            }
            else if (control == "iterationcount") {
                const int iterations =  EWOMS_GET_PARAM(TypeTag, int, TimeStepControlTargetIterations); // 30
                const double decayrate = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlDecayRate); // 0.75
                const double growthrate = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlGrowthRate); // 1.25
                timeStepControl_ = TimeStepControlType(new SimpleIterationCountTimeStepControl(iterations, decayrate, growthrate));
            }
            else if (control == "newtoniterationcount") {
                const int iterations =  EWOMS_GET_PARAM(TypeTag, int, TimeStepControlTargetNewtonIterations); // 8
                const double decayrate = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlDecayRate); // 0.75
                const double growthrate = EWOMS_GET_PARAM(TypeTag, double, TimeStepControlGrowthRate); // 1.25
                timeStepControl_ = TimeStepControlType(new SimpleIterationCountTimeStepControl(iterations, decayrate, growthrate));
                useNewtonIteration_ = true;
            }
            else if (control == "hardcoded") {
                const std::string filename = EWOMS_GET_PARAM(TypeTag, std::string, TimeStepControlFileName); // "timesteps"
                timeStepControl_ = TimeStepControlType(new HardcodedTimeStepControl(filename));

            }
            else
                OPM_THROW(std::runtime_error,"Unsupported time step control selected "<< control);

            // make sure growth factor is something reasonable
            assert(growthFactor_ >= 1.0);
        }

        template <class ProblemType>
        std::set<std::string> consistentlyFailingWells(const std::vector<ProblemType>& sr)
        {
            // If there are wells that cause repeated failures, we
            // close them, and restart the un-chopped timestep.
            std::ostringstream msg;
            msg << "    Excessive chopping detected in report step "
                << sr.back().report_step << ", substep " << sr.back().current_step << "\n";

            std::set<std::string> failing_wells;

            // return empty set if no report exists
            // well failures in assembly is not yet registred
            if(sr.back().report.empty())
                return failing_wells;

            const auto& wfs = sr.back().report.back().wellFailures();
            for (const auto& wf : wfs) {
                msg << "        Well that failed: " << wf.wellName() << "\n";
            }
            msg.flush();
            OpmLog::debug(msg.str());

            // Check the last few step reports.
            const int num_steps = 3;
            const int rep_step = sr.back().report_step;
            const int sub_step = sr.back().current_step;
            const int sr_size = sr.size();
            if (sr_size >= num_steps) {
                for (const auto& wf : wfs) {
                    failing_wells.insert(wf.wellName());
                }
                for (int s = 1; s < num_steps; ++s) {
                    const auto& srep = sr[sr_size - 1 - s];
                    // Report must be from same report step and substep, otherwise we have
                    // not chopped/retried enough times on this step.
                    if (srep.report_step != rep_step || srep.current_step != sub_step) {
                        break;
                    }
                    // Get the failing wells for this step, that also failed all other steps.
                    std::set<std::string> failing_wells_step;
                    for (const auto& wf : srep.report.back().wellFailures()) {
                        if (failing_wells.count(wf.wellName()) > 0) {
                            failing_wells_step.insert(wf.wellName());
                        }
                    }
                    failing_wells.swap(failing_wells_step);
                }
            }
            return failing_wells;
        }

        typedef std::unique_ptr<TimeStepControlInterface> TimeStepControlType;

        TimeStepControlType timeStepControl_; //!< time step control object
        double restartFactor_;               //!< factor to multiply time step with when solver fails to converge
        double growthFactor_;                //!< factor to multiply time step when solver recovered from failed convergence
        double maxGrowth_;                   //!< factor that limits the maximum growth of a time step
        double maxTimeStep_;                //!< maximal allowed time step size in days
        double minTimeStep_;                //!< minimal allowed time step size before throwing
        bool ignoreConvergenceFailure_;     //!< continue instead of stop when minimum time step is reached
        int solverRestartMax_;        //!< how many restart of solver are allowed
        bool solverVerbose_;           //!< solver verbosity
        bool timestepVerbose_;         //!< timestep verbosity
        double suggestedNextTimestep_;      //!< suggested size of next timestep
        bool fullTimestepInitially_;        //!< beginning with the size of the time step from data file
        double timestepAfterEvent_;         //!< suggested size of timestep after an event
        bool useNewtonIteration_;           //!< use newton iteration count for adaptive time step control
        double minTimeStepBeforeShuttingProblematicWells_; //! < shut problematic wells when time step size in days are less than this
    };
}

#endif // OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP
