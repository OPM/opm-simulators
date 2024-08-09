/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_HPP

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 8)
#include <dune/istl/istlexception.hh>
#else
#include <dune/istl/ilu.hh>
#endif

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/input/eclipse/Schedule/Tuning.hpp>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>

#include <opm/simulators/utils/phaseUsageFromDeck.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Properties::TTag {

struct FlowTimeSteppingParameters {
  using InheritsFrom = std::tuple<EclTimeSteppingParameters>;
};

}

namespace Opm::Parameters {

template<class TypeTag, class MyTypeTag>
struct SolverContinueOnConvergenceFailure { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct SolverMaxRestarts { using type = Properties::UndefinedProperty; };
template<class TypeTag, class MyTypeTag>

struct SolverVerbosity { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepVerbosity { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct InitialTimeStepInDays { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct FullTimeStepInitially { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControl { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlTolerance { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlTargetIterations { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlTargetNewtonIterations { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlDecayRate { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlGrowthRate { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlDecayDampingFactor { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlGrowthDampingFactor { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TimeStepControlFileName { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MinTimeStepBeforeShuttingProblematicWellsInDays { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MinTimeStepBasedOnNewtonIterations { using type = Properties::UndefinedProperty; };

template<class TypeTag>
struct SolverContinueOnConvergenceFailure<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr bool value = false; };

template<class TypeTag>
struct SolverMaxRestarts<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr int value = 10; };

template<class TypeTag>
struct SolverVerbosity<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr int value = 1; };

template<class TypeTag>
struct TimeStepVerbosity<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr int value = 1; };

template<class TypeTag>
struct InitialTimeStepInDays<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct FullTimeStepInitially<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr bool value = false; };

template<class TypeTag>
struct TimeStepControl<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr auto value = "pid+newtoniteration"; };

template<class TypeTag>
struct TimeStepControlTolerance<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1e-1;
};

template<class TypeTag>
struct TimeStepControlTargetIterations<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr int value = 30; };

template<class TypeTag>
struct TimeStepControlTargetNewtonIterations<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr int value = 8; };

template<class TypeTag>
struct TimeStepControlDecayRate<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 0.75;
};

template<class TypeTag>
struct TimeStepControlGrowthRate<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1.25;
};

template<class TypeTag>
struct TimeStepControlDecayDampingFactor<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct TimeStepControlGrowthDampingFactor<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 3.2;
};

template<class TypeTag>
struct TimeStepControlFileName<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{ static constexpr auto value = "timesteps"; };

template<class TypeTag>
struct MinTimeStepBeforeShuttingProblematicWellsInDays<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 0.01;
};

template<class TypeTag>
struct MinTimeStepBasedOnNewtonIterations<TypeTag, Properties::TTag::FlowTimeSteppingParameters>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = 0.0;
};

} // namespace Opm::Parameters

namespace Opm {

struct StepReport;

namespace detail {

void logTimer(const AdaptiveSimulatorTimer& substepTimer);

std::set<std::string> consistentlyFailingWells(const std::vector<StepReport>& sr);

}

    // AdaptiveTimeStepping
    //---------------------
    template<class TypeTag>
    class AdaptiveTimeStepping
    {
        template <class Solver>
        class SolutionTimeErrorSolverWrapper : public RelativeChangeInterface
        {
            const Solver& solver_;
        public:
            SolutionTimeErrorSolverWrapper(const Solver& solver)
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
        AdaptiveTimeStepping() = default;

        //! \brief contructor taking parameter object
        AdaptiveTimeStepping(const UnitSystem& unitSystem,
                             const double max_next_tstep = -1.0,
                             const bool terminalOutput = true)
            : timeStepControl_()
            , restartFactor_(Parameters::get<TypeTag, Parameters::SolverRestartFactor>()) // 0.33
            , growthFactor_(Parameters::get<TypeTag, Parameters::SolverGrowthFactor>()) // 2.0
            , maxGrowth_(Parameters::get<TypeTag, Parameters::SolverMaxGrowth>()) // 3.0
            , maxTimeStep_(Parameters::get<TypeTag, Parameters::SolverMaxTimeStepInDays>() * 24 * 60 * 60) // 365.25
            , minTimeStep_(unitSystem.to_si(UnitSystem::measure::time, Parameters::get<TypeTag, Parameters::SolverMinTimeStep>())) // 1e-12;
            , ignoreConvergenceFailure_(Parameters::get<TypeTag, Parameters::SolverContinueOnConvergenceFailure>()) // false;
            , solverRestartMax_(Parameters::get<TypeTag, Parameters::SolverMaxRestarts>()) // 10
            , solverVerbose_(Parameters::get<TypeTag, Parameters::SolverVerbosity>() > 0 && terminalOutput) // 2
            , timestepVerbose_(Parameters::get<TypeTag, Parameters::TimeStepVerbosity>() > 0 && terminalOutput) // 2
            , suggestedNextTimestep_((max_next_tstep <= 0 ? Parameters::get<TypeTag, Parameters::InitialTimeStepInDays>() : max_next_tstep) * 24 * 60 * 60) // 1.0
            , fullTimestepInitially_(Parameters::get<TypeTag, Parameters::FullTimeStepInitially>()) // false
            , timestepAfterEvent_(Parameters::get<TypeTag, Parameters::TimeStepAfterEventInDays>() * 24 * 60 * 60) // 1e30
            , useNewtonIteration_(false)
            , minTimeStepBeforeShuttingProblematicWells_(Parameters::get<TypeTag, Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day)

        {
            init_(unitSystem);
        }



        //! \brief contructor taking parameter object
        //! \param tuning Pointer to ecl TUNING keyword
        //! \param timeStep current report step
        AdaptiveTimeStepping(double max_next_tstep,
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
            , solverRestartMax_(Parameters::get<TypeTag, Parameters::SolverMaxRestarts>()) // 10
            , solverVerbose_(Parameters::get<TypeTag, Parameters::SolverVerbosity>() > 0 && terminalOutput) // 2
            , timestepVerbose_(Parameters::get<TypeTag, Parameters::TimeStepVerbosity>() > 0 && terminalOutput) // 2
            , suggestedNextTimestep_(max_next_tstep <= 0 ? Parameters::get<TypeTag, Parameters::InitialTimeStepInDays>() * 24 * 60 * 60 : max_next_tstep) // 1.0
            , fullTimestepInitially_(Parameters::get<TypeTag, Parameters::FullTimeStepInitially>()) // false
            , timestepAfterEvent_(tuning.TMAXWC) // 1e30
            , useNewtonIteration_(false)
            , minTimeStepBeforeShuttingProblematicWells_(Parameters::get<TypeTag, Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day)
        {
            init_(unitSystem);
        }

        static void registerParameters()
        {
            registerEclTimeSteppingParameters<TypeTag>();
            // TODO: make sure the help messages are correct (and useful)
            Parameters::registerParam<TypeTag, Parameters::SolverContinueOnConvergenceFailure>
                ("Continue instead of stop when minimum solver time step is reached");
            Parameters::registerParam<TypeTag, Parameters::SolverMaxRestarts>
                ("The maximum number of breakdowns before a substep is given up and "
                 "the simulator is terminated");
            Parameters::registerParam<TypeTag, Parameters::SolverVerbosity>
                ("Specify the \"chattiness\" of the non-linear solver itself");
            Parameters::registerParam<TypeTag, Parameters::TimeStepVerbosity>
                ("Specify the \"chattiness\" during the time integration");
            Parameters::registerParam<TypeTag, Parameters::InitialTimeStepInDays>
                ("The size of the initial time step in days");
            Parameters::registerParam<TypeTag, Parameters::FullTimeStepInitially>
                ("Always attempt to finish a report step using a single substep");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControl>
                ("The algorithm used to determine time-step sizes. "
                 "Valid options are: "
                 "'pid' (default), "
                 "'pid+iteration', "
                 "'pid+newtoniteration', "
                 "'iterationcount', "
                "'newtoniterationcount' "
                "and 'hardcoded'");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlTolerance>
                ("The tolerance used by the time step size control algorithm");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlTargetIterations>
                ("The number of linear iterations which the time step control scheme "
                 "should aim for (if applicable)");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlTargetNewtonIterations>
                ("The number of Newton iterations which the time step control scheme "
                 "should aim for (if applicable)");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlDecayRate>
                ("The decay rate of the time step size of the number of "
                 "target iterations is exceeded");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlGrowthRate>
                ("The growth rate of the time step size of the number of "
                 "target iterations is undercut");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlDecayDampingFactor>
                ("The decay rate of the time step decrease when the "
                 "target iterations is exceeded");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlGrowthDampingFactor>
                ("The growth rate of the time step increase when the "
                 "target iterations is undercut");
            Parameters::registerParam<TypeTag, Parameters::TimeStepControlFileName>
                ("The name of the file which contains the hardcoded time steps sizes");
            Parameters::registerParam<TypeTag, Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>
                ("The minimum time step size in days for which problematic wells are not shut");
            Parameters::registerParam<TypeTag, Parameters::MinTimeStepBasedOnNewtonIterations>
                ("The minimum time step size (in days for field and metric unit and hours for lab unit) "
                 "can be reduced to based on newton iteration counts");
        }

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps
            \param tuningUpdater Function used to update TUNING parameters before each
                                 time step. ACTIONX might change tuning.
        */
        template <class Solver>
        SimulatorReport step(const SimulatorTimer& simulatorTimer,
                             Solver& solver,
                             const bool isEvent,
                             const std::vector<int>* fipnum = nullptr,
                             const std::function<bool()> tuningUpdater = [](){return false;})
        {
            // Maybe update tuning
            tuningUpdater();
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

            auto& simulator = solver.model().simulator();
            auto& problem = simulator.problem();

            // create adaptive step timer with previously used sub step size
            AdaptiveSimulatorTimer substepTimer(simulatorTimer, suggestedNextTimestep_, maxTimeStep_);

            // counter for solver restarts
            int restarts = 0;

            // sub step time loop
            while (!substepTimer.done()) {
                // Maybe update tuning
                // get current delta t
                auto oldValue = suggestedNextTimestep_;
                if (tuningUpdater()) {
                    // Use provideTimeStepEstimate to make we sure don't simulate longer than the report step is.
                    substepTimer.provideTimeStepEstimate(suggestedNextTimestep_);
                    suggestedNextTimestep_ = oldValue;
                }
                const double dt = substepTimer.currentStepLength();
                if (timestepVerbose_) {
                    detail::logTimer(substepTimer);
                }

                SimulatorReportSingle substepReport;
                std::string causeOfFailure;
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
                simulator.problem().setSubStepReport(substepReport);

                report += substepReport;

                bool continue_on_uncoverged_solution = ignoreConvergenceFailure_ &&
                                                       !substepReport.converged  &&
                                                       dt <= minTimeStep_;

                if (continue_on_uncoverged_solution && solverVerbose_) {
                    const auto msg = fmt::format(
                        "Solver failed to converge but timestep "
                        "{} is smaller or equal to {}\n"
                        "which is the minimum threshold given "
                        "by option --solver-min-time-step\n",
                        dt, minTimeStep_
                    );
                    OpmLog::problem(msg);
                }

                if (substepReport.converged || continue_on_uncoverged_solution) {

                    // advance by current dt
                    ++substepTimer;

                    // create object to compute the time error, simply forwards the call to the model
                    SolutionTimeErrorSolverWrapper<Solver> relativeChange(solver);

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

                        problem.writeOutput(simulatorTimer);

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
                        const auto msg = fmt::format(
                            "Solver failed to converge after cutting timestep {} times.",
                            restarts
                        );
                        if (solverVerbose_) {
                            OpmLog::error(msg);
                        }
                        // Use throw directly to prevent file and line
                        throw TimeSteppingBreakdown{msg};
                    }

                    // The new, chopped timestep.
                    const double newTimeStep = restartFactor_ * dt;


                    // If we have restarted (i.e. cut the timestep) too
                    // much, we have failed and throw an exception.
                    if (newTimeStep < minTimeStep_) {
                        const auto msg = fmt::format(
                            "Solver failed to converge after cutting timestep to {}\n"
                            "which is the minimum threshold given by option --solver-min-time-step\n",
                            minTimeStep_
                        );
                        if (solverVerbose_) {
                            OpmLog::error(msg);
                        }
                        // Use throw directly to prevent file and line
                        throw TimeSteppingBreakdown{msg};
                    }

                    // Define utility function for chopping timestep.
                    auto chopTimestep = [&]() {
                        substepTimer.provideTimeStepEstimate(newTimeStep);
                        if (solverVerbose_) {
                            const auto msg = fmt::format(
                                "{}\nTimestep chopped to {} days\n",
                                causeOfFailure,
                                std::to_string(unit::convert::to(substepTimer.currentStepLength(), unit::day))
                            );
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
                        std::set<std::string> failing_wells = detail::consistentlyFailingWells(solver.model().stepReports());
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
                problem.setNextTimeStepSize(substepTimer.currentStepLength());
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
            updateNEXTSTEP(max_next_tstep);
            timestepAfterEvent_ = tuning.TMAXWC;
        }

        void updateNEXTSTEP(double max_next_tstep)
        {
             // \Note Only update next suggested step if TSINIT was explicitly set in TUNING or NEXTSTEP is active. 
            if (max_next_tstep > 0) {
                suggestedNextTimestep_ = max_next_tstep;
            }
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(timeStepControlType_);
            switch (timeStepControlType_) {
            case TimeStepControlType::HardCodedTimeStep:
                allocAndSerialize<HardcodedTimeStepControl>(serializer);
                break;
            case TimeStepControlType::PIDAndIterationCount:
                allocAndSerialize<PIDAndIterationCountTimeStepControl>(serializer);
                break;
            case TimeStepControlType::SimpleIterationCount:
                allocAndSerialize<SimpleIterationCountTimeStepControl>(serializer);
                break;
            case TimeStepControlType::PID:
                allocAndSerialize<PIDTimeStepControl>(serializer);
                break;
            }
            serializer(restartFactor_);
            serializer(growthFactor_);
            serializer(maxGrowth_);
            serializer(maxTimeStep_);
            serializer(minTimeStep_);
            serializer(ignoreConvergenceFailure_);
            serializer(solverRestartMax_);
            serializer(solverVerbose_);
            serializer(timestepVerbose_);
            serializer(suggestedNextTimestep_);
            serializer(fullTimestepInitially_);
            serializer(timestepAfterEvent_);
            serializer(useNewtonIteration_);
            serializer(minTimeStepBeforeShuttingProblematicWells_);
        }

        static AdaptiveTimeStepping<TypeTag> serializationTestObjectHardcoded()
        {
            return serializationTestObject_<HardcodedTimeStepControl>();
        }

        static AdaptiveTimeStepping<TypeTag> serializationTestObjectPID()
        {
            return serializationTestObject_<PIDTimeStepControl>();
        }

        static AdaptiveTimeStepping<TypeTag> serializationTestObjectPIDIt()
        {
            return serializationTestObject_<PIDAndIterationCountTimeStepControl>();
        }

        static AdaptiveTimeStepping<TypeTag> serializationTestObjectSimple()
        {
            return serializationTestObject_<SimpleIterationCountTimeStepControl>();
        }

        bool operator==(const AdaptiveTimeStepping<TypeTag>& rhs)
        {
            if (timeStepControlType_ != rhs.timeStepControlType_ ||
                (timeStepControl_ && !rhs.timeStepControl_) ||
                (!timeStepControl_ && rhs.timeStepControl_)) {
                return false;
            }

            bool result = false;
            switch (timeStepControlType_) {
            case TimeStepControlType::HardCodedTimeStep:
                result = castAndComp<HardcodedTimeStepControl>(rhs);
                break;
            case TimeStepControlType::PIDAndIterationCount:
                result = castAndComp<PIDAndIterationCountTimeStepControl>(rhs);
                break;
            case TimeStepControlType::SimpleIterationCount:
                result = castAndComp<SimpleIterationCountTimeStepControl>(rhs);
                break;
            case TimeStepControlType::PID:
                result = castAndComp<PIDTimeStepControl>(rhs);
                break;
            }

            return result &&
                   this->restartFactor_ == rhs.restartFactor_ &&
                   this->growthFactor_ == rhs.growthFactor_ &&
                   this->maxGrowth_ == rhs.maxGrowth_ &&
                   this->maxTimeStep_ == rhs.maxTimeStep_ &&
                   this->minTimeStep_ == rhs.minTimeStep_ &&
                   this->ignoreConvergenceFailure_ == rhs.ignoreConvergenceFailure_ &&
                   this->solverRestartMax_== rhs.solverRestartMax_ &&
                   this->solverVerbose_ == rhs.solverVerbose_ &&
                   this->fullTimestepInitially_ == rhs.fullTimestepInitially_ &&
                   this->timestepAfterEvent_ == rhs.timestepAfterEvent_ &&
                   this->useNewtonIteration_ == rhs.useNewtonIteration_ &&
                   this->minTimeStepBeforeShuttingProblematicWells_ ==
                       rhs.minTimeStepBeforeShuttingProblematicWells_;
        }

    private:
        template<class Controller>
        static AdaptiveTimeStepping<TypeTag> serializationTestObject_()
        {
            AdaptiveTimeStepping<TypeTag> result;

            result.restartFactor_ = 1.0;
            result.growthFactor_ = 2.0;
            result.maxGrowth_ = 3.0;
            result.maxTimeStep_ = 4.0;
            result.minTimeStep_ = 5.0;
            result.ignoreConvergenceFailure_ = true;
            result.solverRestartMax_ = 6;
            result.solverVerbose_ = true;
            result.timestepVerbose_ = true;
            result.suggestedNextTimestep_ = 7.0;
            result.fullTimestepInitially_ = true;
            result.useNewtonIteration_ = true;
            result.minTimeStepBeforeShuttingProblematicWells_ = 9.0;
            result.timeStepControlType_ = Controller::Type;
            result.timeStepControl_ = std::make_unique<Controller>(Controller::serializationTestObject());

            return result;
        }
        template<class T, class Serializer>
        void allocAndSerialize(Serializer& serializer)
        {
            if (!serializer.isSerializing()) {
                timeStepControl_ = std::make_unique<T>();
            }
            serializer(*static_cast<T*>(timeStepControl_.get()));
        }

        template<class T>
        bool castAndComp(const AdaptiveTimeStepping<TypeTag>& Rhs) const
        {
            const T* lhs = static_cast<const T*>(timeStepControl_.get());
            const T* rhs = static_cast<const T*>(Rhs.timeStepControl_.get());
            return *lhs == *rhs;
        }

    protected:
        void init_(const UnitSystem& unitSystem)
        {
            // valid are "pid" and "pid+iteration"
            std::string control = Parameters::get<TypeTag, Parameters::TimeStepControl>(); // "pid"

            const double tol =  Parameters::get<TypeTag, Parameters::TimeStepControlTolerance>(); // 1e-1
            if (control == "pid") {
                timeStepControl_ = std::make_unique<PIDTimeStepControl>(tol);
                timeStepControlType_ = TimeStepControlType::PID;
            }
            else if (control == "pid+iteration") {
                const int iterations =  Parameters::get<TypeTag, Parameters::TimeStepControlTargetIterations>(); // 30
                const double decayDampingFactor = Parameters::get<TypeTag, Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
                const double growthDampingFactor = Parameters::get<TypeTag, Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
                timeStepControl_ = std::make_unique<PIDAndIterationCountTimeStepControl>(iterations, decayDampingFactor, growthDampingFactor, tol);
                timeStepControlType_ = TimeStepControlType::PIDAndIterationCount;
            }
            else if (control == "pid+newtoniteration") {
                const int iterations =  Parameters::get<TypeTag, Parameters::TimeStepControlTargetNewtonIterations>(); // 8
                const double decayDampingFactor = Parameters::get<TypeTag, Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
                const double growthDampingFactor = Parameters::get<TypeTag, Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
                const double nonDimensionalMinTimeStepIterations = Parameters::get<TypeTag, Parameters::MinTimeStepBasedOnNewtonIterations>(); // 0.0 by default
                // the min time step can be reduced by the newton iteration numbers
                double minTimeStepReducedByIterations = unitSystem.to_si(UnitSystem::measure::time, nonDimensionalMinTimeStepIterations);
                timeStepControl_ = std::make_unique<PIDAndIterationCountTimeStepControl>(iterations, decayDampingFactor,
                                                                                         growthDampingFactor, tol, minTimeStepReducedByIterations);
                timeStepControlType_ = TimeStepControlType::PIDAndIterationCount;
                useNewtonIteration_ = true;
            }
            else if (control == "iterationcount") {
                const int iterations =  Parameters::get<TypeTag, Parameters::TimeStepControlTargetIterations>(); // 30
                const double decayrate = Parameters::get<TypeTag, Parameters::TimeStepControlDecayRate>(); // 0.75
                const double growthrate = Parameters::get<TypeTag, Parameters::TimeStepControlGrowthRate>(); // 1.25
                timeStepControl_ = std::make_unique<SimpleIterationCountTimeStepControl>(iterations, decayrate, growthrate);
                timeStepControlType_ = TimeStepControlType::SimpleIterationCount;
            }
            else if (control == "newtoniterationcount") {
                const int iterations =  Parameters::get<TypeTag, Parameters::TimeStepControlTargetNewtonIterations>(); // 8
                const double decayrate = Parameters::get<TypeTag, Parameters::TimeStepControlDecayRate>(); // 0.75
                const double growthrate = Parameters::get<TypeTag, Parameters::TimeStepControlGrowthRate>(); // 1.25
                timeStepControl_ = std::make_unique<SimpleIterationCountTimeStepControl>(iterations, decayrate, growthrate);
                useNewtonIteration_ = true;
                timeStepControlType_ = TimeStepControlType::SimpleIterationCount;
            }
            else if (control == "hardcoded") {
                const std::string filename = Parameters::get<TypeTag, Parameters::TimeStepControlFileName>(); // "timesteps"
                timeStepControl_ = std::make_unique<HardcodedTimeStepControl>(filename);
                timeStepControlType_ = TimeStepControlType::HardCodedTimeStep;
            }
            else
                OPM_THROW(std::runtime_error,
                          "Unsupported time step control selected " + control);

            // make sure growth factor is something reasonable
            assert(growthFactor_ >= 1.0);
        }

        using TimeStepController = std::unique_ptr<TimeStepControlInterface>;

        TimeStepControlType timeStepControlType_; //!< type of time step control object
        TimeStepController timeStepControl_; //!< time step control object
        double restartFactor_;               //!< factor to multiply time step with when solver fails to converge
        double growthFactor_;                //!< factor to multiply time step when solver recovered from failed convergence
        double maxGrowth_;                   //!< factor that limits the maximum growth of a time step
        double maxTimeStep_;                 //!< maximal allowed time step size in days
        double minTimeStep_;                 //!< minimal allowed time step size before throwing
        bool ignoreConvergenceFailure_;      //!< continue instead of stop when minimum time step is reached
        int solverRestartMax_;               //!< how many restart of solver are allowed
        bool solverVerbose_;                 //!< solver verbosity
        bool timestepVerbose_;               //!< timestep verbosity
        double suggestedNextTimestep_;       //!< suggested size of next timestep
        bool fullTimestepInitially_;         //!< beginning with the size of the time step from data file
        double timestepAfterEvent_;          //!< suggested size of timestep after an event
        bool useNewtonIteration_;            //!< use newton iteration count for adaptive time step control
        double minTimeStepBeforeShuttingProblematicWells_; //! < shut problematic wells when time step size in days are less than this
    };
}

#endif // OPM_ADAPTIVE_TIME_STEPPING_HPP
