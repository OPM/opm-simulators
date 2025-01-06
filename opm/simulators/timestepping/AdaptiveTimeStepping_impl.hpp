/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_IMPL_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_IMPL_HPP

// Improve IDE experience
#ifndef OPM_ADAPTIVE_TIME_STEPPING_HPP
#include <config.h>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#endif

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace Opm {

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>::
AdaptiveTimeStepping(const UnitSystem& unitSystem,
                     const double max_next_tstep,
                     const bool terminalOutput)
    : timeStepControl_()
    , restartFactor_(Parameters::Get<Parameters::SolverRestartFactor<Scalar>>()) // 0.33
    , growthFactor_(Parameters::Get<Parameters::SolverGrowthFactor<Scalar>>()) // 2.0
    , maxGrowth_(Parameters::Get<Parameters::SolverMaxGrowth<Scalar>>()) // 3.0
    , maxTimeStep_(Parameters::Get<Parameters::SolverMaxTimeStepInDays<Scalar>>() * 24 * 60 * 60) // 365.25
    , minTimeStep_(unitSystem.to_si(UnitSystem::measure::time, Parameters::Get<Parameters::SolverMinTimeStep<Scalar>>())) // 1e-12;
    , ignoreConvergenceFailure_(Parameters::Get<Parameters::SolverContinueOnConvergenceFailure>()) // false;
    , solverRestartMax_(Parameters::Get<Parameters::SolverMaxRestarts>()) // 10
    , solverVerbose_(Parameters::Get<Parameters::SolverVerbosity>() > 0 && terminalOutput) // 2
    , timestepVerbose_(Parameters::Get<Parameters::TimeStepVerbosity>() > 0 && terminalOutput) // 2
    , suggestedNextTimestep_((max_next_tstep <= 0 ? Parameters::Get<Parameters::InitialTimeStepInDays>() : max_next_tstep) * 24 * 60 * 60) // 1.0
    , fullTimestepInitially_(Parameters::Get<Parameters::FullTimeStepInitially>()) // false
    , timestepAfterEvent_(Parameters::Get<Parameters::TimeStepAfterEventInDays<Scalar>>() * 24 * 60 * 60) // 1e30
    , useNewtonIteration_(false)
    , minTimeStepBeforeShuttingProblematicWells_(Parameters::Get<Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day)

{
    init_(unitSystem);
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>::
AdaptiveTimeStepping(double max_next_tstep,
                     const Tuning& tuning,
                     const UnitSystem& unitSystem,
                     const bool terminalOutput)
    : timeStepControl_()
    , restartFactor_(tuning.TSFCNV)
    , growthFactor_(tuning.TFDIFF)
    , maxGrowth_(tuning.TSFMAX)
    , maxTimeStep_(tuning.TSMAXZ) // 365.25
    , minTimeStep_(tuning.TSFMIN) // 0.1;
    , ignoreConvergenceFailure_(true)
    , solverRestartMax_(Parameters::Get<Parameters::SolverMaxRestarts>()) // 10
    , solverVerbose_(Parameters::Get<Parameters::SolverVerbosity>() > 0 && terminalOutput) // 2
    , timestepVerbose_(Parameters::Get<Parameters::TimeStepVerbosity>() > 0 && terminalOutput) // 2
    , suggestedNextTimestep_(max_next_tstep <= 0 ? Parameters::Get<Parameters::InitialTimeStepInDays>() * 24 * 60 * 60 : max_next_tstep) // 1.0
    , fullTimestepInitially_(Parameters::Get<Parameters::FullTimeStepInitially>()) // false
    , timestepAfterEvent_(tuning.TMAXWC) // 1e30
    , useNewtonIteration_(false)
    , minTimeStepBeforeShuttingProblematicWells_(Parameters::Get<Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day)
{
    init_(unitSystem);
}

template<class TypeTag>
void AdaptiveTimeStepping<TypeTag>::registerParameters()
{
    registerEclTimeSteppingParameters<Scalar>();
    detail::registerAdaptiveParameters();
}

template<class TypeTag>
template<class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::
step(const SimulatorTimer& simulatorTimer,
     Solver& solver,
     const bool isEvent,
     const std::function<bool(const double, const double, const int)> tuningUpdater)
{
    // Maybe update tuning
    tuningUpdater(simulatorTimer.simulationTimeElapsed(), suggestedNextTimestep_, 0);
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
        if (tuningUpdater(substepTimer.simulationTimeElapsed(),
                          substepTimer.currentStepLength(),
                          substepTimer.currentStepNum())) {
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
                OpmLog::debug("Overall linear iterations used: " +
                              std::to_string(substepReport.total_linear_iterations));
            }
        }
        catch (const TooManyIterations& e) {
            substepReport = solver.failureReport();
            causeOfFailure = "Solver convergence failure - Iteration limit reached";

            logException_(e, solverVerbose_);
            // since linearIterations is < 0 this will restart the solver
        }
        catch (const ConvergenceMonitorFailure& e) {
            substepReport = solver.failureReport();
            causeOfFailure = "Convergence monitor failure";
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
                time::StopWatch perfTimer;
                perfTimer.start();

                problem.writeOutput(true);

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
                        unit::convert::to(substepTimer.currentStepLength(), unit::day)
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
                    // Close all consistently failing wells that are not under group control
                    std::vector<std::string> shut_wells;
                    for (const auto& well : failing_wells) {
                        bool was_shut = solver.model().wellModel().forceShutWellByName(
                                                    well, substepTimer.simulationTimeElapsed(), /*dont_shut_grup_wells =*/ true);
                        if (was_shut) {
                            shut_wells.push_back(well);
                        }
                    }
                    // If no wells are closed we also try to shut wells under group control
                    if (shut_wells.empty()) {
                        for (const auto& well : failing_wells) {
                            bool was_shut = solver.model().wellModel().forceShutWellByName(
                                                    well, substepTimer.simulationTimeElapsed(), /*dont_shut_grup_wells =*/ false);
                            if (was_shut) {
                                shut_wells.push_back(well);
                            }
                        }
                    }
                    // If still no wells are closed we must fall back to chopping again
                    if (shut_wells.empty()) {
                        chopTimestep();
                    } else {
                        substepTimer.provideTimeStepEstimate(dt);
                        if (solverVerbose_) {
                            const std::string msg =
                                fmt::format("\nProblematic well(s) were shut: {}"
                                            "(retrying timestep)\n",
                                            fmt::join(shut_wells, " "));
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

template<class TypeTag>
void AdaptiveTimeStepping<TypeTag>::
updateTUNING(double max_next_tstep, const Tuning& tuning)
{
    restartFactor_ = tuning.TSFCNV;
    growthFactor_ = tuning.TFDIFF;
    maxGrowth_ = tuning.TSFMAX;
    maxTimeStep_ = tuning.TSMAXZ;
    updateNEXTSTEP(max_next_tstep);
    timestepAfterEvent_ = tuning.TMAXWC;
}

template<class TypeTag>
void AdaptiveTimeStepping<TypeTag>::
updateNEXTSTEP(double max_next_tstep)
{
     // \Note Only update next suggested step if TSINIT was explicitly set in TUNING or NEXTSTEP is active.
    if (max_next_tstep > 0) {
        suggestedNextTimestep_ = max_next_tstep;
    }
}

template<class TypeTag>
template<class Serializer>
void AdaptiveTimeStepping<TypeTag>::
serializeOp(Serializer& serializer)
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

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectHardcoded()
{
    return serializationTestObject_<HardcodedTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectPID()
{
    return serializationTestObject_<PIDTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectPIDIt()
{
    return serializationTestObject_<PIDAndIterationCountTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectSimple()
{
    return serializationTestObject_<SimpleIterationCountTimeStepControl>();
}

template<class TypeTag>
bool
AdaptiveTimeStepping<TypeTag>::
operator==(const AdaptiveTimeStepping<TypeTag>& rhs) const
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

    return result
        && this->restartFactor_ == rhs.restartFactor_
        && this->growthFactor_ == rhs.growthFactor_
        && this->maxGrowth_ == rhs.maxGrowth_
        && this->maxTimeStep_ == rhs.maxTimeStep_
        && this->minTimeStep_ == rhs.minTimeStep_
        && this->ignoreConvergenceFailure_ == rhs.ignoreConvergenceFailure_
        && this->solverRestartMax_== rhs.solverRestartMax_
        && this->solverVerbose_ == rhs.solverVerbose_
        && this->fullTimestepInitially_ == rhs.fullTimestepInitially_
        && this->timestepAfterEvent_ == rhs.timestepAfterEvent_
        && this->useNewtonIteration_ == rhs.useNewtonIteration_
        && this->minTimeStepBeforeShuttingProblematicWells_ ==
               rhs.minTimeStepBeforeShuttingProblematicWells_;
}

template<class TypeTag>
template<class Controller>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObject_()
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

template<class TypeTag>
void AdaptiveTimeStepping<TypeTag>::
init_(const UnitSystem& unitSystem)
{
    std::tie(timeStepControlType_,
             timeStepControl_,
             useNewtonIteration_) = detail::createController(unitSystem);

    // make sure growth factor is something reasonable
    if (growthFactor_ < 1.0) {
        OPM_THROW(std::runtime_error,
                  "Growth factor cannot be less than 1.");
    }
}

} // namespace Opm

#endif
