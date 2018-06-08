/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP

#include <iostream>
#include <utility>

#include <opm/grid/utility/StopWatch.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>

namespace Opm {
    // AdaptiveTimeStepping
    //---------------------

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
        AdaptiveTimeSteppingEbos(const ParameterGroup& param,
                                 const bool terminalOutput = true)
            : timeStepControl_()
            , restartFactor_( param.getDefault("solver.restartfactor", double(0.33)))
            , growthFactor_( param.getDefault("solver.growthfactor", double(2)))
            , maxGrowth_( param.getDefault("timestep.control.maxgrowth", double(3.0)))
              // default is 1 year, convert to seconds
            , maxTimeStep_( unit::convert::from(param.getDefault("timestep.max_timestep_in_days", 365.0), unit::day))
            , solverRestartMax_( param.getDefault("solver.restart", int(10)))
            , solverVerbose_( param.getDefault("solver.verbose", bool(true)) && terminalOutput)
            , timestepVerbose_( param.getDefault("timestep.verbose", bool(true)) && terminalOutput)
            , suggestedNextTimestep_( unit::convert::from(param.getDefault("timestep.initial_timestep_in_days", 1.0), unit::day))
            , fullTimestepInitially_( param.getDefault("full_timestep_initially", bool(false)))
            , timestepAfterEvent_( unit::convert::from(param.getDefault("timestep.timestep_in_days_after_event", -1.0), unit::day))
            , useNewtonIteration_(false)
        {
            init_(param);
        }



        //! \brief contructor taking parameter object
        //! \param tuning Pointer to ecl TUNING keyword
        //! \param timeStep current report step
        //! \param param The parameter object
        AdaptiveTimeSteppingEbos(const Tuning& tuning,
                                 size_t timeStep,
                                 const ParameterGroup& param,
                                 const bool terminalOutput = true)
            : timeStepControl_()
            , restartFactor_(tuning.getTSFCNV(timeStep))
            , growthFactor_(tuning.getTFDIFF(timeStep))
            , maxGrowth_(tuning.getTSFMAX(timeStep))
              // default is 1 year, convert to seconds
            , maxTimeStep_(tuning.getTSMAXZ(timeStep))
            , solverRestartMax_(param.getDefault("solver.restart", int(10)))
            , solverVerbose_(param.getDefault("solver.verbose", bool(true)) && terminalOutput)
            , timestepVerbose_(param.getDefault("timestep.verbose", bool(true)) && terminalOutput)
            , suggestedNextTimestep_(tuning.getTSINIT(timeStep))
            , fullTimestepInitially_(param.getDefault("full_timestep_initially", bool(false)))
            , timestepAfterEvent_(tuning.getTMAXWC(timeStep))
            , useNewtonIteration_(false)
        {
            init_(param);
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
            auto phaseUsage = Opm::phaseUsageFromDeck(ebosSimulator.vanguard().eclState());

            // create adaptive step timer with previously used sub step size
            AdaptiveSimulatorTimer substepTimer(simulatorTimer, suggestedNextTimestep_, maxTimeStep_);

            // reset the statistics for the failed substeps
            failureReport_ = SimulatorReport();

            // counter for solver restarts
            int restarts = 0;

            // sub step time loop
            while (!substepTimer.done()) {
                // get current delta t
                const double dt = substepTimer.currentStepLength() ;
                if (timestepVerbose_) {
                    std::ostringstream ss;
                    ss <<"\nTime step " << substepTimer.currentStepNum() << ", stepsize "
                       << unit::convert::to(substepTimer.currentStepLength(), unit::day) << " days.";
                    OpmLog::info(ss.str());
                }

                SimulatorReport substepReport;
                std::string causeOfFailure = "";
                try {
                    substepReport = solver.step(substepTimer);
                    report += substepReport;

                    if (solverVerbose_) {
                        // report number of linear iterations
                        OpmLog::debug("Overall linear iterations used: " + std::to_string(substepReport.total_linear_iterations));
                    }
                }
                catch (const Opm::TooManyIterations& e) {
                    substepReport += solver.failureReport();
                    causeOfFailure = "Solver convergence failure - Iteration limit reached";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const Opm::LinearSolverProblem& e) {
                    substepReport += solver.failureReport();
                    causeOfFailure = "Linear solver convergence failure";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const Opm::NumericalIssue& e) {
                    substepReport += solver.failureReport();
                    causeOfFailure = "Solver convergence failure - Numerical problem encountered";

                    logException_(e, solverVerbose_);
                    // since linearIterations is < 0 this will restart the solver
                }
                catch (const std::runtime_error& e) {
                    substepReport += solver.failureReport();

                    logException_(e, solverVerbose_);
                    // also catch linear solver not converged
                }
                catch (const Dune::ISTLError& e) {
                    substepReport += solver.failureReport();

                    logException_(e, solverVerbose_);
                    // also catch errors in ISTL AMG that occur when time step is too large
                }
                catch (const Dune::MatrixBlockError& e) {
                    substepReport += solver.failureReport();

                    logException_(e, solverVerbose_);
                    // this can be thrown by ISTL's ILU0 in block mode, yet is not an ISTLError
                }

                if (substepReport.converged) {
                    // advance by current dt
                    ++substepTimer;

                    // create object to compute the time error, simply forwards the call to the model
                    SolutionTimeErrorSolverWrapperEbos<Solver> relativeChange(solver);

                    // compute new time step estimate
                    const int iterations = useNewtonIteration_ ? substepReport.total_newton_iterations
                        : substepReport.total_linear_iterations;
                    double dtEstimate = timeStepControl_->computeTimeStepSize(dt, iterations, relativeChange,
                                                                               substepTimer.simulationTimeElapsed());

                    // limit the growth of the timestep size by the growth factor
                    dtEstimate = std::min(dtEstimate, double(maxGrowth_ * dt));

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
                        if (fipnum) {
                            solver.computeFluidInPlace(*fipnum);
                        }
                        Opm::time::StopWatch perfTimer;
                        perfTimer.start();

                        // The writeOutput expects a local data::solution vector and a local data::well vector.
                        auto localWellData = solver.model().wellModel().wellState().report(phaseUsage, Opm::UgGridHelpers::globalCell(ebosSimulator.vanguard().grid()));
                        ebosProblem.writeOutput(localWellData,
                                                substepTimer.simulationTimeElapsed(),
                                                /*isSubstep=*/true,
                                                substepReport.total_time,
                                                /*nextStepSize=*/-1.0);

                        report.output_write_time += perfTimer.secsSinceStart();
                    }

                    // set new time step length
                    substepTimer.provideTimeStepEstimate(dtEstimate);

                    report.converged = substepTimer.done();
                    substepTimer.setLastStepFailed(false);

                }
                else { // in case of no convergence (linearIterations < 0)
                    substepTimer.setLastStepFailed(true);

                    failureReport_ += substepReport;

                    // increase restart counter
                    if (restarts >= solverRestartMax_) {
                        const auto msg = std::string("Solver failed to converge after cutting timestep ")
                            + std::to_string(restarts) + " times.";
                        if (solverVerbose_) {
                            OpmLog::error(msg);
                        }
                        OPM_THROW_NOLOG(Opm::NumericalIssue, msg);
                    }

                    const double newTimeStep = restartFactor_ * dt;
                    // we need to revise this
                    substepTimer.provideTimeStepEstimate(newTimeStep);
                    if (solverVerbose_) {
                        std::string msg;
                        msg = causeOfFailure + "\nTimestep chopped to "
                            + std::to_string(unit::convert::to(substepTimer.currentStepLength(), unit::day)) + " days\n";
                        OpmLog::problem(msg);
                    }

                    ++restarts;
                }
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
        const SimulatorReport& failureReport() const
        { return failureReport_; };

        double suggestedNextStep() const
        { return suggestedNextTimestep_; }

        void setSuggestedNextStep(const double x)
        { suggestedNextTimestep_ = x; }

        void updateTUNING(const Tuning& tuning, size_t timeStep)
        {
            restartFactor_ = tuning.getTSFCNV(timeStep);
            growthFactor_ = tuning.getTFDIFF(timeStep);
            maxGrowth_ = tuning.getTSFMAX(timeStep);
            maxTimeStep_ = tuning.getTSMAXZ(timeStep);
            suggestedNextTimestep_ = tuning.getTSINIT(timeStep);
            timestepAfterEvent_ = tuning.getTMAXWC(timeStep);
        }


    protected:
        void init_(const ParameterGroup& param)
        {
            // valid are "pid" and "pid+iteration"
            std::string control = param.getDefault("timestep.control", std::string("pid"));
            // iterations is the accumulation of all linear iterations over all newton steops per time step
            const int defaultTargetIterations = 30;
            const int defaultTargetNewtonIterations = 8;

            const double tol = param.getDefault("timestep.control.tol", double(1e-1));
            if (control == "pid") {
                timeStepControl_ = TimeStepControlType(new PIDTimeStepControl(tol));
            }
            else if (control == "pid+iteration") {
                const int iterations = param.getDefault("timestep.control.targetiteration", defaultTargetIterations);
                timeStepControl_ = TimeStepControlType(new PIDAndIterationCountTimeStepControl(iterations, tol));
            }
            else if (control == "pid+newtoniteration") {
                const int iterations = param.getDefault("timestep.control.targetiteration", defaultTargetNewtonIterations);
                timeStepControl_ = TimeStepControlType(new PIDAndIterationCountTimeStepControl(iterations, tol));
                useNewtonIteration_ = true;
            }
            else if (control == "iterationcount") {
                const int iterations = param.getDefault("timestep.control.targetiteration", defaultTargetIterations);
                const double decayrate = param.getDefault("timestep.control.decayrate",  double(0.75));
                const double growthrate = param.getDefault("timestep.control.growthrate", double(1.25));
                timeStepControl_ = TimeStepControlType(new SimpleIterationCountTimeStepControl(iterations, decayrate, growthrate));
            }
            else if (control == "hardcoded") {
                const std::string filename = param.getDefault("timestep.control.filename", std::string("timesteps"));
                timeStepControl_ = TimeStepControlType(new HardcodedTimeStepControl(filename));

            }
            else
                OPM_THROW(std::runtime_error,"Unsupported time step control selected "<< control);

            // make sure growth factor is something reasonable
            assert(growthFactor_ >= 1.0);
        }


        typedef std::unique_ptr<TimeStepControlInterface> TimeStepControlType;

        SimulatorReport failureReport_;       //!< statistics for the failed substeps of the last timestep
        TimeStepControlType timeStepControl_; //!< time step control object
        double restartFactor_;               //!< factor to multiply time step with when solver fails to converge
        double growthFactor_;                //!< factor to multiply time step when solver recovered from failed convergence
        double maxGrowth_;                   //!< factor that limits the maximum growth of a time step
        double maxTimeStep_;                //!< maximal allowed time step size
        const int solverRestartMax_;        //!< how many restart of solver are allowed
        const bool solverVerbose_;           //!< solver verbosity
        const bool timestepVerbose_;         //!< timestep verbosity
        double suggestedNextTimestep_;      //!< suggested size of next timestep
        bool fullTimestepInitially_;        //!< beginning with the size of the time step from data file
        double timestepAfterEvent_;         //!< suggested size of timestep after an event
        bool useNewtonIteration_;           //!< use newton iteration count for adaptive time step control
    };
}

#endif // OPM_ADAPTIVE_TIME_STEPPING_EBOS_HPP
