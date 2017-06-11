/*
  Copyright 2014 IRIS AS

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
#ifndef OPM_ADAPTIVETIMESTEPPING_IMPL_HEADER_INCLUDED
#define OPM_ADAPTIVETIMESTEPPING_IMPL_HEADER_INCLUDED

#include <iostream>
#include <string>
#include <utility>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <dune/istl/istlexception.hh>
#include <dune/istl/ilu.hh> // For MatrixBlockException
#include <opm/parser/eclipse/EclipseState/Schedule/Tuning.hpp>

namespace Opm {

    namespace detail
    {
        template <class Solver, class State>
        class SolutionTimeErrorSolverWrapper : public RelativeChangeInterface
        {
            const Solver& solver_;
            const State&  previous_;
            const State&  current_;
        public:
            SolutionTimeErrorSolverWrapper( const Solver& solver,
                                            const State&  previous,
                                            const State&  current )
              : solver_( solver ),
                previous_( previous ),
                current_( current )
            {}

            /// return || u^n+1 - u^n || / || u^n+1 ||
            double relativeChange() const
            {
                return solver_.model().relativeChange( previous_, current_ );
            }
        };

        template<class E>
        void logException(const E& exception, bool verbose)
        {
            if( verbose )
            {
                std::string message;
                message = "Caught Exception: ";
                message += exception.what();
                OpmLog::debug(message);
            }
        }
    }

    // AdaptiveTimeStepping
    //---------------------

    AdaptiveTimeStepping::AdaptiveTimeStepping( const Tuning& tuning,
                                                size_t time_step,
                                                const ParameterGroup& param,
                                                const bool terminal_output )
        : timeStepControl_()
        , restart_factor_( tuning.getTSFCNV(time_step) )
        , growth_factor_(tuning.getTFDIFF(time_step) )
        , max_growth_( tuning.getTSFMAX(time_step) )
          // default is 1 year, convert to seconds
        , max_time_step_( tuning.getTSMAXZ(time_step) )
        , solver_restart_max_( param.getDefault("solver.restart", int(10) ) )
        , solver_verbose_( param.getDefault("solver.verbose", bool(true) ) && terminal_output )
        , timestep_verbose_( param.getDefault("timestep.verbose", bool(true) ) && terminal_output )
        , suggested_next_timestep_( tuning.getTSINIT(time_step) )
        , full_timestep_initially_( param.getDefault("full_timestep_initially", bool(false) ) )
        , timestep_after_event_( tuning.getTMAXWC(time_step))
        , use_newton_iteration_(false)
    {
        init(param);

    }

    AdaptiveTimeStepping::AdaptiveTimeStepping( const ParameterGroup& param,
                                                const bool terminal_output )
        : timeStepControl_()
        , restart_factor_( param.getDefault("solver.restartfactor", double(0.33) ) )
        , growth_factor_( param.getDefault("solver.growthfactor", double(2) ) )
        , max_growth_( param.getDefault("timestep.control.maxgrowth", double(3.0) ) )
          // default is 1 year, convert to seconds
        , max_time_step_( unit::convert::from(param.getDefault("timestep.max_timestep_in_days", 365.0 ), unit::day) )
        , solver_restart_max_( param.getDefault("solver.restart", int(10) ) )
        , solver_verbose_( param.getDefault("solver.verbose", bool(true) ) && terminal_output )
        , timestep_verbose_( param.getDefault("timestep.verbose", bool(true) ) && terminal_output )
        , suggested_next_timestep_( unit::convert::from(param.getDefault("timestep.initial_timestep_in_days", -1.0 ), unit::day) )
        , full_timestep_initially_( param.getDefault("full_timestep_initially", bool(false) ) )
        , timestep_after_event_( unit::convert::from(param.getDefault("timestep.timestep_in_days_after_event", -1.0 ), unit::day))
        , use_newton_iteration_(false)
    {
        init(param);
    }

    void AdaptiveTimeStepping::
    init(const ParameterGroup& param)
    {
        // valid are "pid" and "pid+iteration"
        std::string control = param.getDefault("timestep.control", std::string("pid") );
        // iterations is the accumulation of all linear iterations over all newton steops per time step
        const int defaultTargetIterations = 30;
        const int defaultTargetNewtonIterations = 8;

        const double tol = param.getDefault("timestep.control.tol", double(1e-1) );
        if( control == "pid" ) {
            timeStepControl_ = TimeStepControlType( new PIDTimeStepControl( tol ) );
        }
        else if ( control == "pid+iteration" )
        {
            const int iterations   = param.getDefault("timestep.control.targetiteration", defaultTargetIterations );
            timeStepControl_ = TimeStepControlType( new PIDAndIterationCountTimeStepControl( iterations, tol ) );
        }
        else if ( control == "pid+newtoniteration" )
        {
            const int iterations   = param.getDefault("timestep.control.targetiteration", defaultTargetNewtonIterations );
            timeStepControl_ = TimeStepControlType( new PIDAndIterationCountTimeStepControl( iterations, tol ) );
            use_newton_iteration_ = true;
        }
        else if ( control == "iterationcount" )
        {
            const int iterations    = param.getDefault("timestep.control.targetiteration", defaultTargetIterations );
            const double decayrate  = param.getDefault("timestep.control.decayrate",  double(0.75) );
            const double growthrate = param.getDefault("timestep.control.growthrate", double(1.25) );
            timeStepControl_ = TimeStepControlType( new SimpleIterationCountTimeStepControl( iterations, decayrate, growthrate ) );
        } else if ( control == "hardcoded") {
            const std::string filename    = param.getDefault("timestep.control.filename", std::string("timesteps"));
            timeStepControl_ = TimeStepControlType( new HardcodedTimeStepControl( filename ) );

        }
        else
            OPM_THROW(std::runtime_error,"Unsupported time step control selected "<< control );

        // make sure growth factor is something reasonable
        assert( growth_factor_ >= 1.0 );
    }



    template <class Solver, class State, class WellState>
    SimulatorReport AdaptiveTimeStepping::
    step( const SimulatorTimer& simulatorTimer, Solver& solver, State& state, WellState& well_state, const bool event )
    {
        return stepImpl( simulatorTimer, solver, state, well_state, event, nullptr, nullptr );
    }

    template <class Solver, class State, class WellState, class Output>
    SimulatorReport AdaptiveTimeStepping::
    step( const SimulatorTimer& simulatorTimer,
          Solver& solver, State& state, WellState& well_state,
          const bool event,
          Output& outputWriter,
          const std::vector<int>* fipnum)
    {
        return stepImpl( simulatorTimer, solver, state, well_state, event, &outputWriter, fipnum );
    }


    // implementation of the step method
    template <class Solver, class State, class WState, class Output >
    SimulatorReport AdaptiveTimeStepping::
    stepImpl( const SimulatorTimer& simulatorTimer,
              Solver& solver, State& state, WState& well_state,
              const bool event,
              Output* outputWriter,
              const std::vector<int>* fipnum)
    {
        SimulatorReport report;
        const double timestep = simulatorTimer.currentStepLength();

        // init last time step as a fraction of the given time step
        if( suggested_next_timestep_ < 0 ) {
            suggested_next_timestep_ = restart_factor_ * timestep;
        }

        if (full_timestep_initially_) {
            suggested_next_timestep_ = timestep;
        }

        // use seperate time step after event
        if (event && timestep_after_event_ > 0) {
            suggested_next_timestep_ = timestep_after_event_;
        }

        // create adaptive step timer with previously used sub step size
        AdaptiveSimulatorTimer substepTimer( simulatorTimer, suggested_next_timestep_, max_time_step_ );

        // copy states in case solver has to be restarted (to be revised)
        State  last_state( state );
        WState last_well_state( well_state );

        // reset the statistics for the failed substeps
        failureReport_ = SimulatorReport();

        // counter for solver restarts
        int restarts = 0;

        // sub step time loop
        while( ! substepTimer.done() )
        {
            // get current delta t
            const double dt = substepTimer.currentStepLength() ;
            if( timestep_verbose_ )
            {
                std::ostringstream ss;
                ss <<"  Substep " << substepTimer.currentStepNum() << ", stepsize "
                   << unit::convert::to(substepTimer.currentStepLength(), unit::day) << " days.";
                OpmLog::info(ss.str());
            }

            SimulatorReport substepReport;
            std::string cause_of_failure = "";
            try {
                substepReport = solver.step( substepTimer, state, well_state);
                report += substepReport;

                if( solver_verbose_ ) {
                    // report number of linear iterations
                    OpmLog::note("Overall linear iterations used: " + std::to_string(substepReport.total_linear_iterations));
                }
            }
            catch (const Opm::TooManyIterations& e) {
                substepReport += solver.failureReport();
                cause_of_failure = "Solver convergence failure - Iteration limit reached";

                detail::logException(e, solver_verbose_);
                // since linearIterations is < 0 this will restart the solver
            }
            catch (const Opm::LinearSolverProblem& e) {
                substepReport += solver.failureReport();

                detail::logException(e, solver_verbose_);
                // since linearIterations is < 0 this will restart the solver
            }
            catch (const Opm::NumericalProblem& e) {
                substepReport += solver.failureReport();

                detail::logException(e, solver_verbose_);
                // since linearIterations is < 0 this will restart the solver
            }
            catch (const std::runtime_error& e) {
                substepReport += solver.failureReport();

                detail::logException(e, solver_verbose_);
                // also catch linear solver not converged
            }
            catch (const Dune::ISTLError& e) {
                substepReport += solver.failureReport();

                detail::logException(e, solver_verbose_);
                // also catch errors in ISTL AMG that occur when time step is too large
            }
            catch (const Dune::MatrixBlockError& e) {
                substepReport += solver.failureReport();

                detail::logException(e, solver_verbose_);
                // this can be thrown by ISTL's ILU0 in block mode, yet is not an ISTLError
            }

            if( substepReport.converged )
            {
                // advance by current dt
                ++substepTimer;

                // create object to compute the time error, simply forwards the call to the model
                detail::SolutionTimeErrorSolverWrapper< Solver, State >
                    relativeChange( solver, last_state, state );

                // compute new time step estimate
                const int iterations = use_newton_iteration_ ? substepReport.total_newton_iterations
                                                             : substepReport.total_linear_iterations;
                double dtEstimate = timeStepControl_->computeTimeStepSize( dt, iterations, relativeChange,
                                                                           substepTimer.simulationTimeElapsed());

                // limit the growth of the timestep size by the growth factor
                dtEstimate = std::min( dtEstimate, double(max_growth_ * dt) );

                // further restrict time step size growth after convergence problems
                if( restarts > 0 ) {
                    dtEstimate = std::min( growth_factor_ * dt, dtEstimate );
                    // solver converged, reset restarts counter
                    restarts = 0;
                }

                if( timestep_verbose_ )
                {
                    std::ostringstream ss;
                    ss << "    Substep summary: ";
                    if (substepReport.total_well_iterations != 0) {
                        ss << "well its = " << std::setw(2) << substepReport.total_well_iterations << ", ";
                    }
                    ss << "newton its = " << std::setw(2) << substepReport.total_newton_iterations << ", "
                       << "linearizations = "  << std::setw(2) << substepReport.total_linearizations
                       << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << substepReport.assemble_time << " sec), "
                       << "linear its = " << std::setw(3) << substepReport.total_linear_iterations
                       << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << substepReport.linear_solve_time << " sec)";
                    OpmLog::info(ss.str());
                }

                // write data if outputWriter was provided
                // if the time step is done we do not need
                // to write it as this will be done by the simulator
                // anyway.
                if( outputWriter && !substepTimer.done() ) {
                    if (fipnum) {
                        solver.computeFluidInPlace(state, *fipnum);
                    }
                    Opm::time::StopWatch perfTimer;
                    perfTimer.start();
                    bool substep = true;
                    const auto& physicalModel = solver.model();
                    outputWriter->writeTimeStep( substepTimer, state, well_state, physicalModel, substep);
                    report.output_write_time += perfTimer.secsSinceStart();
                }

                // set new time step length
                substepTimer.provideTimeStepEstimate( dtEstimate );

                // update states
                last_state      = state ;
                last_well_state = well_state;

                report.converged = substepTimer.done();
                substepTimer.setLastStepFailed(false);

            }
            else // in case of no convergence (linearIterations < 0)
            {
                substepTimer.setLastStepFailed(true);

                failureReport_ += substepReport;

                // increase restart counter
                if( restarts >= solver_restart_max_ ) {
                    const auto msg = std::string("Solver failed to converge after cutting timestep ")
                        + std::to_string(restarts) + " times.";
                    if (solver_verbose_) {
                        OpmLog::error(msg);
                    }
                    OPM_THROW_NOLOG(Opm::NumericalProblem, msg);
                }

                const double newTimeStep = restart_factor_ * dt;
                // we need to revise this
                substepTimer.provideTimeStepEstimate( newTimeStep );
                if( solver_verbose_ ) {
                    std::string msg;
                    msg = cause_of_failure + "\n         Timestep chopped to "
                        + std::to_string(unit::convert::to( substepTimer.currentStepLength(), unit::day )) + " days\n";
                    OpmLog::problem(msg);
                }
                // reset states
                state      = last_state;
                well_state = last_well_state;

                ++restarts;
            }
        }


        // store estimated time step for next reportStep
        suggested_next_timestep_ = substepTimer.currentStepLength();
        if( timestep_verbose_ )
        {
            std::ostringstream ss;
            substepTimer.report(ss);
            ss << "Suggested next step size = " << unit::convert::to( suggested_next_timestep_, unit::day ) << " (days)" << std::endl;
            OpmLog::note(ss.str());
        }

        if( ! std::isfinite( suggested_next_timestep_ ) ) { // check for NaN
            suggested_next_timestep_ = timestep;
        }
        return report;
    }
}

#endif
