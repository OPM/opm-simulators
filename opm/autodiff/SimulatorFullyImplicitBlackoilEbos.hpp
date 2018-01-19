/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
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

#include <opm/autodiff/BlackoilOutputEbos.hpp>
#include <opm/autodiff/IterationReport.hpp>
#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilWellModel.hpp>
#include <opm/autodiff/moduleVersion.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/core/utility/StopWatch.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <dune/common/unused.hh>

namespace Opm {


/// a simulator for the blackoil model
template<class TypeTag>
class SimulatorFullyImplicitBlackoilEbos
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables)  PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector)    SolutionVector ;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

    typedef WellStateFullyImplicitBlackoil WellState;
    typedef BlackoilState ReservoirState;
    typedef BlackoilOutputEbos<TypeTag> OutputWriter;
    typedef BlackoilModelEbos<TypeTag> Model;
    typedef BlackoilModelParameters ModelParameters;
    typedef NonlinearSolver<Model> Solver;
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
    /// \param[in] has_disgas    true for dissolved gas option
    /// \param[in] has_vapoil    true for vaporized oil option
    /// \param[in] eclipse_state the object which represents an internalized ECL deck
    /// \param[in] output_writer
    /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
    SimulatorFullyImplicitBlackoilEbos(Simulator& ebosSimulator,
                                       const ParameterGroup& param,
                                       NewtonIterationBlackoilInterface& linsolver,
                                       const bool has_disgas,
                                       const bool has_vapoil,
                                       OutputWriter& output_writer)
        : ebosSimulator_(ebosSimulator),
          param_(param),
          model_param_(param),
          solver_param_(param),
          solver_(linsolver),
          phaseUsage_(phaseUsageFromDeck(eclState())),
          has_disgas_(has_disgas),
          has_vapoil_(has_vapoil),
          terminal_output_(param.getDefault("output_terminal", true)),
          output_writer_(output_writer),
          is_parallel_run_( false )
    {
#if HAVE_MPI
        if ( solver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
        {
            const ParallelISTLInformation& info =
                boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation());
            // Only rank 0 does print to std::cout
            terminal_output_ = terminal_output_ && ( info.communicator().rank() == 0 );
            is_parallel_run_ = ( info.communicator().size() > 1 );
        }
#endif
    }

    /// Run the simulation.
    /// This will run succesive timesteps until timer.done() is true. It will
    /// modify the reservoir and well states.
    /// \param[in,out] timer       governs the requested reporting timesteps
    /// \param[in,out] state       state of reservoir: pressure, fluxes
    /// \return                    simulation report, with timing data
    SimulatorReport run(SimulatorTimer& timer)
    {

        ReservoirState dummy_state(0,0,0);

        WellState prev_well_state;

        ExtraData extra;

        failureReport_ = SimulatorReport();

        if (output_writer_.isRestart()) {
            // This is a restart, populate WellState
            ReservoirState stateInit(Opm::UgGridHelpers::numCells(grid()),
                                     Opm::UgGridHelpers::numFaces(grid()),
                                     phaseUsage_.num_phases);
            output_writer_.initFromRestartFile(phaseUsage_, grid(), stateInit, prev_well_state, extra);
        }

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();

        // adaptive time stepping
        const auto& events = schedule().getEvents();
        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        const bool useTUNING = param_.getDefault("use_TUNING", false);
        if( param_.getDefault("timestep.adaptive", true ) )
        {
            if (useTUNING) {
                adaptiveTimeStepping.reset( new AdaptiveTimeStepping( schedule().getTuning(), timer.currentStepNum(), param_, terminal_output_ ) );
            } else {
                adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_, terminal_output_ ) );
            }

            if (output_writer_.isRestart()) {
                if (extra.suggested_step > 0.0) {
                    adaptiveTimeStepping->setSuggestedNextStep(extra.suggested_step);
                }
            }
        }

        SimulatorReport report;
        SimulatorReport stepReport;

        WellModel well_model(ebosSimulator_, model_param_, terminal_output_);
        if (output_writer_.isRestart()) {
            well_model.setRestartWellState(prev_well_state); // Neccessary for perfect restarts
        }

        WellState wellStateDummy; //not used. Only passed to make the old interfaces happy

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            if ( terminal_output_ )
            {
                std::ostringstream ss;
                timer.report(ss);
                OpmLog::debug(ss.str());
            }

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            well_model.beginReportStep(timer.currentStepNum());

            auto solver = createSolver(well_model);

            // write the inital state at the report stage
            if (timer.initialStep()) {
                Dune::Timer perfTimer;
                perfTimer.start();

                // No per cell data is written for initial step, but will be
                // for subsequent steps, when we have started simulating
                output_writer_.writeTimeStep( timer, dummy_state, well_model.wellState(), solver->model() );

                report.output_write_time += perfTimer.stop();
            }

            if( terminal_output_ )
            {
                std::ostringstream step_msg;
                boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
                step_msg.imbue(std::locale(std::locale::classic(), facet));
                step_msg << "\nReport step " << std::setw(2) <<timer.currentStepNum()
                         << "/" << timer.numSteps()
                         << " at day " << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day)
                         << "/" << (double)unit::convert::to(timer.totalTime(), unit::day)
                         << ", date = " << timer.currentDateTime();
                OpmLog::info(step_msg.str());
            }

            solver->model().beginReportStep();

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are too large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            if( adaptiveTimeStepping ) {
                if (useTUNING) {
                    if(events.hasEvent(ScheduleEvents::TUNING_CHANGE,timer.currentStepNum())) {
                        adaptiveTimeStepping->updateTUNING(schedule().getTuning(), timer.currentStepNum());
                    }
                }

                bool event = events.hasEvent(ScheduleEvents::NEW_WELL, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::INJECTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE, timer.currentStepNum());
                stepReport = adaptiveTimeStepping->step( timer, *solver, dummy_state, wellStateDummy, event, output_writer_, nullptr );
                report += stepReport;
                failureReport_ += adaptiveTimeStepping->failureReport();
            }
            else {
                // solve for complete report step
                stepReport = solver->step(timer, dummy_state, wellStateDummy);
                report += stepReport;
                failureReport_ += solver->failureReport();

                if( terminal_output_ )
                {
                    //stepReport.briefReport();
                    std::ostringstream iter_msg;
                    iter_msg << "Stepsize " << (double)unit::convert::to(timer.currentStepLength(), unit::day);
                    if (solver->wellIterations() != 0) {
                        iter_msg << " days well iterations = " << solver->wellIterations() << ", ";
                    }
                    iter_msg << "non-linear iterations = " << solver->nonlinearIterations()
                             << ", total linear iterations = " << solver->linearIterations()
                             << "\n";
                    OpmLog::info(iter_msg.str());
                }
            }

            solver->model().endReportStep();
            well_model.endReportStep();

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // update timing.
            report.solver_time += solver_timer.secsSinceStart();

            // Increment timer, remember well state.
            ++timer;


            if (terminal_output_ )
            {
                if (!timer.initialStep()) {
                    const std::string version = moduleVersionName();
                    outputTimestampFIP(timer, version);
                }
            }

            // write simulation state at the report stage
            Dune::Timer perfTimer;
            perfTimer.start();
            const double nextstep = adaptiveTimeStepping ? adaptiveTimeStepping->suggestedNextStep() : -1.0;            
            output_writer_.writeTimeStep( timer, dummy_state, well_model.wellState(), solver->model(), false, nextstep, report);
            report.output_write_time += perfTimer.stop();

            if (terminal_output_ )
            {
                std::string msg =
                    "Time step took " + std::to_string(solver_timer.secsSinceStart()) + " seconds; "
                    "total solver time " + std::to_string(report.solver_time) + " seconds.";
                OpmLog::debug(msg);
            }

        }

        // Stop timer and create timing report
        total_timer.stop();
        report.total_time = total_timer.secsSinceStart();
        report.converged = true;
        return report;
    }

    /** \brief Returns the simulator report for the failed substeps of the simulation.
     */
    const SimulatorReport& failureReport() const { return failureReport_; };

    const Grid& grid() const
    { return ebosSimulator_.gridManager().grid(); }

protected:

    std::unique_ptr<Solver> createSolver(WellModel& well_model)
    {
        auto model = std::unique_ptr<Model>(new Model(ebosSimulator_,
                                                      model_param_,
                                                      well_model,
                                                      solver_,
                                                      terminal_output_));

        return std::unique_ptr<Solver>(new Solver(solver_param_, std::move(model)));
    }

    void outputTimestampFIP(const SimulatorTimer& timer, const std::string version)
    {
        std::ostringstream ss;
        boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d %b %Y");
        ss.imbue(std::locale(std::locale::classic(), facet));
        ss << "\n                              **************************************************************************\n"
        << "  Balance  at" << std::setw(10) << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day) << "  Days"
        << " *" << std::setw(30) << eclState().getTitle() << "                                          *\n"
        << "  Report " << std::setw(4) << timer.reportStepNum() << "    " << timer.currentDateTime()
        << "  *                                             Flow  version " << std::setw(11) << version << "  *\n"
        << "                              **************************************************************************\n";
        OpmLog::note(ss.str());
    }

    const EclipseState& eclState() const
    { return ebosSimulator_.gridManager().eclState(); }


    const Schedule& schedule() const
    { return ebosSimulator_.gridManager().schedule(); }


    // Data.
    Simulator& ebosSimulator_;

    typedef typename Solver::SolverParameters SolverParameters;

    SimulatorReport failureReport_;

    const ParameterGroup param_;
    ModelParameters model_param_;
    SolverParameters solver_param_;

    // Observed objects.
    NewtonIterationBlackoilInterface& solver_;
    PhaseUsage phaseUsage_;
    // Misc. data
    const bool has_disgas_;
    const bool has_vapoil_;
    bool       terminal_output_;
    // output_writer
    OutputWriter& output_writer_;

    // Whether this a parallel simulation or not
    bool is_parallel_run_;

};

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
