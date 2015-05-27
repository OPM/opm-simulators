/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS
  Copyright 2014 STATOIL ASA.

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

namespace Opm
{
    template<class T>
    SimulatorFullyImplicitBlackoilPolymer<T>::
    SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
                                          const Grid& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const PolymerPropsAd& polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity,
                                          const bool has_disgas,
                                          const bool has_vapoil,
                                          const bool has_polymer,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
                                          const std::vector<double>& threshold_pressures_by_face)
    : BaseType(param,
               grid,
               geo,
               props,
               rock_comp_props,
               linsolver,
               gravity,
               has_disgas,
               has_vapoil,
               eclipse_state,
               output_writer,
               threshold_pressures_by_face)
        , polymer_props_(polymer_props)
        , has_polymer_(has_polymer)
        , deck_(deck)
    {
        // Misc init.
        const int num_cells = AutoDiffGrid::numCells(grid);
        BaseType::allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            BaseType::allcells_[cell] = cell;
        }
#if HAVE_MPI
        if ( BaseType::terminal_output_ ) {
            if ( BaseType::solver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(BaseType::solver_.parallelInformation());
                // Only rank 0 does print to std::cout
                BaseType::terminal_output_= (info.communicator().rank()==0);
            }
        }
#endif
    }




    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoilPolymer<T>::run(SimulatorTimer& timer,
                                                                  PolymerBlackoilState& state)
    {
        WellStateFullyImplicitBlackoilPolymer prev_well_state;

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = BaseType::output_writer_.outputDirectory() + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        typedef T Grid;
        typedef BlackoilPolymerModel<Grid> Model;
        typedef typename Model::ModelParameters ModelParams;
        ModelParams modelParams( BaseType::param_ );
        typedef NewtonSolver<Model> Solver;
        typedef typename Solver::SolverParameters SolverParams;
        SolverParams solverParams( BaseType::param_ );

        //adaptive time stepping
        //        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        //        if( BaseType::param_.getDefault("timestep.adaptive", bool(false) ) )
        //        {
        //            adaptiveTimeStepping.reset( new AdaptiveTimeStepping( BaseType::param_ ) );
        //        }

        // init output writer
        BaseType::output_writer_.writeInit( timer );

        std::string restorefilename = BaseType::param_.getDefault("restorefile", std::string("") );
        if( ! restorefilename.empty() )
        {
            // -1 means that we'll take the last report step that was written
            const int desiredRestoreStep = BaseType::param_.getDefault("restorestep", int(-1) );
            BaseType::output_writer_.restore( timer, state.blackoilState(), prev_well_state, restorefilename, desiredRestoreStep );
        }

        unsigned int totalNewtonIterations = 0;
        unsigned int totalLinearIterations = 0;

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            step_timer.start();
            if ( BaseType::terminal_output_ )
            {
                timer.report(std::cout);
            }

            // Create wells and well state.
            WellsManager wells_manager(BaseType::eclipse_state_,
                                       timer.currentStepNum(),
                                       Opm::UgGridHelpers::numCells(BaseType::grid_),
                                       Opm::UgGridHelpers::globalCell(BaseType::grid_),
                                       Opm::UgGridHelpers::cartDims(BaseType::grid_),
                                       Opm::UgGridHelpers::dimensions(BaseType::grid_),
                                       Opm::UgGridHelpers::cell2Faces(BaseType::grid_),
                                       Opm::UgGridHelpers::beginFaceCentroids(BaseType::grid_),
                                       BaseType::props_.permeability());
            const Wells* wells = wells_manager.c_wells();
            WellStateFullyImplicitBlackoilPolymer well_state;
            well_state.init(wells, state.blackoilState(), prev_well_state);

            // compute polymer inflow
            std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
            if (deck_->hasKeyword("WPOLYMER")) {
                if (wells_manager.c_wells() == 0) {
                    OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
                }
                polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, BaseType::eclipse_state_, *wells, Opm::UgGridHelpers::numCells(BaseType::grid_), timer.currentStepNum()));
            } else {
                polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                                1.0*Opm::unit::day,
                                                                0.0));
            }
            std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(BaseType::grid_));
            polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(), 
                                                timer.simulationTimeElapsed() + timer.currentStepLength(),
                                                polymer_inflow_c);
            well_state.polymerInflow() = polymer_inflow_c;

            // write simulation state at the report stage
            BaseType::output_writer_.writeTimeStep( timer, state.blackoilState(), well_state );

            // Max oil saturation (for VPPARS), hysteresis update.
            BaseType::props_.updateSatOilMax(state.saturation());
            BaseType::props_.updateSatHyst(state.saturation(), BaseType::allcells_);

            // Compute reservoir volumes for RESV controls.
            BaseType::computeRESV(timer.currentStepNum(), wells, state.blackoilState(), well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            Model model(modelParams, BaseType::grid_, BaseType::props_, BaseType::geo_, BaseType::rock_comp_props_, polymer_props_, wells, BaseType::solver_, BaseType::has_disgas_, BaseType::has_vapoil_, has_polymer_, BaseType::terminal_output_);
            if (!BaseType::threshold_pressures_by_face_.empty()) {
                model.setThresholdPressures(BaseType::threshold_pressures_by_face_);
            }
            Solver solver(solverParams, model);

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are to large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            //            if( adaptiveTimeStepping ) {
            //                adaptiveTimeStepping->step( timer, solver, state, well_state,  BaseType::output_writer_ );
            //            } else {
                // solve for complete report step
            solver.step(timer.currentStepLength(), state, well_state);
                //            }

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of Newton and Linear Iterations
            totalNewtonIterations += solver.newtonIterations();
            totalLinearIterations += solver.linearIterations();

            // Report timing.
            const double st = solver_timer.secsSinceStart();

            if ( BaseType::terminal_output_ )
            {
                std::cout << "Fully implicit solver took: " << st << " seconds." << std::endl;
            }

            stime += st;
            if ( BaseType::output_writer_.output() ) {
                SimulatorReport step_report;
                step_report.pressure_time = st;
                step_report.total_time =  step_timer.secsSinceStart();
                step_report.reportParam(tstep_os);
            }

            // Increment timer, remember well state.
            ++timer;
            prev_well_state = well_state;
        }

        // Write final simulation state.
        BaseType::output_writer_.writeTimeStep( timer, state.blackoilState(), prev_well_state );

        // Stop timer and create timing report
        total_timer.stop();
        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        report.total_newton_iterations = totalNewtonIterations;
        report.total_linear_iterations = totalLinearIterations;
        return report;
    }
} // namespace Opm
