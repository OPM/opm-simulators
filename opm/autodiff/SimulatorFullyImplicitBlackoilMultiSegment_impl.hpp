/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS
  Copyright 2015 Andreas Lauser

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

    template <class GridT>
    auto SimulatorFullyImplicitBlackoilMultiSegment<GridT>::
    createSolver(const WellModel& well_model)
        -> std::unique_ptr<Solver>
    {
        typedef typename Traits::Model Model;

        auto model = std::unique_ptr<Model>(new Model(model_param_,
                                                      grid_,
                                                      props_,
                                                      geo_,
                                                      rock_comp_props_,
                                                      well_model,
                                                      solver_,
                                                      eclipse_state_,
                                                      schedule_,
                                                      summary_config_,
                                                      has_disgas_,
                                                      has_vapoil_,
                                                      terminal_output_));

        if (!Base::threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(Base::threshold_pressures_by_face_);
        }

        return std::unique_ptr<ThisType::Solver>(new Solver(Base::solver_param_, std::move(model)));

    }

    template <class GridT>
    SimulatorReport SimulatorFullyImplicitBlackoilMultiSegment<GridT>::run(SimulatorTimer& timer,
                                                                           ReservoirState& state)
    {
        WellState prev_well_state;

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_writer_.outputDirectory() + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        // adaptive time stepping
        const auto& events = schedule_->getEvents();
        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        if( param_.getDefault("timestep.adaptive", true ) )
        {
            adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_, terminal_output_ ) );
        }

        std::string restorefilename = param_.getDefault("restorefile", std::string("") );
        if( ! restorefilename.empty() )
        {
            // -1 means that we'll take the last report step that was written
            const int desiredRestoreStep = param_.getDefault("restorestep", int(-1) );
            output_writer_.restore( timer,
                                    state,
                                    prev_well_state,
                                    restorefilename,
                                    desiredRestoreStep );
        }

        unsigned int totalNonlinearIterations = 0;
        unsigned int totalLinearIterations = 0;
        DynamicListEconLimited dynamic_list_econ_limited;

        bool ooip_computed = false;
        std::vector<int> fipnum_global = eclipse_state_->get3DProperties().getIntGridProperty("FIPNUM").getData();
        //Get compressed cell fipnum.
        std::vector<int> fipnum(AutoDiffGrid::numCells(grid_));
        if (fipnum_global.empty()) {
            std::fill(fipnum.begin(), fipnum.end(), 0);
        } else {
            for (size_t c = 0; c < fipnum.size(); ++c) {
                fipnum[c] = fipnum_global[AutoDiffGrid::globalCell(grid_)[c]];
            }
        }
        std::vector<std::vector<double> > OOIP;

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            step_timer.start();
            if ( terminal_output_ )
            {
                timer.report(std::cout);
            }

            // Create wells and well state.
            WellsManager wells_manager(*eclipse_state_,
                                       *schedule_,
                                       timer.currentStepNum(),
                                       Opm::UgGridHelpers::numCells(grid_),
                                       Opm::UgGridHelpers::globalCell(grid_),
                                       Opm::UgGridHelpers::cartDims(grid_),
                                       Opm::UgGridHelpers::dimensions(grid_),
                                       Opm::UgGridHelpers::cell2Faces(grid_),
                                       Opm::UgGridHelpers::beginFaceCentroids(grid_),
                                       dynamic_list_econ_limited,
                                       is_parallel_run_,
                                       // We need to pass the optionaly arguments
                                       // as we get the following error otherwise
                                       // with c++ (Debian 4.9.2-10) 4.9.2 and -std=c++11
                                       // converting to ‘const std::unordered_set<std::basic_string<char> >’ from initializer list would use explicit constructor
                                       Base::defunct_well_names_);
            const Wells* wells = wells_manager.c_wells();
            WellState well_state;
            // well_state.init(wells, state, prev_well_state);

            const auto wells_ecl = schedule_->getWells(timer.currentStepNum());
            const int current_time_step = timer.currentStepNum();

            const WellModel well_model(wells, &(wells_manager.wellCollection()), wells_ecl, current_time_step);

            well_state.init(well_model, state, prev_well_state, wells);

            // give the polymer and surfactant simulators the chance to do their stuff
            Base::asImpl().handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // write the inital state at the report stage
            if (timer.initialStep()) {
                // No per cell data is written for initial step, but will be
                // for subsequent steps, when we have started simulating
                output_writer_.writeTimeStepWithoutCellProperties( timer, state, well_state, {}, {} );
            }

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            Base::asImpl().computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            auto solver = createSolver(well_model);

            // Compute orignal FIP;
            if (!ooip_computed) {
                OOIP = solver->computeFluidInPlace(state, fipnum);
                Base::FIPUnitConvert(eclipse_state_->getUnits(), OOIP);
                ooip_computed = true;
            }

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are too large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            if( adaptiveTimeStepping ) {
                bool event = events.hasEvent(ScheduleEvents::NEW_WELL, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::INJECTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE, timer.currentStepNum());
                adaptiveTimeStepping->step( timer, *solver, state, well_state, event, output_writer_);
            }
            else {
                // solve for complete report step
                solver->step(timer, state, well_state);
            }

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of nonlinear and linear Iterations
            totalNonlinearIterations += solver->nonlinearIterations();
            totalLinearIterations += solver->linearIterations();

            // Report timing.
            const double st = solver_timer.secsSinceStart();

            // Compute current FIP.
            std::vector<std::vector<double> > COIP;
            COIP = solver->computeFluidInPlace(state, fipnum);
            std::vector<double> OOIP_totals = Base::FIPTotals(OOIP, state);
            std::vector<double> COIP_totals = Base::FIPTotals(COIP, state);

            //Convert to correct units
            Base::FIPUnitConvert(eclipse_state_->getUnits(), COIP);
            Base::FIPUnitConvert(eclipse_state_->getUnits(), OOIP_totals);
            Base::FIPUnitConvert(eclipse_state_->getUnits(), COIP_totals);

            if ( terminal_output_ )
            {
                Base::outputFluidInPlace(OOIP_totals, COIP_totals,eclipse_state_->getUnits(), 0);
                for (size_t reg = 0; reg < OOIP.size(); ++reg) {
                    Base::outputFluidInPlace(OOIP[reg], COIP[reg], eclipse_state_->getUnits(), reg+1);
                }
            }

            if ( terminal_output_ )
            {
                std::cout << "Fully implicit solver took: " << st << " seconds." << std::endl;
            }

            stime += st;
            if ( output_writer_.output() ) {
                SimulatorReport step_report;
                step_report.solver_time = st;
                step_report.total_time =  step_timer.secsSinceStart();
                step_report.reportParam(tstep_os);
            }

            // Increment timer, remember well state.
            ++timer;


            // write simulation state at the report stage
            const auto& physicalModel = solver->model();
            output_writer_.writeTimeStep( timer, state, well_state, physicalModel );

            prev_well_state = well_state;
        }

        // Stop timer and create timing report
        total_timer.stop();
        SimulatorReport report;
        report.total_time = total_timer.secsSinceStart();
        report.solver_time = stime;
        report.total_newton_iterations = totalNonlinearIterations;
        report.total_linear_iterations = totalLinearIterations;
        return report;
    }

} // namespace Opm
