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
    createSolver(const Wells* wells, std::vector<WellMultiSegmentConstPtr>& wells_multisegment)
        -> std::unique_ptr<Solver>
    {
        typedef typename Traits::Model Model;

        auto model = std::unique_ptr<Model>(new Model(model_param_,
                                                      grid_,
                                                      props_,
                                                      geo_,
                                                      rock_comp_props_,
                                                      wells,
                                                      solver_,
                                                      eclipse_state_,
                                                      has_disgas_,
                                                      has_vapoil_,
                                                      terminal_output_,
                                                      wells_multisegment));

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
        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        if( param_.getDefault("timestep.adaptive", true ) )
        {
            adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_, terminal_output_ ) );
        }

        // init output writer
        output_writer_.writeInit( timer );

        std::string restorefilename = param_.getDefault("restorefile", std::string("") );
        if( ! restorefilename.empty() )
        {
            // -1 means that we'll take the last report step that was written
            const int desiredRestoreStep = param_.getDefault("restorestep", int(-1) );
            output_writer_.restore( timer, state, prev_well_state, restorefilename, desiredRestoreStep );
        }

        unsigned int totalNonlinearIterations = 0;
        unsigned int totalLinearIterations = 0;

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            step_timer.start();
            if ( terminal_output_ )
            {
                timer.report(std::cout);
            }

            // Create wells and well state.
            WellsManager wells_manager(eclipse_state_,
                                       timer.currentStepNum(),
                                       Opm::UgGridHelpers::numCells(grid_),
                                       Opm::UgGridHelpers::globalCell(grid_),
                                       Opm::UgGridHelpers::cartDims(grid_),
                                       Opm::UgGridHelpers::dimensions(grid_),
                                       Opm::UgGridHelpers::cell2Faces(grid_),
                                       Opm::UgGridHelpers::beginFaceCentroids(grid_),
                                       props_.permeability(),
                                       is_parallel_run_);
            const Wells* wells = wells_manager.c_wells();
            WellState well_state;
            // well_state.init(wells, state, prev_well_state);

            const std::vector<WellConstPtr>& wells_ecl = eclipse_state_->getSchedule()->getWells(timer.currentStepNum());
            std::vector<WellMultiSegmentConstPtr> wells_multisegment;
            wells_multisegment.reserve(wells_ecl.size());
            for (size_t i = 0; i < wells_ecl.size(); ++i) {
                if (wells_ecl[i]->getStatus(timer.currentStepNum()) == WellCommon::SHUT) {
                    continue;
                }
                wells_multisegment.push_back(std::make_shared<WellMultiSegment>(wells_ecl[i], timer.currentStepNum(), wells));
            }

#if 0
            // for DEBUGGING OUTPUT
            std::cout << " the number of the wells from EclipseState " << wells_ecl.size() << std::endl;
            for (size_t i = 0; i < wells_ecl.size(); ++i) {
                std::cout << " well name " << wells_ecl[i]->name() << std::endl;
                std::cout << " segment wells " << wells_ecl[i]->isMultiSegment() << std::endl;
            }

            std::cout << " output all the well information for debugging " << std::endl;

            int nw = wells_multisegment.size();
            for (int w = 0; w < nw; ++w) {
                WellMultiSegmentConstPtr well = wells_multisegment[w];
                std::cout << " well name " << well->name() << std::endl;
                std::cout << " is mutli-segmented ? : " << well->isMultiSegmented() << std::endl;
                std::cout << " number of segments : " << well->numberOfSegments() << std::endl;
                std::cout << " number of perforations : " << well->numberOfPerforations() << std::endl;

                std::cout << " beginning outputing segment informations " << std::endl;
                for (int s = 0; s < well->numberOfSegments(); ++s) {
                    std::cout << "    segment number : " << s << std::endl;
                    int n_perf_segment = well->segmentPerforations()[s].size();
                    std::cout << "    number of perforations for this segment " << n_perf_segment << std::endl;
                    std::cout << "    the depth of the segment " << well->segmentDepth()[s] << std::endl;
                    std::cout << "    the length of the segment " << well->segmentLength()[s] << std::endl;
                    std::cout << "    the volume of the segment " << well->segmentVolume()[s] << std::endl;
                    std::cout << "    the roughness of the segment " << well->segmentRoughness()[s] << std::endl;
                    std::cout << "    the cross area of the segment " << well->segmentCrossArea()[s] << std::endl;
                    std::cout << "    its outletSegment " << well->outletSegment()[s] << std::endl;
                    std::cout << "    the number of the inlet segments " << well->inletSegments()[s].size() << std::endl;
                    std::cout << "    its inlet segments are  ";
                    for (int inlet = 0; inlet < well->inletSegments()[s].size(); ++inlet) {
                        std::cout << well->inletSegments()[s][inlet] << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "    its perforations infromations " << std::endl;
                    for (int perf = 0; perf < n_perf_segment; ++perf) {
                        int perf_number = well->segmentPerforations()[s][perf];
                        std::cout << "      perforation " << perf_number;
                        std::cout << "   peforation depth " << well->perfDepth()[perf_number] << std::endl;;
                    }
                    std::cout << std::endl;
                }
                std::cout << " output all the mapping informations " << std::endl;
                std::cout << " matrix s2p " << std::endl;
                std::cout << well->wellOps().s2p << std::endl;


                std::cout << " maxtrix p2s " << std::endl;
                std::cout << well->wellOps().p2s << std::endl;

                std::cout << " matrix p2s_average " << std::endl;
                std::cout << well->wellOps().p2s_average << std::endl;

                std::cout << " maxtrix s2s_gather " << std::endl;
                std::cout << well->wellOps().s2s_gather << std::endl;

                std::cout << " maxtrix p2s_gather " << std::endl;
                std::cout << well->wellOps().p2s_gather << std::endl;

                std::cout << " s2s_inlets " << std::endl;
                std::cout << well->wellOps().s2s_inlets << std::endl;

                std::cout << " s2s_outlet " << std::endl;
                std::cout << well->wellOps().s2s_outlet << std::endl;

                std::cout << " output well information for well " << well->name() << " done!!!! " << std::endl;
            }
            std::cin.ignore();
#endif
            // DEBUGGING OUTPUT is DONE


            well_state.init(wells_multisegment, state, prev_well_state);

            // give the polymer and surfactant simulators the chance to do their stuff
            Base::asImpl().handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // write simulation state at the report stage
            output_writer_.writeTimeStep( timer, state, well_state );

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            Base::asImpl().computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            auto solver = createSolver(wells, wells_multisegment);

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are too large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            if( adaptiveTimeStepping ) {
                adaptiveTimeStepping->step( timer, *solver, state, well_state,  output_writer_ );
            }
            else {
                // solve for complete report step
                solver->step(timer.currentStepLength(), state, well_state);
            }

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of nonlinear and linear Iterations
            totalNonlinearIterations += solver->nonlinearIterations();
            totalLinearIterations += solver->linearIterations();

            // Report timing.
            const double st = solver_timer.secsSinceStart();

            if ( terminal_output_ )
            {
                std::cout << "Fully implicit solver took: " << st << " seconds." << std::endl;
            }

            stime += st;
            if ( output_writer_.output() ) {
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
        output_writer_.writeTimeStep( timer, state, prev_well_state );

        // Stop timer and create timing report
        total_timer.stop();
        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        report.total_newton_iterations = totalNonlinearIterations;
        report.total_linear_iterations = totalLinearIterations;
        return report;
    }

} // namespace Opm
