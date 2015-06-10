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

#include <algorithm>

namespace Opm
{

    template <class Implementation>
    SimulatorBase<Implementation>::SimulatorBase(const parameter::ParameterGroup& param,
                                                 const Grid& grid,
                                                 const DerivedGeology& geo,
                                                 BlackoilPropsAdInterface& props,
                                                 const RockCompressibility* rock_comp_props,
                                                 NewtonIterationBlackoilInterface& linsolver,
                                                 const double* gravity,
                                                 const bool has_disgas,
                                                 const bool has_vapoil,
                                                 std::shared_ptr<EclipseState> eclipse_state,
                                                 OutputWriter& output_writer,
                                                 const std::vector<double>& threshold_pressures_by_face)
        : param_(param),
          model_param_(param),
          solver_param_(param),
          grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          gravity_(gravity),
          geo_(geo),
          solver_(linsolver),
          has_disgas_(has_disgas),
          has_vapoil_(has_vapoil),
          terminal_output_(param.getDefault("output_terminal", true)),
          eclipse_state_(eclipse_state),
          output_writer_(output_writer),
          rateConverter_(props_, std::vector<int>(AutoDiffGrid::numCells(grid_), 0)),
          threshold_pressures_by_face_(threshold_pressures_by_face)
    {
        // Misc init.
        const int num_cells = AutoDiffGrid::numCells(grid);
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
#if HAVE_MPI
        if ( terminal_output_ ) {
            if ( solver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation());
                // Only rank 0 does print to std::cout
                terminal_output_ = ( info.communicator().rank() == 0 );
                is_parallel_run_ = ( info.communicator().size() > 1 );
            }
        }
#endif
    }

    template <class Implementation>
    SimulatorReport SimulatorBase<Implementation>::run(SimulatorTimer& timer,
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
            adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_, solver_.parallelInformation() ) );
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

        unsigned int totalNewtonIterations = 0;
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
            well_state.init(wells, state, prev_well_state);

            // give the polymer and surfactant simulators the chance to do their stuff
            asImpl().handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // write simulation state at the report stage
            output_writer_.writeTimeStep( timer, state, well_state );

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            asImpl().computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            auto solver = asImpl().createSolver(wells);

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

            // accumulate the number of Newton and Linear Iterations
            totalNewtonIterations += solver->newtonIterations();
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
        report.total_newton_iterations = totalNewtonIterations;
        report.total_linear_iterations = totalLinearIterations;
        return report;
    }

    namespace SimFIBODetails {
        typedef std::unordered_map<std::string, WellConstPtr> WellMap;

        inline WellMap
        mapWells(const std::vector<WellConstPtr>& wells)
        {
            WellMap wmap;

            for (std::vector<WellConstPtr>::const_iterator
                     w = wells.begin(), e = wells.end();
                 w != e; ++w)
            {
                wmap.insert(std::make_pair((*w)->name(), *w));
            }

            return wmap;
        }

        inline int
        resv_control(const WellControls* ctrl)
        {
            int i, n = well_controls_get_num(ctrl);

            bool match = false;
            for (i = 0; (! match) && (i < n); ++i) {
                match = well_controls_iget_type(ctrl, i) == RESERVOIR_RATE;
            }

            if (! match) { i = 0; }

            return i - 1; // -1 if no match, undo final "++" otherwise
        }

        inline bool
        is_resv(const Wells& wells,
                const int    w)
        {
            return (0 <= resv_control(wells.ctrls[w]));
        }

        inline bool
        is_resv(const WellMap&     wmap,
                const std::string& name,
                const std::size_t  step)
        {
            bool match = false;

            WellMap::const_iterator i = wmap.find(name);

            if (i != wmap.end()) {
                WellConstPtr wp = i->second;

                match = (wp->isProducer(step) &&
                         wp->getProductionProperties(step)
                         .hasProductionControl(WellProducer::RESV))
                    ||  (wp->isInjector(step) &&
                         wp->getInjectionProperties(step)
                         .hasInjectionControl(WellInjector::RESV));
            }

            return match;
        }

        inline std::vector<int>
        resvWells(const Wells*      wells,
                  const std::size_t step,
                  const WellMap&    wmap)
        {
            std::vector<int> resv_wells;
            if( wells )
            {
                for (int w = 0, nw = wells->number_of_wells; w < nw; ++w) {
                    if (is_resv(*wells, w) ||
                        ((wells->name[w] != 0) &&
                         is_resv(wmap, wells->name[w], step)))
                    {
                        resv_wells.push_back(w);
                    }
                }
            }

            return resv_wells;
        }

        inline void
        historyRates(const PhaseUsage&               pu,
                     const WellProductionProperties& p,
                     std::vector<double>&            rates)
        {
            assert (! p.predictionMode);
            assert (rates.size() ==
                    std::vector<double>::size_type(pu.num_phases));

            if (pu.phase_used[ BlackoilPhases::Aqua ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Aqua ];

                rates[i] = p.WaterRate;
            }

            if (pu.phase_used[ BlackoilPhases::Liquid ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Liquid ];

                rates[i] = p.OilRate;
            }

            if (pu.phase_used[ BlackoilPhases::Vapour ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Vapour ];

                rates[i] = p.GasRate;
            }
        }
    } // namespace SimFIBODetails

    template <class Implementation>
    void SimulatorBase<Implementation>::handleAdditionalWellInflow(SimulatorTimer& /* timer */,
                                                                   WellsManager& /* wells_manager */,
                                                                   WellState& /* well_state */,
                                                                   const Wells* /* wells */)
    { }

    template <class Implementation>
    auto SimulatorBase<Implementation>::createSolver(const Wells* wells)
        -> std::unique_ptr<Solver>
    {
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
                                                      terminal_output_));

        if (!threshold_pressures_by_face_.empty()) {
            model->setThresholdPressures(threshold_pressures_by_face_);
        }

        return std::unique_ptr<Solver>(new Solver(solver_param_, std::move(model)));
    }

    template <class Implementation>
    void SimulatorBase<Implementation>::computeRESV(const std::size_t               step,
                                                    const Wells*                    wells,
                                                    const BlackoilState&            x,
                                                    WellState& xw)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const std::vector<WellConstPtr>& w_ecl = eclipse_state_->getSchedule()->getWells(step);
        const WellMap& wmap = SimFIBODetails::mapWells(w_ecl);

        const std::vector<int>& resv_wells = SimFIBODetails::resvWells(wells, step, wmap);

        if (! resv_wells.empty()) {
            const PhaseUsage&                    pu = props_.phaseUsage();
            const std::vector<double>::size_type np = props_.numPhases();

            rateConverter_.defineState(x);

            std::vector<double> distr (np);
            std::vector<double> hrates(np);
            std::vector<double> prates(np);

            for (std::vector<int>::const_iterator
                     rp = resv_wells.begin(), e = resv_wells.end();
                 rp != e; ++rp)
            {
                WellControls* ctrl = wells->ctrls[*rp];
                const bool is_producer = wells->type[*rp] == PRODUCER;

                // RESV control mode, all wells
                {
                    const int rctrl = SimFIBODetails::resv_control(ctrl);

                    if (0 <= rctrl) {
                        const std::vector<double>::size_type off = (*rp) * np;

                        if (is_producer) {
                            // Convert to positive rates to avoid issues
                            // in coefficient calculations.
                            std::transform(xw.wellRates().begin() + (off + 0*np),
                                           xw.wellRates().begin() + (off + 1*np),
                                           prates.begin(), std::negate<double>());
                        } else {
                            std::copy(xw.wellRates().begin() + (off + 0*np),
                                      xw.wellRates().begin() + (off + 1*np),
                                      prates.begin());
                        }

                        const int fipreg = 0; // Hack.  Ignore FIP regions.
                        rateConverter_.calcCoeff(prates, fipreg, distr);

                        well_controls_iset_distr(ctrl, rctrl, & distr[0]);
                    }
                }

                // RESV control, WCONHIST wells.  A bit of duplicate
                // work, regrettably.
                if (is_producer && wells->name[*rp] != 0) {
                    WellMap::const_iterator i = wmap.find(wells->name[*rp]);

                    if (i != wmap.end()) {
                        WellConstPtr wp = i->second;

                        const WellProductionProperties& p =
                            wp->getProductionProperties(step);

                        if (! p.predictionMode) {
                            // History matching (WCONHIST/RESV)
                            SimFIBODetails::historyRates(pu, p, hrates);

                            const int fipreg = 0; // Hack.  Ignore FIP regions.
                            rateConverter_.calcCoeff(hrates, fipreg, distr);

                            // WCONHIST/RESV target is sum of all
                            // observed phase rates translated to
                            // reservoir conditions.  Recall sign
                            // convention: Negative for producers.
                            const double target =
                                - std::inner_product(distr.begin(), distr.end(),
                                                     hrates.begin(), 0.0);

                            well_controls_clear(ctrl);
                            well_controls_assert_number_of_phases(ctrl, int(np));

                            const int ok_resv =
                                well_controls_add_new(RESERVOIR_RATE, target,
                                                      & distr[0], ctrl);

                            // For WCONHIST/RESV the BHP limit is set to 1 atm.
                            // TODO: Make it possible to modify the BHP limit using
                            // the WELTARG keyword
                            const int ok_bhp =
                                well_controls_add_new(BHP, unit::convert::from(1.0, unit::atm),
                                                      NULL, ctrl);

                            if (ok_resv != 0 && ok_bhp != 0) {
                                xw.currentControls()[*rp] = 0;
                                well_controls_set_current(ctrl, 0);
                            }
                        }
                    }
                }
            }
        }
    }
} // namespace Opm
