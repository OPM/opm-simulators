/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014-2016 IRIS AS
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

#include <utility>
#include <functional>
#include <algorithm>
#include <locale>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/core/utility/initHydroCarbonState.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/wells/DynamicListEconLimited.hpp>

namespace Opm
{

    template <class Implementation>
    SimulatorBase<Implementation>::SimulatorBase(const parameter::ParameterGroup& param,
                                                 const Grid& grid,
                                                 DerivedGeology& geo,
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
          threshold_pressures_by_face_(threshold_pressures_by_face),
          is_parallel_run_( false )
    {
        // Misc init.
        const int num_cells = AutoDiffGrid::numCells(grid);
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
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

    template <class Implementation>
    SimulatorReport SimulatorBase<Implementation>::run(SimulatorTimer& timer,
                                                       ReservoirState& state)
    {
        WellState prev_well_state;


        if (output_writer_.isRestart()) {
            // This is a restart, populate WellState and ReservoirState state objects from restart file
            output_writer_.initFromRestartFile(props_.phaseUsage(), props_.permeability(), grid_, state, prev_well_state);
            initHydroCarbonState(state, props_.phaseUsage(), Opm::UgGridHelpers::numCells(grid_), has_disgas_, has_vapoil_);
        }

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_writer_.outputDirectory() + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        const auto& schedule = eclipse_state_->getSchedule();
        const auto& events = schedule->getEvents();

        // adaptive time stepping
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

        unsigned int totalLinearizations = 0;
        unsigned int totalNonlinearIterations = 0;
        unsigned int totalLinearIterations = 0;
        bool is_well_potentials_computed = param_.getDefault("compute_well_potentials", false );
        std::vector<double> well_potentials;
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
        std::vector<V> OOIP;
        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            step_timer.start();
            if ( terminal_output_ )
            {
                std::ostringstream ss;
                timer.report(ss);
                OpmLog::note(ss.str());
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
                                       dynamic_list_econ_limited,
                                       is_parallel_run_,
                                       well_potentials);
            const Wells* wells = wells_manager.c_wells();
            WellState well_state;
            well_state.init(wells, state, prev_well_state);

            // give the polymer and surfactant simulators the chance to do their stuff
            asImpl().handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // write the inital state at the report stage
            if (timer.initialStep()) {
                // No per cell data is written for initial step, but will be
                // for subsequent steps, when we have started simulating
                output_writer_.writeTimeStepWithoutCellProperties( timer, state, well_state );
            }

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            asImpl().computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            const WellModel well_model(wells);

            std::unique_ptr<Solver> solver = asImpl().createSolver(well_model);

            // Compute orignal FIP;
            if (!ooip_computed) {
                OOIP = solver->computeFluidInPlace(state, fipnum);
                FIPUnitConvert(eclipse_state_->getUnits(), OOIP);
                ooip_computed = true;
            }

            if( terminal_output_ )
            {
                std::ostringstream step_msg;
                boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
                step_msg.imbue(std::locale(std::locale::classic(), facet));
                step_msg << "\nTime step " << std::setw(4) <<timer.currentStepNum()
                         << " at day " << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day)
                         << "/" << (double)unit::convert::to(timer.totalTime(), unit::day)
                         << ", date = " << timer.currentDateTime()
                         << "\n";
                OpmLog::info(step_msg.str());
            }

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are too large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            if( adaptiveTimeStepping ) {
                adaptiveTimeStepping->step( timer, *solver, state, well_state, output_writer_ );
            }
            else {
                // solve for complete report step
                solver->step(timer, state, well_state);

                if( terminal_output_ )
                {
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

            // update the derived geology (transmissibilities, pore volumes, etc) if the
            // has geology changed for the next report step
            const int nextTimeStepIdx = timer.currentStepNum() + 1;
            if (nextTimeStepIdx < timer.numSteps()
                && events.hasEvent(ScheduleEvents::GEO_MODIFIER, nextTimeStepIdx)) {
                // bring the contents of the keywords to the current state of the SCHEDULE
                // section
                //
                // TODO (?): handle the parallel case (maybe this works out of the box)
                DeckConstPtr miniDeck = schedule->getModifierDeck(nextTimeStepIdx);
                eclipse_state_->applyModifierDeck(*miniDeck);
                geo_.update(grid_, props_, eclipse_state_, gravity_);
            }

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of nonlinear and linear Iterations
            totalLinearizations += solver->linearizations();
            totalNonlinearIterations += solver->nonlinearIterations();
            totalLinearIterations += solver->linearIterations();

            // Report timing.
            const double st = solver_timer.secsSinceStart();
            // Compute current FIP.
            std::vector<V> COIP;
            COIP = solver->computeFluidInPlace(state, fipnum);
            FIPUnitConvert(eclipse_state_->getUnits(), COIP);
            V OOIP_totals = FIPTotals(OOIP, state);
            V COIP_totals = FIPTotals(COIP, state);
            outputFluidInPlace(OOIP_totals, COIP_totals,eclipse_state_->getUnits(), 0);
            for (size_t reg = 0; reg < OOIP.size(); ++reg) {
                outputFluidInPlace(OOIP[reg], COIP[reg], eclipse_state_->getUnits(), reg+1);
            }


            // accumulate total time
            stime += st;
            
            if ( terminal_output_ )
            {
                std::string msg;
                msg = "Fully implicit solver took: " + std::to_string(st) + " seconds. Total solver time taken: " + std::to_string(stime) + " seconds.";
                OpmLog::note(msg);
            }

            if ( output_writer_.output() ) {
                SimulatorReport step_report;
                step_report.pressure_time = st;
                step_report.total_time =  step_timer.secsSinceStart();
                step_report.reportParam(tstep_os);
            }

            // Increment timer, remember well state.
            ++timer;

            // write simulation state at the report stage
            const auto& physicalModel = solver->model();
            output_writer_.writeTimeStep( timer, state, well_state, physicalModel );

            prev_well_state = well_state;
            // The well potentials are only computed if they are needed
            // For now thay are only used to determine default guide rates for group controlled wells
            if ( is_well_potentials_computed ) {
                asImpl().computeWellPotentials(wells, well_state, well_potentials);
            }

            asImpl().updateListEconLimited(solver, eclipse_state_->getSchedule(), timer.currentStepNum(), wells,
                                           well_state, dynamic_list_econ_limited);
        }

        // Stop timer and create timing report
        total_timer.stop();
        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        report.total_linearizations = totalLinearizations;
        report.total_newton_iterations = totalNonlinearIterations;
        report.total_linear_iterations = totalLinearIterations;
        return report;
    }

    namespace SimFIBODetails {
        typedef std::unordered_map<std::string, const Well* > WellMap;

        inline WellMap
        mapWells(const std::vector< const Well* >& wells)
        {
            WellMap wmap;

            for (std::vector< const Well* >::const_iterator
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
                const Well* wp = i->second;

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
    auto SimulatorBase<Implementation>::createSolver(const WellModel& well_model)
        -> std::unique_ptr<Solver>
    {
        auto model = std::unique_ptr<Model>(new Model(model_param_,
                                                      grid_,
                                                      props_,
                                                      geo_,
                                                      rock_comp_props_,
                                                      well_model,
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
    void SimulatorBase<Implementation>::computeWellPotentials(const Wells* wells,
                                                              const WellState& xw,
                                                              std::vector<double>& well_potentials)
    {
        const int nw = wells->number_of_wells;
        const int np = wells->number_of_phases;
        well_potentials.clear();
        well_potentials.resize(nw*np,0.0);
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                for (int phase = 0; phase < np; ++phase) {
                    well_potentials[w*np + phase] += xw.wellPotentials()[perf*np + phase];
                }
            }
        }
    }

    template <class Implementation>
    void SimulatorBase<Implementation>::computeRESV(const std::size_t               step,
                                                    const Wells*                    wells,
                                                    const BlackoilState&            x,
                                                    WellState& xw)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const auto w_ecl = eclipse_state_->getSchedule()->getWells(step);
        const WellMap& wmap = SimFIBODetails::mapWells(w_ecl);

        const std::vector<int>& resv_wells = SimFIBODetails::resvWells(wells, step, wmap);

        const std::size_t number_resv_wells        = resv_wells.size();
        std::size_t       global_number_resv_wells = number_resv_wells;
#if HAVE_MPI
        if ( solver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
        {
            const auto& info =
                boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation());
            global_number_resv_wells = info.communicator().sum(global_number_resv_wells);
            if ( global_number_resv_wells )
            {
                // At least one process has resv wells. Therefore rate converter needs
                // to calculate averages over regions that might cross process
                // borders. This needs to be done by all processes and therefore
                // outside of the next if statement.
                rateConverter_.defineState(x, boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation()));
            }
        }
        else
#endif
        {
            if ( global_number_resv_wells )
            {
                rateConverter_.defineState(x);
            }
        }

        if (! resv_wells.empty()) {
            const PhaseUsage&                    pu = props_.phaseUsage();
            const std::vector<double>::size_type np = props_.numPhases();

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
                        const auto* wp = i->second;

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

                            static const double invalid_alq = -std::numeric_limits<double>::max();
                            static const int invalid_vfp = -std::numeric_limits<int>::max();

                            const int ok_resv =
                                well_controls_add_new(RESERVOIR_RATE, target,
                                                      invalid_alq, invalid_vfp,
                                                      & distr[0], ctrl);

                            // For WCONHIST the BHP limit is set to 1 atm.
                            // or a value specified using WELTARG
                            double bhp_limit = (p.BHPLimit > 0) ? p.BHPLimit : unit::convert::from(1.0, unit::atm);
                            const int ok_bhp =
                                well_controls_add_new(BHP, bhp_limit,
                                                      invalid_alq, invalid_vfp,
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

        if( wells )
        {
            for (int w = 0, nw = wells->number_of_wells; w < nw; ++w) {
                WellControls* ctrl = wells->ctrls[w];
                const bool is_producer = wells->type[w] == PRODUCER;
                if (!is_producer && wells->name[w] != 0) {
                    WellMap::const_iterator i = wmap.find(wells->name[w]);
                    if (i != wmap.end()) {
                        const auto* wp = i->second;
                        const WellInjectionProperties& injector = wp->getInjectionProperties(step);
                        if (!injector.predictionMode) {
                            //History matching WCONINJEH
                            static const double invalid_alq = -std::numeric_limits<double>::max();
                            static const int invalid_vfp = -std::numeric_limits<int>::max();
                            // For WCONINJEH the BHP limit is set to a large number
                            // or a value specified using WELTARG
                            double bhp_limit = (injector.BHPLimit > 0) ? injector.BHPLimit : std::numeric_limits<double>::max();
                            const int ok_bhp =
                                well_controls_add_new(BHP, bhp_limit,
                                                      invalid_alq, invalid_vfp,
                                                      NULL, ctrl);
                            if (!ok_bhp) {
                                OPM_THROW(std::runtime_error, "Failed to add well control.");
                            }
                        }
                    }
                }
            }
        }
    }


    template <class Implementation>
    void
    SimulatorBase<Implementation>::FIPUnitConvert(const UnitSystem& units,
                                                  std::vector<V>& fip)
    {
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            for (size_t i = 0; i < fip.size(); ++i) {
                fip[i][0] = unit::convert::to(fip[i][0], unit::stb);
                fip[i][1] = unit::convert::to(fip[i][1], unit::stb); 
                fip[i][2] = unit::convert::to(fip[i][2], 1000*unit::cubic(unit::feet));
                fip[i][3] = unit::convert::to(fip[i][3], 1000*unit::cubic(unit::feet));
                fip[i][4] = unit::convert::to(fip[i][4], unit::stb);
                fip[i][5] = unit::convert::to(fip[i][5], unit::stb);
                fip[i][6] = unit::convert::to(fip[i][6], unit::psia);
            }
        }
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            for (size_t i = 0; i < fip.size(); ++i) {
                fip[i][6] = unit::convert::to(fip[i][6], unit::barsa);
            }
        }
    }


    template <class Implementation>
    V
    SimulatorBase<Implementation>::FIPTotals(const std::vector<V>& fip, const ReservoirState& state)
    {
        V totals(V::Zero(7));
        for (int i = 0; i < 5; ++i) {
            for (size_t reg = 0; reg < fip.size(); ++reg) {
                totals[i] += fip[reg][i];
            }
        }
        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const int np = state.numPhases();
        const PhaseUsage& pu = props_.phaseUsage();
        const DataBlock s = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        const V so = pu.phase_used[BlackoilPhases::Liquid] ? V(s.col(BlackoilPhases::Liquid)) : V::Zero(nc);
        const V sg = pu.phase_used[BlackoilPhases::Vapour] ? V(s.col(BlackoilPhases::Vapour)) : V::Zero(nc);
        const V hydrocarbon = so + sg;
        const V p = Eigen::Map<const V>(& state.pressure()[0], nc);
        if ( ! is_parallel_run_ )
        {
            totals[5] = geo_.poreVolume().sum();
            totals[6] = unit::convert::to((p * geo_.poreVolume() * hydrocarbon).sum() / ((geo_.poreVolume() * hydrocarbon).sum()), unit::barsa);
        }
        else
        {
#if HAVE_MPI
            const auto & pinfo =
                boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation());
            auto operators = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<double>(),
                                             Opm::Reduction::makeGlobalSumFunctor<double>(),
                                             Opm::Reduction::makeGlobalSumFunctor<double>());
            auto pav_nom   = p * geo_.poreVolume() * hydrocarbon;
            auto pav_denom = geo_.poreVolume() * hydrocarbon;

            // using ref cref to prevent copying
            auto inputs = std::make_tuple(std::cref(geo_.poreVolume()),
                                          std::cref(pav_nom), std::cref(pav_denom));
            std::tuple<double, double, double> results(0.0, 0.0, 0.0);

            pinfo.computeReduction(inputs, operators, results);
            using std::get;
            totals[5] = get<0>(results);
            totals[6] = unit::convert::to(get<1>(results)/get<2>(results),
                                          unit::barsa);
#else
            // This should never happen!
            OPM_THROW(std::logic_error, "HAVE_MPI should be defined if we are running in parallel");
#endif
        }

        return totals;
    }



    template <class Implementation>
    void
    SimulatorBase<Implementation>::outputFluidInPlace(const V& oip, const V& cip, const UnitSystem& units, const int reg)
    {
        std::ostringstream ss;
        if (!reg) {
            ss << "                                                  ===================================================\n"
               << "                                                  :                   Field Totals                  :\n";
        } else {
            ss << "                                                  ===================================================\n"
               << "                                                  :        FIPNUM report region  "
               << std::setw(2) << reg << "                 :\n";
        }
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << "                                                  :      PAV  =" << std::setw(14) << cip[6] << " BARSA                 :\n"
               << std::fixed << std::setprecision(0)
               << "                                                  :      PORV =" << std::setw(14) << cip[5] << "   RM3                 :\n";
            if (!reg) {
                ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
                   << "                                                  : Porv volumes are taken at reference conditions  :\n";
            }
            ss << "                         :--------------- Oil    SM3 ---------------:-- Wat    SM3 --:--------------- Gas    SM3 ---------------:\n";
        }
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << "                                                  :      PAV  =" << std::setw(14) << cip[6] << "  PSIA                 :\n"
               << std::fixed << std::setprecision(0)
               << "                                                  :      PORV =" << std::setw(14) << cip[5] << "   RB                  :\n";
            if (!reg) {
                ss << "                                                  : Pressure is weighted by hydrocarbon pore voulme :\n"
                   << "                                                  : Pore volumes are taken at reference conditions  :\n";
            }
            ss << "                         :--------------- Oil    STB ---------------:-- Wat    STB --:--------------- Gas   MSCF ---------------:\n";
        }
        ss << "                         :      Liquid        Vapour        Total   :      Total     :      Free        Dissolved       Total   :" << "\n"
           << ":------------------------:------------------------------------------:----------------:------------------------------------------:" << "\n"
           << ":Currently   in place    :" << std::setw(14) << cip[1] << std::setw(14) << cip[4] << std::setw(14) << (cip[1]+cip[4]) << ":"
           << std::setw(13) << cip[0] << "   :" << std::setw(14) << (cip[2]) << std::setw(14) << cip[3] << std::setw(14) << (cip[2] + cip[3]) << ":\n"
           << ":------------------------:------------------------------------------:----------------:------------------------------------------:\n"
           << ":Originally  in place    :" << std::setw(14) << oip[1] << std::setw(14) << oip[4] << std::setw(14) << (oip[1]+oip[4]) << ":"
           << std::setw(13) << oip[0] << "   :" << std::setw(14) << oip[2] << std::setw(14) << oip[3] << std::setw(14) << (oip[2] + oip[3]) << ":\n"
           << ":========================:==========================================:================:==========================================:\n";
        OpmLog::note(ss.str());
    }


    template <class Implementation>
    void
    SimulatorBase<Implementation>::
    updateListEconLimited(const std::unique_ptr<Solver>& solver,
                          ScheduleConstPtr schedule,
                          const int current_step,
                          const Wells* wells,
                          const WellState& well_state,
                          DynamicListEconLimited& list_econ_limited) const
    {

        solver->model().wellModel().updateListEconLimited(schedule, current_step, wells,
                                                          well_state, list_econ_limited);
    }

} // namespace Opm
