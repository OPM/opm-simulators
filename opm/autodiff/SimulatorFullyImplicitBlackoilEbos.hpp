/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILEBOS_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILEBOS_HEADER_INCLUDED

//#include <opm/autodiff/SimulatorBase.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutputEbos.hpp>
#include <opm/autodiff/IterationReport.hpp>
#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoilDense.hpp>
#include <opm/autodiff/StandardWellsDense.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/SimFIBODetails.hpp>

#include <opm/core/simulator/AdaptiveTimeStepping.hpp>
#include <opm/core/utility/initHydroCarbonState.hpp>
#include <opm/core/utility/StopWatch.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

namespace Opm {

class SimulatorFullyImplicitBlackoilEbos;
//class StandardWellsDense<FluidSystem>;

/// a simulator for the blackoil model
class SimulatorFullyImplicitBlackoilEbos
{
public:
    typedef typename TTAG(EclFlowProblem) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;

    typedef WellStateFullyImplicitBlackoilDense WellState;
    typedef BlackoilState ReservoirState;
    typedef BlackoilOutputWriterEbos OutputWriter;
    typedef BlackoilModelEbos Model;
    typedef BlackoilModelParameters ModelParameters;
    typedef NonlinearSolver<Model> Solver;
    typedef StandardWellsDense<FluidSystem, BlackoilIndices> WellModel;


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
    /// \param[in] geo           derived geological properties
    /// \param[in] props         fluid and rock properties
    /// \param[in] rock_comp_props     if non-null, rock compressibility properties
    /// \param[in] linsolver     linear solver
    /// \param[in] gravity       if non-null, gravity vector
    /// \param[in] has_disgas    true for dissolved gas option
    /// \param[in] has_vapoil    true for vaporized oil option
    /// \param[in] eclipse_state the object which represents an internalized ECL deck
    /// \param[in] output_writer
    /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
    SimulatorFullyImplicitBlackoilEbos(Simulator& ebosSimulator,
                                       const parameter::ParameterGroup& param,
                                       DerivedGeology& geo,
                                       BlackoilPropsAdInterface& props,
                                       const RockCompressibility* rock_comp_props,
                                       NewtonIterationBlackoilInterface& linsolver,
                                       const double* gravity,
                                       const bool has_disgas,
                                       const bool has_vapoil,
                                       std::shared_ptr<EclipseState> eclipse_state,
                                       BlackoilOutputWriterEbos& output_writer,
                                       const std::vector<double>& threshold_pressures_by_face)
        : ebosSimulator_(ebosSimulator),
          param_(param),
          model_param_(param),
          solver_param_(param),
          props_(props),
          rock_comp_props_(rock_comp_props),
          gravity_(gravity),
          geo_(geo),
          solver_(linsolver),
          has_disgas_(has_disgas),
          has_vapoil_(has_vapoil),
          terminal_output_(param.getDefault("output_terminal", true)),
          output_writer_(output_writer),
          threshold_pressures_by_face_(threshold_pressures_by_face),
          is_parallel_run_( false )
    {

        // Misc init.
        const int num_cells = AutoDiffGrid::numCells(grid());
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }

        rateConverter_.reset(new RateConverterType(props_, std::vector<int>(AutoDiffGrid::numCells(grid()), 0)));

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
    SimulatorReport run(SimulatorTimer& timer,
                        ReservoirState& state)
    {
        WellState prev_well_state;

        if (output_writer_.isRestart()) {
            // This is a restart, populate WellState and ReservoirState state objects from restart file
            output_writer_.initFromRestartFile(props_.phaseUsage(), props_.permeability(), grid(), state, prev_well_state);
            initHydroCarbonState(state, props_.phaseUsage(), Opm::UgGridHelpers::numCells(grid()), has_disgas_, has_vapoil_);
        }

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

        std::string restorefilename = param_.getDefault("restorefile", std::string("") );
        if( ! restorefilename.empty() )
        {
            // -1 means that we'll take the last report step that was written
            //const int desiredRestoreStep = param_.getDefault("restorestep", int(-1) );

            //            output_writer_.restore( timer,
            //                                    state,
            //                                    prev_well_state,
            //                                    restorefilename,
            //                                    desiredRestoreStep );
        }

        unsigned int totalLinearizations = 0;
        unsigned int totalNonlinearIterations = 0;
        unsigned int totalLinearIterations = 0;
        bool is_well_potentials_computed = param_.getDefault("compute_well_potentials", false );
        std::vector<double> well_potentials;
        DynamicListEconLimited dynamic_list_econ_limited;

        bool ooip_computed = false;
        std::vector<int> fipnum_global = eclState().get3DProperties().getIntGridProperty("FIPNUM").getData();
        //Get compressed cell fipnum.
        std::vector<int> fipnum(Opm::UgGridHelpers::numCells(grid()));
        if (fipnum_global.empty()) {
            std::fill(fipnum.begin(), fipnum.end(), 0);
        } else {
            for (size_t c = 0; c < fipnum.size(); ++c) {
                fipnum[c] = fipnum_global[Opm::UgGridHelpers::globalCell(grid())[c]];
            }
        }
        std::vector<std::vector<double>> OOIP;

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
            WellsManager wells_manager(eclState(),
                                       timer.currentStepNum(),
                                       Opm::UgGridHelpers::numCells(grid()),
                                       Opm::UgGridHelpers::globalCell(grid()),
                                       Opm::UgGridHelpers::cartDims(grid()),
                                       Opm::UgGridHelpers::dimensions(grid()),
                                       Opm::UgGridHelpers::cell2Faces(grid()),
                                       Opm::UgGridHelpers::beginFaceCentroids(grid()),
                                       props_.permeability(),
                                       dynamic_list_econ_limited,
                                       is_parallel_run_,
                                       well_potentials );

            const Wells* wells = wells_manager.c_wells();
            WellState well_state;
            well_state.init(wells, state, prev_well_state);

            // give the polymer and surfactant simulators the chance to do their stuff
            handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // Compute reservoir volumes for RESV controls.
            computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            const std::vector<double> pv(geo_.poreVolume().data(), geo_.poreVolume().data() + geo_.poreVolume().size());
            const WellModel well_model(wells, model_param_, terminal_output_, pv);

            auto solver = createSolver(well_model);

            // Compute orignal FIP;
            if (!timer.initialStep() && !ooip_computed) {
                OOIP = solver->computeFluidInPlace(fipnum);
                FIPUnitConvert(eclState().getUnits(), OOIP);
                ooip_computed = true;
            }

            // write the inital state at the report stage
            if (timer.initialStep()) {
                //output_writer_.writeTimeStep( timer, state, well_state, solver->model() );
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

            solver->model().beginReportStep();

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

            solver->model().endReportStep();

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of nonlinear and linear Iterations
            totalLinearizations += solver->linearizations();
            totalNonlinearIterations += solver->nonlinearIterations();
            totalLinearIterations += solver->linearIterations();

            // Report timing.
            const double st = solver_timer.secsSinceStart();

            // Compute current FIP.
            std::vector<std::vector<double>> COIP;
            COIP = solver->computeFluidInPlace(fipnum);
            FIPUnitConvert(eclState().getUnits(), COIP);
            std::vector<double> OOIP_totals = FIPTotals(OOIP, state);
            std::vector<double> COIP_totals = FIPTotals(COIP, state);

            if ( terminal_output_ )
            {
                outputFluidInPlace(OOIP_totals, COIP_totals,eclState().getUnits(), 0);
                for (size_t reg = 0; reg < OOIP.size(); ++reg) {
                    outputFluidInPlace(OOIP[reg], COIP[reg], eclState().getUnits(), reg+1);
                }
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
            output_writer_.writeTimeStep( timer, state, well_state, solver->model() );

            prev_well_state = well_state;
            // The well potentials are only computed if they are needed
            // For now thay are only used to determine default guide rates for group controlled wells
            if ( is_well_potentials_computed ) {
                computeWellPotentials(wells, well_state, well_potentials);
            }

            updateListEconLimited(solver, eclState().getSchedule(), timer.currentStepNum(), wells,
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

    const Grid& grid() const
    { return ebosSimulator_.gridManager().grid(); }

protected:
    void handleAdditionalWellInflow(SimulatorTimer& timer,
                                    WellsManager& wells_manager,
                                    WellState& well_state,
                                    const Wells* wells)
    { }

    std::unique_ptr<Solver> createSolver(const WellModel& well_model)
    {
        auto model = std::unique_ptr<Model>(new Model(ebosSimulator_,
                                                      model_param_,
                                                      props_,
                                                      geo_,
                                                      rock_comp_props_,
                                                      well_model,
                                                      solver_,
                                                      terminal_output_));

        return std::unique_ptr<Solver>(new Solver(solver_param_, std::move(model)));
    }

    void computeRESV(const std::size_t step,
                     const Wells* wells,
                     const BlackoilState& x,
                     WellState& xw)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const auto w_ecl = eclState().getSchedule().getWells(step);
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
                rateConverter_->defineState(x, boost::any_cast<const ParallelISTLInformation&>(solver_.parallelInformation()));
            }
        }
        else
#endif
        {
            if ( global_number_resv_wells )
            {
                rateConverter_->defineState(x);
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
                        rateConverter_->calcCoeff(prates, fipreg, distr);

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
                            rateConverter_->calcCoeff(hrates, fipreg, distr);

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


    void computeWellPotentials(const Wells*                    wells,
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


    void updateListEconLimited(const std::unique_ptr<Solver>& solver,
                               const Schedule& schedule,
                               const int current_step,
                               const Wells* wells,
                               const WellState& well_state,
                               DynamicListEconLimited& list_econ_limited) const
    {
        solver->model().wellModel().updateListEconLimited(schedule, current_step, wells,
                                                          well_state, list_econ_limited);
    }

    void FIPUnitConvert(const UnitSystem& units,
                        std::vector<std::vector<double>>& fip)
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


    std::vector<double> FIPTotals(const std::vector<std::vector<double>>& fip, const ReservoirState& state)
    {
        std::vector<double> totals(6,0.0);
        for (int i = 0; i < 5; ++i) {
            for (size_t reg = 0; reg < fip.size(); ++reg) {
                totals[i] += fip[reg][i];
            }
        }
        const int numCells = Opm::AutoDiffGrid::numCells(grid());
        const auto& pv = geo_.poreVolume();
        double pv_hydrocarbon_sum = 0.0;
        double p_pv_hydrocarbon_sum = 0.0;

        for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            const auto& intQuants = *ebosSimulator_.model().cachedIntensiveQuantities(cellIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const double& p = fs.pressure(FluidSystem::oilPhaseIdx).value();
            const double hydrocarbon = fs.saturation(FluidSystem::oilPhaseIdx).value() + fs.saturation(FluidSystem::gasPhaseIdx).value();
            if ( ! is_parallel_run_ )
            {
                totals[5] += pv[cellIdx];
                pv_hydrocarbon_sum += pv[cellIdx] * hydrocarbon;
                p_pv_hydrocarbon_sum += p * pv[cellIdx] * hydrocarbon;
            }
            else {
                OPM_THROW(std::logic_error, "FIP not yet implemented for MPI");
            }
        }
        totals[6] = unit::convert::to( (p_pv_hydrocarbon_sum / pv_hydrocarbon_sum), unit::barsa);
        return totals;
    }



    void outputFluidInPlace(const std::vector<double>& oip, const std::vector<double>& cip, const UnitSystem& units, const int reg)
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


    const EclipseState& eclState() const
    { return *ebosSimulator_.gridManager().eclState(); }

    EclipseState& eclState()
    { return *ebosSimulator_.gridManager().eclState(); }

    // Data.
    Simulator& ebosSimulator_;

    typedef RateConverter::
    SurfaceToReservoirVoidage< BlackoilPropsAdInterface,
                               std::vector<int> > RateConverterType;
    typedef typename Solver::SolverParameters SolverParameters;

    const parameter::ParameterGroup param_;
    ModelParameters model_param_;
    SolverParameters solver_param_;

    // Observed objects.
    BlackoilPropsAdInterface& props_;
    const RockCompressibility* rock_comp_props_;
    const double* gravity_;
    // Solvers
    DerivedGeology& geo_;
    NewtonIterationBlackoilInterface& solver_;
    // Misc. data
    std::vector<int> allcells_;
    const bool has_disgas_;
    const bool has_vapoil_;
    bool       terminal_output_;
    // output_writer
    OutputWriter& output_writer_;
    std::unique_ptr<RateConverterType> rateConverter_;
    // Threshold pressures.
    std::vector<double> threshold_pressures_by_face_;
    // Whether this a parallel simulation or not
    bool is_parallel_run_;

};

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
