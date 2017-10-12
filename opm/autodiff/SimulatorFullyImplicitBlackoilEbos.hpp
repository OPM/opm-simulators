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

#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/autodiff/IterationReport.hpp>
#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilWellModel.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/SimFIBODetails.hpp>
#include <opm/autodiff/moduleVersion.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/core/utility/initHydroCarbonState.hpp>
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
    typedef BlackoilOutputWriter OutputWriter;
    typedef BlackoilModelEbos<TypeTag> Model;
    typedef BlackoilModelParameters ModelParameters;
    typedef NonlinearSolver<Model> Solver;
    typedef BlackoilWellModel<TypeTag> WellModel;
    typedef RateConverter::SurfaceToReservoirVoidage<FluidSystem, std::vector<int> > RateConverterType;


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
                                       const EclipseState& /* eclState */,
                                       OutputWriter& output_writer,
                                       const std::unordered_set<std::string>& defunct_well_names)
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
          rateConverter_(createRateConverter_()),
          defunct_well_names_( defunct_well_names ),
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
    SimulatorReport run(SimulatorTimer& timer,
                        ReservoirState& state)
    {
        WellState prev_well_state;

        ExtraData extra;

        failureReport_ = SimulatorReport();
        extractLegacyDepth_();

        // communicate the initial solution to ebos
        if (timer.initialStep()) {
            convertInput(/*iterationIdx=*/0, state, ebosSimulator_ );
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
        }

        if (output_writer_.isRestart()) {
            // This is a restart, populate WellState and ReservoirState state objects from restart file
            output_writer_.initFromRestartFile(phaseUsage_, grid(), state, prev_well_state, extra);
            initHydroCarbonState(state, phaseUsage_, Opm::UgGridHelpers::numCells(grid()), has_disgas_, has_vapoil_);
            initHysteresisParams(state);
            // communicate the restart solution to ebos
            convertInput(/*iterationIdx=*/0, state, ebosSimulator_ );
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
        }

        // Sync the overlap region of the inital solution. It was generated
        // from the ReservoirState which has wrong values in the ghost region
        // for some models (SPE9, Norne, Model 2)
        ebosSimulator_.model().syncOverlap();

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_writer_.outputDirectory() + "/step_timing.txt";
        std::ofstream tstep_os;

        if ( output_writer_.output() && output_writer_.isIORank() )
        {
            tstep_os.open(tstep_filename.c_str());
        }

        const auto& schedule = eclState().getSchedule();

        // adaptive time stepping
        const auto& events = schedule.getEvents();
        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        if( param_.getDefault("timestep.adaptive", true ) )
        {

            if (param_.getDefault("use_TUNING", false)) {
                adaptiveTimeStepping.reset( new AdaptiveTimeStepping( schedule.getTuning(), timer.currentStepNum(), param_, terminal_output_ ) );
            } else {
                adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_, terminal_output_ ) );
            }

            if (output_writer_.isRestart()) {
                if (extra.suggested_step > 0.0) {
                    adaptiveTimeStepping->setSuggestedNextStep(extra.suggested_step);
                }
            }
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
            initHydroCarbonState(state, phaseUsage_, Opm::UgGridHelpers::numCells(grid()), has_disgas_, has_vapoil_);
            initHysteresisParams(state);
            // communicate the restart solution to ebos
            convertInput(0, state, ebosSimulator_);
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
        }

        DynamicListEconLimited dynamic_list_econ_limited;
        SimulatorReport report;
        SimulatorReport stepReport;

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
        std::vector<std::vector<double>> originalFluidInPlace;
        std::vector<double> originalFluidInPlaceTotals;

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
                                       dynamic_list_econ_limited,
                                       is_parallel_run_,
                                       defunct_well_names_ );
            const Wells* wells = wells_manager.c_wells();
            WellState well_state;

            // The well state initialize bhp with the cell pressure in the top cell.
            // We must therefore provide it with updated cell pressures
            size_t nc = Opm::UgGridHelpers::numCells(grid());
            std::vector<double> cellPressures(nc, 0.0);
            const auto& gridView = ebosSimulator_.gridManager().gridView();
            ElementContext elemCtx(ebosSimulator_);
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (auto elemIt = gridView.template begin</*codim=*/0>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity) {
                    continue;
                }

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const double p = fs.pressure(FluidSystem::oilPhaseIdx).value();
                cellPressures[cellIdx] = p;
            }
            well_state.init(wells, cellPressures, prev_well_state, phaseUsage_);

            // give the polymer and surfactant simulators the chance to do their stuff
            handleAdditionalWellInflow(timer, wells_manager, well_state, wells);

            // Compute reservoir volumes for RESV controls.
            computeRESV(timer.currentStepNum(), wells, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            const auto& wells_ecl = eclState().getSchedule().getWells(timer.currentStepNum());
            extractLegacyCellPvtRegionIndex_();
            WellModel well_model(wells, &(wells_manager.wellCollection()), wells_ecl, model_param_,
                                 rateConverter_, terminal_output_, timer.currentStepNum(), legacyCellPvtRegionIdx_);

	    // handling MS well related
            if (model_param_.use_multisegment_well_) { // if we use MultisegmentWell model
                for (const auto& well : wells_ecl) {
		    // TODO: this is acutally not very accurate, because sometimes a deck just claims a MS well
		    // while keep the well shut. More accurately, we should check if the well exisits in the Wells
		    // structure here
                    if (well->isMultiSegment(timer.currentStepNum()) ) { // there is one well is MS well
                        well_state.initWellStateMSWell(wells, wells_ecl, timer.currentStepNum(), phaseUsage_, prev_well_state);
                        break;
                    }
                }
            }


            auto solver = createSolver(well_model);

            std::vector<std::vector<double>> currentFluidInPlace;
            std::vector<double> currentFluidInPlaceTotals;

            // Compute orignal fluid in place if this has not been done yet
            if (originalFluidInPlace.empty()) {
                originalFluidInPlace = solver->computeFluidInPlace(fipnum);
                originalFluidInPlaceTotals = FIPTotals(originalFluidInPlace);
                FIPUnitConvert(eclState().getUnits(), originalFluidInPlace);
                FIPUnitConvert(eclState().getUnits(), originalFluidInPlaceTotals);

                currentFluidInPlace = originalFluidInPlace;
                currentFluidInPlaceTotals = originalFluidInPlaceTotals;
            }

            // write the inital state at the report stage
            if (timer.initialStep()) {
                Dune::Timer perfTimer;
                perfTimer.start();

                if (terminal_output_) {
                    outputFluidInPlace(originalFluidInPlaceTotals, currentFluidInPlaceTotals,eclState().getUnits(), 0);
                    for (size_t reg = 0; reg < originalFluidInPlace.size(); ++reg) {
                        outputFluidInPlace(originalFluidInPlace[reg], currentFluidInPlace[reg], eclState().getUnits(), reg+1);
                    }
                }

                // No per cell data is written for initial step, but will be
                // for subsequent steps, when we have started simulating
                output_writer_.writeTimeStep( timer, state, well_state, solver->model() );

                report.output_write_time += perfTimer.stop();
            }

            if( terminal_output_ )
            {
                std::ostringstream step_msg;
                boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
                step_msg.imbue(std::locale(std::locale::classic(), facet));
                step_msg << "\nTime step " << std::setw(4) <<timer.currentStepNum()
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
                bool event = events.hasEvent(ScheduleEvents::NEW_WELL, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::INJECTION_UPDATE, timer.currentStepNum()) ||
                        events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE, timer.currentStepNum());
                stepReport = adaptiveTimeStepping->step( timer, *solver, state, well_state, event, output_writer_,
                                                         output_writer_.requireFIPNUM() ? &fipnum : nullptr );
                report += stepReport;
                failureReport_ += adaptiveTimeStepping->failureReport();
            }
            else {
                // solve for complete report step
                stepReport = solver->step(timer, state, well_state);
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

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // update timing.
            report.solver_time += solver_timer.secsSinceStart();

            if ( output_writer_.output() && output_writer_.isIORank() )
            {
                stepReport.reportParam(tstep_os);
            }

            // We don't need the reservoir state anymore. It is just passed around to avoid
            // code duplication. Pass empty state instead.
            if (timer.initialStep()) {
                ReservoirState stateTrivial(0,0,0);
                state = stateTrivial;
            }

            // Increment timer, remember well state.
            ++timer;

            // Compute current fluid in place.
            currentFluidInPlace = solver->computeFluidInPlace(fipnum);
            currentFluidInPlaceTotals = FIPTotals(currentFluidInPlace);

            const std::string version = moduleVersionName();

            FIPUnitConvert(eclState().getUnits(), currentFluidInPlace);
            FIPUnitConvert(eclState().getUnits(), currentFluidInPlaceTotals);

            if (terminal_output_ )
            {
                outputTimestampFIP(timer, version);
                outputFluidInPlace(originalFluidInPlaceTotals, currentFluidInPlaceTotals,eclState().getUnits(), 0);
                for (size_t reg = 0; reg < originalFluidInPlace.size(); ++reg) {
                    outputFluidInPlace(originalFluidInPlace[reg], currentFluidInPlace[reg], eclState().getUnits(), reg+1);
                }

                std::string msg;
                msg =
                    "Time step took " + std::to_string(solver_timer.secsSinceStart()) + " seconds; "
                    "total solver time " + std::to_string(report.solver_time) + " seconds.";
                OpmLog::note(msg);
            }

            // write simulation state at the report stage
            Dune::Timer perfTimer;
            perfTimer.start();
            const double nextstep = adaptiveTimeStepping ? adaptiveTimeStepping->suggestedNextStep() : -1.0;
            output_writer_.writeTimeStep( timer, state, well_state, solver->model(), false, nextstep, report);
            report.output_write_time += perfTimer.stop();

            prev_well_state = well_state;

            updateListEconLimited(solver, eclState().getSchedule(), timer.currentStepNum(), wells,
                                  well_state, dynamic_list_econ_limited);
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
    void handleAdditionalWellInflow(SimulatorTimer& /*timer*/,
                                    WellsManager& /* wells_manager */,
                                    WellState& /* well_state */,
                                    const Wells* /* wells */)
    {
    }

    std::unique_ptr<Solver> createSolver(WellModel& well_model)
    {
        const auto& gridView = ebosSimulator_.gridView();
        const PhaseUsage& phaseUsage = phaseUsage_;
        const std::vector<bool> activePhases = detail::activePhases(phaseUsage);
        const double gravity = ebosSimulator_.problem().gravity()[2];

        // calculate the number of elements of the compressed sequential grid. this needs
        // to be done in two steps because the dune communicator expects a reference as
        // argument for sum()
        int globalNumCells = gridView.size(/*codim=*/0);
        globalNumCells = gridView.comm().sum(globalNumCells);

        well_model.init(phaseUsage,
                        activePhases,
                        gravity,
                        legacyDepth_,
                        globalNumCells,
                        grid());
        auto model = std::unique_ptr<Model>(new Model(ebosSimulator_,
                                                      model_param_,
                                                      well_model,
                                                      rateConverter_,
                                                      solver_,
                                                      terminal_output_));

        return std::unique_ptr<Solver>(new Solver(solver_param_, std::move(model)));
    }

    void computeRESV(const std::size_t step,
                     const Wells* wells,
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
                rateConverter_.template defineState<ElementContext>(ebosSimulator_);
            }
        }
        else
#endif
        {
            if ( global_number_resv_wells )
            {
                rateConverter_.template defineState<ElementContext>(ebosSimulator_);
            }
        }

        if (! resv_wells.empty()) {
            const PhaseUsage&                    pu = phaseUsage_;
            const std::vector<double>::size_type np = phaseUsage_.num_phases;

            std::vector<double> distr (np);
            std::vector<double> hrates(np);
            std::vector<double> prates(np);

            for (std::vector<int>::const_iterator
                     rp = resv_wells.begin(), e = resv_wells.end();
                 rp != e; ++rp)
            {
                WellControls* ctrl = wells->ctrls[*rp];
                const bool is_producer = wells->type[*rp] == PRODUCER;
                const int well_cell_top = wells->well_cells[wells->well_connpos[*rp]];
                const auto& eclProblem = ebosSimulator_.problem();
                const int pvtreg = eclProblem.pvtRegionIndex(well_cell_top);

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
                        rateConverter_.calcCoeff(fipreg, pvtreg, distr);

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
                            rateConverter_.calcCoeff(fipreg, pvtreg, distr);

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
        for (size_t i = 0; i < fip.size(); ++i) {
            FIPUnitConvert(units, fip[i]);
        }
    }


    void FIPUnitConvert(const UnitSystem& units,
                        std::vector<double>& fip)
    {
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            fip[0] = unit::convert::to(fip[0], unit::stb);
            fip[1] = unit::convert::to(fip[1], unit::stb);
            fip[2] = unit::convert::to(fip[2], 1000*unit::cubic(unit::feet));
            fip[3] = unit::convert::to(fip[3], 1000*unit::cubic(unit::feet));
            fip[4] = unit::convert::to(fip[4], unit::stb);
            fip[5] = unit::convert::to(fip[5], unit::stb);
            fip[6] = unit::convert::to(fip[6], unit::psia);
        }
        else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            fip[6] = unit::convert::to(fip[6], unit::barsa);
        }
        else {
            OPM_THROW(std::runtime_error, "Unsupported unit type for fluid in place output.");
        }
    }


    std::vector<double> FIPTotals(const std::vector<std::vector<double>>& fip)
    {
        std::vector<double> totals(7,0.0);
        for (int i = 0; i < 5; ++i) {
            for (size_t reg = 0; reg < fip.size(); ++reg) {
                totals[i] += fip[reg][i];
            }
        }

        const auto& gridView = ebosSimulator_.gridManager().gridView();
        const auto& comm = gridView.comm();
        double pv_hydrocarbon_sum = 0.0;
        double p_pv_hydrocarbon_sum = 0.0;

        ElementContext elemCtx(ebosSimulator_);
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (auto elemIt = gridView.template begin</*codim=*/0>();
             elemIt != elemEndIt;
             ++elemIt)
        {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const double p = fs.pressure(FluidSystem::oilPhaseIdx).value();
            const double hydrocarbon = fs.saturation(FluidSystem::oilPhaseIdx).value() + fs.saturation(FluidSystem::gasPhaseIdx).value();

            // calculate the pore volume of the current cell. Note that the
            // porosity returned by the intensive quantities is defined as the
            // ratio of pore space to total cell volume and includes all pressure
            // dependent (-> rock compressibility) and static modifiers (MULTPV,
            // MULTREGP, NTG, PORV, MINPV and friends). Also note that because of
            // this, the porosity returned by the intensive quantities can be
            // outside of the physical range [0, 1] in pathetic cases.
            const double pv =
                ebosSimulator_.model().dofTotalVolume(cellIdx)
                * intQuants.porosity().value();

            totals[5] += pv;
            pv_hydrocarbon_sum += pv*hydrocarbon;
            p_pv_hydrocarbon_sum += p*pv*hydrocarbon;
        }

        pv_hydrocarbon_sum = comm.sum(pv_hydrocarbon_sum);
        p_pv_hydrocarbon_sum = comm.sum(p_pv_hydrocarbon_sum);
        totals[5] = comm.sum(totals[5]);
        totals[6] = (p_pv_hydrocarbon_sum / pv_hydrocarbon_sum);

        return totals;
    }


    void outputTimestampFIP(SimulatorTimer& timer, const std::string version)
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
                ss << "                                                  : Pressure is weighted by hydrocarbon pore volume :\n"
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
    { return ebosSimulator_.gridManager().eclState(); }

    void extractLegacyCellPvtRegionIndex_()
    {
        const auto& grid = ebosSimulator_.gridManager().grid();
        const auto& eclProblem = ebosSimulator_.problem();
        const unsigned numCells = grid.size(/*codim=*/0);

        legacyCellPvtRegionIdx_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            legacyCellPvtRegionIdx_[cellIdx] =
                eclProblem.pvtRegionIndex(cellIdx);
        }
    }

    void initHysteresisParams(ReservoirState& state) {
        const int num_cells = Opm::UgGridHelpers::numCells(grid());

        typedef std::vector<double> VectorType;

        const VectorType& somax = state.getCellData( "SOMAX" );

        for (int cellIdx = 0; cellIdx < num_cells; ++cellIdx) {
            ebosSimulator_.model().setMaxOilSaturation(somax[cellIdx], cellIdx);
        }

        if (ebosSimulator_.problem().materialLawManager()->enableHysteresis()) {
            auto matLawManager = ebosSimulator_.problem().materialLawManager();

            VectorType& pcSwMdc_ow = state.getCellData( "PCSWMDC_OW" );
            VectorType& krnSwMdc_ow = state.getCellData( "KRNSWMDC_OW" );

            VectorType& pcSwMdc_go = state.getCellData( "PCSWMDC_GO" );
            VectorType& krnSwMdc_go = state.getCellData( "KRNSWMDC_GO" );

            for (int cellIdx = 0; cellIdx < num_cells; ++cellIdx) {
                matLawManager->setOilWaterHysteresisParams(
                        pcSwMdc_ow[cellIdx],
                        krnSwMdc_ow[cellIdx],
                        cellIdx);
                matLawManager->setGasOilHysteresisParams(
                        pcSwMdc_go[cellIdx],
                        krnSwMdc_go[cellIdx],
                        cellIdx);
            }
        }
    }

    void extractLegacyDepth_()
    {
        const auto& grid = ebosSimulator_.gridManager().grid();
        const unsigned numCells = grid.size(/*codim=*/0);

        legacyDepth_.resize(numCells);
        for (unsigned cellIdx = 0; cellIdx < numCells; ++cellIdx) {
            legacyDepth_[cellIdx] =
                grid.cellCenterDepth(cellIdx);
        }
    }

    // Used to convert initial Reservoirstate to primary variables in the SolutionVector
    void convertInput( const int iterationIdx,
                       const ReservoirState& reservoirState,
                       Simulator& simulator ) const
    {
        SolutionVector& solution = simulator.model().solution( 0 /* timeIdx */ );
        const Opm::PhaseUsage pu = phaseUsage_;

        const std::vector<bool> active = detail::activePhases(pu);
        bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
        bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);

        const int numCells = reservoirState.numCells();
        const int numPhases = phaseUsage_.num_phases;
        const auto& oilPressure = reservoirState.pressure();
        const auto& saturations = reservoirState.saturation();
        const auto& rs          = reservoirState.gasoilratio();
        const auto& rv          = reservoirState.rv();
        for( int cellIdx = 0; cellIdx<numCells; ++cellIdx )
        {
            // set non-switching primary variables
            PrimaryVariables& cellPv = solution[ cellIdx ];
            // set water saturation
            if ( active[Water] ) {
                cellPv[BlackoilIndices::waterSaturationIdx] = saturations[cellIdx*numPhases + pu.phase_pos[Water]];
            }

            if (has_solvent) {
                cellPv[BlackoilIndices::solventSaturationIdx] = reservoirState.getCellData( reservoirState.SSOL )[cellIdx];
            }

            if (has_polymer) {
                cellPv[BlackoilIndices::polymerConcentrationIdx] = reservoirState.getCellData( reservoirState.POLYMER )[cellIdx];
            }


            // set switching variable and interpretation
            if ( active[Gas] ) {
                if( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::OilOnly && has_disgas_ )
                {
                    cellPv[BlackoilIndices::compositionSwitchIdx] = rs[cellIdx];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[cellIdx];
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_po_Rs );
                }
                else if( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::GasOnly && has_vapoil_ )
                {
                    // this case (-> gas only with vaporized oil in the gas) is
                    // relatively expensive as it requires to compute the capillary
                    // pressure in order to get the gas phase pressure. (the reason why
                    // ebos uses the gas pressure here is that it makes the common case
                    // of the primary variable switching code fast because to determine
                    // whether the oil phase appears one needs to compute the Rv value
                    // for the saturated gas phase and if this is not available as a
                    // primary variable, it needs to be computed.) luckily for here, the
                    // gas-only case is not too common, so the performance impact of this
                    // is limited.
                    typedef Opm::SimpleModularFluidState<double,
                            /*numPhases=*/3,
                            /*numComponents=*/3,
                            FluidSystem,
                            /*storePressure=*/false,
                            /*storeTemperature=*/false,
                            /*storeComposition=*/false,
                            /*storeFugacity=*/false,
                            /*storeSaturation=*/true,
                            /*storeDensity=*/false,
                            /*storeViscosity=*/false,
                            /*storeEnthalpy=*/false> SatOnlyFluidState;
                    SatOnlyFluidState fluidState;
                    if ( active[Water] ) {
                        fluidState.setSaturation(FluidSystem::waterPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Water]]);
                    }
                    else {
                        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
                    }
                    fluidState.setSaturation(FluidSystem::oilPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Oil]]);
                    fluidState.setSaturation(FluidSystem::gasPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Gas]]);

                    double pC[/*numPhases=*/3] = { 0.0, 0.0, 0.0 };
                    const MaterialLawParams& matParams = simulator.problem().materialLawParams(cellIdx);
                    MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                    double pg = oilPressure[cellIdx] + (pC[FluidSystem::gasPhaseIdx] - pC[FluidSystem::oilPhaseIdx]);

                    cellPv[BlackoilIndices::compositionSwitchIdx] = rv[cellIdx];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = pg;
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_pg_Rv );
                }
                else
                {
                    assert( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::GasAndOil);
                    cellPv[BlackoilIndices::compositionSwitchIdx] = saturations[cellIdx*numPhases + pu.phase_pos[Gas]];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[ cellIdx ];
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_po_Sg );
                }
            } else {
                // for oil-water case oil pressure should be used as primary variable
                cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[cellIdx];
            }
        }

        // store the solution at the beginning of the time step
        if( iterationIdx == 0 )
        {
            simulator.model().solution( 1 /* timeIdx */ ) = solution;
        }
    }

    RateConverterType createRateConverter_() {
        RateConverterType rate_converter(phaseUsage_,
                                         std::vector<int>(AutoDiffGrid::numCells(grid()), 0)); // FIP = 0
        return rate_converter;
    }


    // Data.
    Simulator& ebosSimulator_;

    std::vector<int> legacyCellPvtRegionIdx_;
    std::vector<double> legacyDepth_;
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
    RateConverterType rateConverter_;
    // The names of wells that should be defunct
    // (e.g. in a parallel run when they are handeled by
    // a different process)
    std::unordered_set<std::string> defunct_well_names_;

    // Whether this a parallel simulation or not
    bool is_parallel_run_;

};

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
