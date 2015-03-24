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

//#include <opm/polymer/fullyimplicit/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/polymer/fullyimplicit/SimulatorFullyImplicitBlackoilPolymer.hpp>
#include <opm/polymer/fullyimplicit/FullyImplicitBlackoilPolymerSolver.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/PolymerInflow.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
//#include <opm/core/simulator/AdaptiveSimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

//#include <opm/core/simulator/AdaptiveTimeStepping.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <cstddef>
#include <cassert>
#include <functional>
#include <memory>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{
    template<class T>
    class SimulatorFullyImplicitBlackoilPolymer<T>::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const Grid& grid,
             const DerivedGeology& geo,
             BlackoilPropsAdInterface& props,
             const PolymerPropsAd& polymer_props,
             const RockCompressibility* rock_comp_props,
             NewtonIterationBlackoilInterface& linsolver,
             const double* gravity,
             bool has_disgas,
             bool has_vapoil,
             bool has_polymer,
             std::shared_ptr<EclipseState> eclipse_state,
             BlackoilOutputWriter& output_writer,
             Opm::DeckConstPtr& deck,
             const std::vector<double>& threshold_pressures_by_face);

        SimulatorReport run(SimulatorTimer& timer,
                            PolymerBlackoilState& state);

    private:
        // Data.
        typedef RateConverter::
        SurfaceToReservoirVoidage< BlackoilPropsAdInterface,
                                   std::vector<int> > RateConverterType;

        const parameter::ParameterGroup param_;

        // Observed objects.
        const Grid& grid_;
        BlackoilPropsAdInterface& props_;
        const PolymerPropsAd& polymer_props_;
        const RockCompressibility* rock_comp_props_;
        const double* gravity_;
        // Solvers
        const DerivedGeology& geo_;
        NewtonIterationBlackoilInterface& solver_;
        // Misc. data
        std::vector<int> allcells_;
        const bool has_disgas_;
        const bool has_vapoil_;
        const bool has_polymer_;
        bool       terminal_output_;
        // eclipse_state
        std::shared_ptr<EclipseState> eclipse_state_;
        // output_writer
        BlackoilOutputWriter& output_writer_;
        Opm::DeckConstPtr& deck_;
        RateConverterType rateConverter_;
        // Threshold pressures.
        std::vector<double> threshold_pressures_by_face_;

        void
        computeRESV(const std::size_t               step,
                    const Wells*                    wells,
                    const BlackoilState&     x,
                    WellStateFullyImplicitBlackoil& xw);
    };




    template<class T>
    SimulatorFullyImplicitBlackoilPolymer<T>::SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
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

    {
        pimpl_.reset(new Impl(param, grid, geo, props, polymer_props, rock_comp_props, linsolver, gravity, has_disgas, has_vapoil, has_polymer,
                              eclipse_state, output_writer, deck, threshold_pressures_by_face));
    }





    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoilPolymer<T>::run(SimulatorTimer& timer,
                                                                  PolymerBlackoilState& state)
    {
        return pimpl_->run(timer, state);
    }


    // \TODO: Treat bcs.
    template<class T>
    SimulatorFullyImplicitBlackoilPolymer<T>::Impl::Impl(const parameter::ParameterGroup& param,
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
        : param_(param),
          grid_(grid),
          props_(props),
          polymer_props_(polymer_props),
          rock_comp_props_(rock_comp_props),
          gravity_(gravity),
          geo_(geo),
          solver_(linsolver),
          has_disgas_(has_disgas),
          has_vapoil_(has_vapoil),
          has_polymer_(has_polymer),
          terminal_output_(param.getDefault("output_terminal", true)),
          eclipse_state_(eclipse_state),
          output_writer_(output_writer),
          deck_(deck),
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
                terminal_output_= (info.communicator().rank()==0);
            }
        }
#endif
    }




    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoilPolymer<T>::Impl::run(SimulatorTimer& timer,
                                                                        PolymerBlackoilState& state)
    {
        WellStateFullyImplicitBlackoil prev_well_state;

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_writer_.outputDirectory() + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        typename FullyImplicitBlackoilPolymerSolver<T>::SolverParameter solverParam( param_ );

        //adaptive time stepping
        //        std::unique_ptr< AdaptiveTimeStepping > adaptiveTimeStepping;
        //        if( param_.getDefault("timestep.adaptive", bool(false) ) )
        //        {
        //            adaptiveTimeStepping.reset( new AdaptiveTimeStepping( param_ ) );
        //        }

        // init output writer
        output_writer_.writeInit( timer );

        std::string restorefilename = param_.getDefault("restorefile", std::string("") );
        if( ! restorefilename.empty() )
        {
            // -1 means that we'll take the last report step that was written
            const int desiredRestoreStep = param_.getDefault("restorestep", int(-1) );
            output_writer_.restore( timer, state.blackoilState(), prev_well_state, restorefilename, desiredRestoreStep );
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
                                       props_.permeability());
            const Wells* wells = wells_manager.c_wells();
            WellStateFullyImplicitBlackoil well_state;
            well_state.init(wells, state.blackoilState(), prev_well_state);

            // compute polymer inflow
            std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
            if (deck_->hasKeyword("WPOLYMER")) {
                if (wells_manager.c_wells() == 0) {
                    OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
                }
                polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, eclipse_state_, *wells, Opm::UgGridHelpers::numCells(grid_), timer.currentStepNum()));
            } else {
                polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                                1.0*Opm::unit::day,
                                                                0.0));
            }
            std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(grid_));
            polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(), 
                                                timer.simulationTimeElapsed() + timer.currentStepLength(),
                                                polymer_inflow_c);
            
            // write simulation state at the report stage
            output_writer_.writeTimeStep( timer, state.blackoilState(), well_state );

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            computeRESV(timer.currentStepNum(), wells, state.blackoilState(), well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            FullyImplicitBlackoilPolymerSolver<T> solver(solverParam, grid_, props_, geo_, rock_comp_props_, polymer_props_, wells, solver_, has_disgas_, has_vapoil_, has_polymer_, terminal_output_);
            if (!threshold_pressures_by_face_.empty()) {
                solver.setThresholdPressures(threshold_pressures_by_face_);
            }

            // If sub stepping is enabled allow the solver to sub cycle
            // in case the report steps are to large for the solver to converge
            //
            // \Note: The report steps are met in any case
            // \Note: The sub stepping will require a copy of the state variables
            //            if( adaptiveTimeStepping ) {
            //                adaptiveTimeStepping->step( timer, solver, state, well_state,  output_writer_ );
            //            } else {
                // solve for complete report step
            solver.step(timer.currentStepLength(), state, well_state, polymer_inflow_c);
                //            }

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // accumulate the number of Newton and Linear Iterations
            totalNewtonIterations += solver.newtonIterations();
            totalLinearIterations += solver.linearIterations();

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
        output_writer_.writeTimeStep( timer, state.blackoilState(), prev_well_state );

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

    template <class T>
    void
    SimulatorFullyImplicitBlackoilPolymer<T>::
    Impl::computeRESV(const std::size_t               step,
                      const Wells*                    wells,
                      const BlackoilState&            x,
                      WellStateFullyImplicitBlackoil& xw)
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
