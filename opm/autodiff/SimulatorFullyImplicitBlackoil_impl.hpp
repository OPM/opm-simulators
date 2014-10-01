/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
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

#include  <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/FullyImplicitBlackoilSolver.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>

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
    class SimulatorFullyImplicitBlackoil<T>::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const Grid& grid,
             const DerivedGeology& geo,
             BlackoilPropsAdInterface& props,
             const RockCompressibility* rock_comp_props,
             NewtonIterationBlackoilInterface& linsolver,
             const double* gravity,
             bool has_disgas,
             bool has_vapoil,
             std::shared_ptr<EclipseState> eclipse_state,
             EclipseWriter& output_writer,
             const std::vector<double>& threshold_pressures_by_face);

        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state);

    private:
        // Data.
        typedef RateConverter::
        SurfaceToReservoirVoidage< BlackoilPropsAdInterface,
                                   std::vector<int> > RateConverterType;

        const parameter::ParameterGroup param_;

        // Parameters for output.
        bool output_;
        bool output_vtk_;
        std::string output_dir_;
        int output_interval_;
        // Observed objects.
        const Grid& grid_;
        BlackoilPropsAdInterface& props_;
        const RockCompressibility* rock_comp_props_;
        const double* gravity_;
        // Solvers
        const DerivedGeology& geo_;
        NewtonIterationBlackoilInterface& solver_;
        // Misc. data
        std::vector<int> allcells_;
        const bool has_disgas_;
        const bool has_vapoil_;
        // eclipse_state
        std::shared_ptr<EclipseState> eclipse_state_;
        // output_writer
        EclipseWriter& output_writer_;
        RateConverterType rateConverter_;
        // Threshold pressures.
        std::vector<double> threshold_pressures_by_face_;

        void
        computeRESV(const std::size_t               step,
                    const Wells*                    wells,
                    const BlackoilState&            x,
                    WellStateFullyImplicitBlackoil& xw);
    };




    template<class T>
    SimulatorFullyImplicitBlackoil<T>::SimulatorFullyImplicitBlackoil(const parameter::ParameterGroup& param,
                                                                   const Grid& grid,
                                                                   const DerivedGeology& geo,
                                                                   BlackoilPropsAdInterface& props,
                                                                   const RockCompressibility* rock_comp_props,
                                                                   NewtonIterationBlackoilInterface& linsolver,
                                                                   const double* gravity,
                                                                   const bool has_disgas,
                                                                   const bool has_vapoil,
                                                                   std::shared_ptr<EclipseState> eclipse_state,
                                                                   EclipseWriter& output_writer,
                                                                   const std::vector<double>& threshold_pressures_by_face)

    {
        pimpl_.reset(new Impl(param, grid, geo, props, rock_comp_props, linsolver, gravity, has_disgas, has_vapoil,
                              eclipse_state, output_writer, threshold_pressures_by_face));
    }





    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoil<T>::run(SimulatorTimer& timer,
                                                        BlackoilState& state)
    {
        return pimpl_->run(timer, state);
    }



    static void outputWellStateMatlab(const Opm::WellStateFullyImplicitBlackoil& well_state,
                                  const int step,
                                  const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["bhp"] = &well_state.bhp();
        dm["wellrates"] = &well_state.wellRates();

        // Write data (not grid) in Matlab format
        for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
            std::ostringstream fname;
            fname << output_dir << "/" << it->first;
            boost::filesystem::path fpath = fname.str();
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error,"Creating directories failed: " << fpath);
            }
            fname << "/" << std::setw(3) << std::setfill('0') << step << ".txt";
            std::ofstream file(fname.str().c_str());
            if (!file) {
                OPM_THROW(std::runtime_error,"Failed to open " << fname.str());
            }
            file.precision(15);
            const std::vector<double>& d = *(it->second);
            std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
        }
    }

#if 0
    static void outputWaterCut(const Opm::Watercut& watercut,
                               const std::string& output_dir)
    {
        // Write water cut curve.
        std::string fname = output_dir  + "/watercut.txt";
        std::ofstream os(fname.c_str());
        if (!os) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname);
        }
        watercut.write(os);
    }

    static void outputWellReport(const Opm::WellReport& wellreport,
                                 const std::string& output_dir)
    {
        // Write well report.
        std::string fname = output_dir  + "/wellreport.txt";
        std::ofstream os(fname.c_str());
        if (!os) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname);
        }
        wellreport.write(os);
    }
#endif


    // \TODO: Treat bcs.
    template<class T>
    SimulatorFullyImplicitBlackoil<T>::Impl::Impl(const parameter::ParameterGroup& param,
                                               const Grid& grid,
                                               const DerivedGeology& geo,
                                               BlackoilPropsAdInterface& props,
                                               const RockCompressibility* rock_comp_props,
                                               NewtonIterationBlackoilInterface& linsolver,
                                               const double* gravity,
                                               const bool has_disgas,
                                               const bool has_vapoil,
                                               std::shared_ptr<EclipseState> eclipse_state,
                                               EclipseWriter& output_writer,
                                               const std::vector<double>& threshold_pressures_by_face)
        : param_(param),
          grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          gravity_(gravity),
          geo_(geo),
          solver_(linsolver),
          has_disgas_(has_disgas),
          has_vapoil_(has_vapoil),
          eclipse_state_(eclipse_state),
          output_writer_(output_writer),
          rateConverter_(props_, std::vector<int>(AutoDiffGrid::numCells(grid_), 0)),
          threshold_pressures_by_face_(threshold_pressures_by_face)
    {
        // For output.
        output_ = param.getDefault("output", true);
        if (output_) {
            output_vtk_ = param.getDefault("output_vtk", true);
            output_dir_ = param.getDefault("output_dir", std::string("output"));
            // Ensure that output dir exists
            boost::filesystem::path fpath(output_dir_);
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }
            output_interval_ = param.getDefault("output_interval", 1);
        }

        // Misc init.
        const int num_cells = AutoDiffGrid::numCells(grid);
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
    }




    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoil<T>::Impl::run(SimulatorTimer& timer,
                                                                 BlackoilState& state)
    {
        WellStateFullyImplicitBlackoil prev_well_state;

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_dir_ + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        double lastSubStep = timer.currentStepLength();

        typename FullyImplicitBlackoilSolver<T>::SolverParameter solverParam( param_ );

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep.
            step_timer.start();
            timer.report(std::cout);

            // Create wells and well state.
            WellsManager wells_manager(eclipse_state_,
                                       timer.currentStepNum(),
                                       Opm::UgGridHelpers::numCells(grid_),
                                       Opm::UgGridHelpers::globalCell(grid_),
                                       Opm::UgGridHelpers::cartDims(grid_),
                                       Opm::UgGridHelpers::dimensions(grid_),
                                       Opm::UgGridHelpers::beginCellCentroids(grid_),
                                       Opm::UgGridHelpers::cell2Faces(grid_),
                                       Opm::UgGridHelpers::beginFaceCentroids(grid_),
                                       props_.permeability());
            const Wells* wells = wells_manager.c_wells();
            WellStateFullyImplicitBlackoil well_state;
            well_state.init(wells, state);
            if (timer.currentStepNum() != 0) {
                // Transfer previous well state to current.
                well_state.partialCopy(prev_well_state, *wells, prev_well_state.numWells());
            }

            // Output state at start of time step.
            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
                outputWellStateMatlab(well_state,timer.currentStepNum(), output_dir_);
            }
            if (output_) {
                if (timer.currentStepNum() == 0) {
                    output_writer_.writeInit(timer);
                }
                output_writer_.writeTimeStep(timer, state, well_state.basicWellState());
            }

            // Max oil saturation (for VPPARS), hysteresis update.
            props_.updateSatOilMax(state.saturation());
            props_.updateSatHyst(state.saturation(), allcells_);

            // Compute reservoir volumes for RESV controls.
            computeRESV(timer.currentStepNum(), wells, state, well_state);

            // Run a multiple steps of the solver depending on the time step control.
            solver_timer.start();

            FullyImplicitBlackoilSolver<T> solver(solverParam, grid_, props_, geo_, rock_comp_props_, *wells, solver_, has_disgas_, has_vapoil_);
            if (!threshold_pressures_by_face_.empty()) {
                solver.setThresholdPressures(threshold_pressures_by_face_);
            }

            const bool subStepping = false;
            if( subStepping )
            {

                // create sub step simulator timer with previously used sub step size
                SubStepSimulatorTimer subStepper( timer, lastSubStep );

                while( ! subStepper.done() )
                {
                    const double dt_new = solver.step(subStepper.currentStepLength(), state, well_state);
                    subStepper.next( dt_new );
                }

                subStepper.report( std::cout );

                // store last small time step for next reportStep
                lastSubStep = subStepper.currentStepLength();

                std::cout << "Last suggested step size = " << lastSubStep << std::endl;
                if( lastSubStep != lastSubStep )
                    lastSubStep = timer.currentStepLength();
            }
            else
                solver.step(timer.currentStepLength(), state, well_state);

            // take time that was used to solve system for this reportStep
            solver_timer.stop();

            // Report timing.
            const double st = solver_timer.secsSinceStart();
            std::cout << "Fully implicit solver took: " << st << " seconds." << std::endl;
            stime += st;
            if (output_) {
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
        if (output_) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
            outputWellStateMatlab(prev_well_state, timer.currentStepNum(), output_dir_);
            output_writer_.writeTimeStep(timer, state, prev_well_state.basicWellState());
        }

        // Stop timer and create timing report
        total_timer.stop();
        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
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
        is_resv_prod(const Wells& wells,
                     const int    w)
        {
            return ((wells.type[w] == PRODUCER) &&
                    (0 <= resv_control(wells.ctrls[w])));
        }

        inline bool
        is_resv_prod(const WellMap&     wmap,
                     const std::string& name,
                     const std::size_t  step)
        {
            bool match = false;

            WellMap::const_iterator i = wmap.find(name);

            if (i != wmap.end()) {
                WellConstPtr wp = i->second;

                match = (wp->isProducer(step) &&
                         wp->getProductionProperties(step)
                         .hasProductionControl(WellProducer::RESV));
            }

            return match;
        }

        inline std::vector<int>
        resvProducers(const Wells&      wells,
                      const std::size_t step,
                      const WellMap&    wmap)
        {
            std::vector<int> resv_prod;

            for (int w = 0, nw = wells.number_of_wells; w < nw; ++w) {
                if (is_resv_prod(wells, w) ||
                    ((wells.name[w] != 0) &&
                     is_resv_prod(wmap, wells.name[w], step)))
                {
                    resv_prod.push_back(w);
                }
            }

            return resv_prod;
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
    SimulatorFullyImplicitBlackoil<T>::
    Impl::computeRESV(const std::size_t               step,
                      const Wells*                    wells,
                      const BlackoilState&            x,
                      WellStateFullyImplicitBlackoil& xw)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const std::vector<WellConstPtr>& w_ecl = eclipse_state_->getSchedule()->getWells(step);
        const WellMap& wmap = SimFIBODetails::mapWells(w_ecl);

        const std::vector<int>& resv_prod =
            SimFIBODetails::resvProducers(*wells, step, wmap);

        if (! resv_prod.empty()) {
            const PhaseUsage&                    pu = props_.phaseUsage();
            const std::vector<double>::size_type np = props_.numPhases();

            rateConverter_.defineState(x);

            std::vector<double> distr (np);
            std::vector<double> hrates(np);
            std::vector<double> prates(np);

            for (std::vector<int>::const_iterator
                     rp = resv_prod.begin(), e = resv_prod.end();
                 rp != e; ++rp)
            {
                WellControls* ctrl = wells->ctrls[*rp];

                // RESV control mode, all wells
                {
                    const int rctrl = SimFIBODetails::resv_control(ctrl);

                    if (0 <= rctrl) {
                        const std::vector<double>::size_type off = (*rp) * np;

                        // Convert to positive rates to avoid issues
                        // in coefficient calculations.
                        std::transform(xw.wellRates().begin() + (off + 0*np),
                                       xw.wellRates().begin() + (off + 1*np),
                                       prates.begin(), std::negate<double>());

                        const int fipreg = 0; // Hack.  Ignore FIP regions.
                        rateConverter_.calcCoeff(prates, fipreg, distr);

                        well_controls_iset_distr(ctrl, rctrl, & distr[0]);
                    }
                }

                // RESV control, WCONHIST wells.  A bit of duplicate
                // work, regrettably.
                if (wells->name[*rp] != 0) {
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

                            const int ok =
                                well_controls_add_new(RESERVOIR_RATE, target,
                                                      & distr[0], ctrl);

                            if (ok != 0) {
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
