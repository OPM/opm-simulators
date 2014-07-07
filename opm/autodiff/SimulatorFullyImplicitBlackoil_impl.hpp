/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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
             ScheduleConstPtr schedule,
             const Grid& grid,
             const DerivedGeology& geo,
             BlackoilPropsAdInterface& props,
             const RockCompressibility* rock_comp_props,
             const Wells* wells,
             NewtonIterationBlackoilInterface& linsolver,
             const double* gravity,
             bool has_disgas,
             bool has_vapoil );

        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state,
                            WellStateFullyImplicitBlackoil& well_state);

    private:
        // Data.

        typedef RateConverter::
        SurfaceToReservoirVoidage< BlackoilPropsAdInterface,
                                   std::vector<int> > RateConverterType;

        // Parameters for output.
        bool output_;
        bool output_vtk_;
        std::string output_dir_;
        int output_interval_;
        // Observed objects.
        ScheduleConstPtr schedule_;
        const Grid& grid_;
        BlackoilPropsAdInterface& props_;
        const RockCompressibility* rock_comp_props_;
        std::shared_ptr<Wells> wells_;
        const double* gravity_;
        // Solvers
        const DerivedGeology &geo_;
        FullyImplicitBlackoilSolver<Grid> solver_;
        // Misc. data
        RateConverterType rateConverter_;
        std::vector<int> allcells_;

        void
        computeRESV(const std::size_t               step,
                    const BlackoilState&            x,
                    WellStateFullyImplicitBlackoil& xw);
    };




    template<class T>
    SimulatorFullyImplicitBlackoil<T>::SimulatorFullyImplicitBlackoil(const parameter::ParameterGroup& param,
                                                                      ScheduleConstPtr schedule,
                                                                   const Grid& grid,
                                                                   const DerivedGeology& geo,
                                                                   BlackoilPropsAdInterface& props,
                                                                   const RockCompressibility* rock_comp_props,
                                                                   const Wells* wells,
                                                                   NewtonIterationBlackoilInterface& linsolver,
                                                                   const double* gravity,
                                                                   const bool has_disgas,
                                                                   const bool has_vapoil )

    {
        pimpl_.reset(new Impl(param, schedule, grid, geo, props, rock_comp_props, wells, linsolver, gravity, has_disgas, has_vapoil));
    }





    template<class T>
    SimulatorReport SimulatorFullyImplicitBlackoil<T>::run(SimulatorTimer& timer,
                                                        BlackoilState& state,
                                                        WellStateFullyImplicitBlackoil& well_state)
    {
        return pimpl_->run(timer, state, well_state);
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
                                                  ScheduleConstPtr schedule,
                                               const Grid& grid,
                                               const DerivedGeology& geo,
                                               BlackoilPropsAdInterface& props,
                                               const RockCompressibility* rock_comp_props,
                                               const Wells* wells,
                                               NewtonIterationBlackoilInterface& linsolver,
                                               const double* gravity,
                                               const bool has_disgas,
                                               const bool has_vapoil)
        : schedule_(schedule)
        , grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          wells_(clone_wells(wells), & destroy_wells),
          gravity_(gravity),
          geo_(geo),
          solver_(param, grid_, props_, geo_, rock_comp_props, *wells_, linsolver, has_disgas, has_vapoil)
        , rateConverter_(props_, std::vector<int>(AutoDiffGrid::numCells(grid_), 0))
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
                                                              BlackoilState& state,
                                                              WellStateFullyImplicitBlackoil& well_state)
    {
        // Initialisation.
        std::vector<double> porevol;
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(AutoDiffGrid::numCells(grid_), AutoDiffGrid::beginCellVolumes(grid_), props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
        } else {
            computePorevolume(AutoDiffGrid::numCells(grid_), AutoDiffGrid::beginCellVolumes(grid_), props_.porosity(), porevol);
        }
        // const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        std::vector<double> initial_porevol = porevol;

        // Main simulation loop.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::vector<double> fractional_flows;
        std::vector<double> well_resflows_phase;
        if (wells_) {
            well_resflows_phase.resize((wells_->number_of_phases)*(wells_->number_of_wells), 0.0);
        }
        std::fstream tstep_os;
        if (output_) {
            std::string filename = output_dir_ + "/step_timing.param";
            tstep_os.open(filename.c_str(), std::fstream::out | std::fstream::app);
        }
        while (!timer.done()) {
            // Report timestep and (optionally) write state to disk.
            step_timer.start();
            timer.report(std::cout);
            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
                outputWellStateMatlab(well_state,timer.currentStepNum(), output_dir_);

            }

            SimulatorReport sreport;
            {
                computeRESV(timer.currentStepNum(), state, well_state);

                // Run solver.
                solver_timer.start();
                solver_.step(timer.currentStepLength(), state, well_state);

                // Stop timer and report.
                solver_timer.stop();
                const double st = solver_timer.secsSinceStart();
                std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;

                stime += st;
                sreport.pressure_time = st;
            }

            // Update pore volumes if rock is compressible.
            if (rock_comp_props_ && rock_comp_props_->isActive()) {
                initial_porevol = porevol;
                computePorevolume(AutoDiffGrid::numCells(grid_), AutoDiffGrid::beginCellVolumes(grid_), props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
            }

            // Hysteresis
            props_.updateSatHyst(state.saturation(), allcells_);

            sreport.total_time =  step_timer.secsSinceStart();
            if (output_) {
                sreport.reportParam(tstep_os);

                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
                outputWellStateMatlab(well_state,timer.currentStepNum(), output_dir_);
                tstep_os.close();
            }

            // advance to next timestep before reporting at this location
            // ++timer; // Commented out since this has temporarily moved to the main() function.
            break; // this is a temporary measure
        }

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
                      const BlackoilState&            x,
                      WellStateFullyImplicitBlackoil& xw)
    {
        typedef SimFIBODetails::WellMap WellMap;

        const std::vector<WellConstPtr>& w_ecl = schedule_->getWells(step);
        const WellMap& wmap = SimFIBODetails::mapWells(w_ecl);

        const std::vector<int>& resv_prod =
            SimFIBODetails::resvProducers(*wells_, step, wmap);

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
                WellControls* ctrl = wells_->ctrls[*rp];

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
                if (wells_->name[*rp] != 0) {
                    WellMap::const_iterator i = wmap.find(wells_->name[*rp]);

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
