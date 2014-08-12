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
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/wells/WellsManager.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <memory>
#include <numeric>
#include <fstream>
#include <iostream>


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
             EclipseWriter& output_writer);

        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state);

    private:
        // Data.
        const parameter::ParameterGroup param_;
        // Parameters for output.
        bool output_;
        bool output_vtk_;
        std::string output_dir_;
        int output_interval_;
        // Parameters for well control
        bool check_well_controls_;
        int max_well_control_iterations_;
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
                                                                   EclipseWriter& output_writer)

    {
        pimpl_.reset(new Impl(param, grid, geo, props, rock_comp_props, linsolver, gravity, has_disgas, has_vapoil,
                              eclipse_state, output_writer));
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
                                               EclipseWriter& output_writer)
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
          output_writer_(output_writer)
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

        // Well control related init.
        check_well_controls_ = param.getDefault("check_well_controls", false);
        max_well_control_iterations_ = param.getDefault("max_well_control_iterations", 10);

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
        WellStateFullyImplicitBlackoil well_state;

        // Create timers and file for writing timing info.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        std::string tstep_filename = output_dir_ + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        // Main simulation loop.
        while (!timer.done()) {
            // Report timestep and (optionally) write state to disk.
            step_timer.start();
            timer.report(std::cout);

            WellsManager wells_manager(eclipse_state_,
                               timer.currentStepNum(),
                                                grid_,
                                 props_.permeability());

            const Wells *wells = wells_manager.c_wells();

            if (timer.currentStepNum() == 0) {
                    well_state.init(wells, state);
                    output_writer_.writeInit(timer);
            } else {
                    // TODO: add a function to update the well_state here.
            }

            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
                outputWellStateMatlab(well_state,timer.currentStepNum(), output_dir_);

            }

            SimulatorReport sreport;

            FullyImplicitBlackoilSolver<T> solver(param_, grid_, props_, geo_, rock_comp_props_, *wells, solver_, has_disgas_, has_vapoil_);

            // Run solver.
            solver_timer.start();
            std::vector<double> initial_pressure = state.pressure();
            solver.step(timer.currentStepLength(), state, well_state);

            // Stop timer and report.
            solver_timer.stop();
            const double st = solver_timer.secsSinceStart();
            std::cout << "Fully implicit solver took: " << st << " seconds." << std::endl;

            stime += st;
            sreport.pressure_time = st;

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
            }

            output_writer_.writeTimeStep(timer, state, well_state.basicWellState());

            ++timer;
        }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }


} // namespace Opm
