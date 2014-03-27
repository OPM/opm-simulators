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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>
#include <iostream>


namespace Opm
{

    class SimulatorFullyImplicitBlackoil::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             BlackoilPropsAdInterface& props,
             const RockCompressibility* rock_comp_props,
             WellsManager& wells_manager,
             LinearSolverInterface& linsolver,
             const double* gravity);

        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state,
                            WellStateFullyImplicitBlackoil& well_state);

    private:
        // Data.

        // Parameters for output.
        bool output_;
        bool output_vtk_;
        std::string output_dir_;
        int output_interval_;
        // Parameters for well control
        bool check_well_controls_;
        int max_well_control_iterations_;
        // Observed objects.
        const UnstructuredGrid& grid_;
        BlackoilPropsAdInterface& props_;
        const RockCompressibility* rock_comp_props_;
        WellsManager& wells_manager_;
        const Wells* wells_;
        const double* gravity_;
        // Solvers
        DerivedGeology geo_;
        FullyImplicitBlackoilSolver solver_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorFullyImplicitBlackoil::SimulatorFullyImplicitBlackoil(const parameter::ParameterGroup& param,
                                                                   const UnstructuredGrid& grid,
                                                                   BlackoilPropsAdInterface& props,
                                                                   const RockCompressibility* rock_comp_props,
                                                                   WellsManager& wells_manager,
                                                                   LinearSolverInterface& linsolver,
                                                                   const double* gravity)

    {
        pimpl_.reset(new Impl(param, grid, props, rock_comp_props, wells_manager, linsolver, gravity));
    }





    SimulatorReport SimulatorFullyImplicitBlackoil::run(SimulatorTimer& timer,
                                                        BlackoilState& state,
                                                        WellStateFullyImplicitBlackoil& well_state)
    {
        return pimpl_->run(timer, state, well_state);
    }



    static void outputStateVtk(const UnstructuredGrid& grid,
                               const Opm::BlackoilState& state,
                               const int step,
                               const std::string& output_dir)
    {
        // Write data in VTK format.
        std::ostringstream vtkfilename;
        vtkfilename << output_dir << "/vtk_files";
        boost::filesystem::path fpath(vtkfilename.str());
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        vtkfilename << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        if (!vtkfile) {
            OPM_THROW(std::runtime_error, "Failed to open " << vtkfilename.str());
        }
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        std::vector<double> cell_velocity;
        Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
        dm["velocity"] = &cell_velocity;
        Opm::writeVtkData(grid, dm, vtkfile);
    }


    static void outputStateMatlab(const UnstructuredGrid& grid,
                                  const Opm::BlackoilState& state,
                                  const int step,
                                  const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        dm["surfvolume"] = &state.surfacevol();
        std::vector<double> cell_velocity;
        Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
        dm["velocity"] = &cell_velocity;

        // Write data (not grid) in Matlab format
        for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
            std::ostringstream fname;
            fname << output_dir << "/" << it->first;
            boost::filesystem::path fpath = fname.str();
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }
            fname << "/" << std::setw(3) << std::setfill('0') << step << ".txt";
            std::ofstream file(fname.str().c_str());
            if (!file) {
                OPM_THROW(std::runtime_error, "Failed to open " << fname.str());
            }
            file.precision(15);
            const std::vector<double>& d = *(it->second);
            std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
        }
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
    SimulatorFullyImplicitBlackoil::Impl::Impl(const parameter::ParameterGroup& param,
                                               const UnstructuredGrid& grid,
                                               BlackoilPropsAdInterface& props,
                                               const RockCompressibility* rock_comp_props,
                                               WellsManager& wells_manager,
                                               LinearSolverInterface& linsolver,
                                               const double* gravity)
        : grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          wells_manager_(wells_manager),
          wells_(wells_manager.c_wells()),
          gravity_(gravity),
          geo_(grid_, props_, gravity_),
          solver_(grid_, props_, geo_, rock_comp_props, *wells_manager.c_wells(), linsolver)
          /*                   param.getDefault("nl_pressure_residual_tolerance", 0.0),
                               param.getDefault("nl_pressure_change_tolerance", 1.0),
                               param.getDefault("nl_pressure_maxiter", 10),
                               gravity,  */
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
        const int num_cells = grid.number_of_cells;
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
    }




    SimulatorReport SimulatorFullyImplicitBlackoil::Impl::run(SimulatorTimer& timer,
                                                              BlackoilState& state,
                                                              WellStateFullyImplicitBlackoil& well_state)
    {
        // Initialisation.
        std::vector<double> porevol;
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
        } else {
            computePorevolume(grid_, props_.porosity(), porevol);
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

            // Solve pressure equation.
            // if (check_well_controls_) {
            //     computeFractionalFlow(props_, allcells_,
            //                           state.pressure(), state.surfacevol(), state.saturation(),
            //                           fractional_flows);
            //     wells_manager_.applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
            // }

            bool well_control_passed = !check_well_controls_;
            int well_control_iteration = 0;
            do {
                // Run solver.
                solver_timer.start();
                std::vector<double> initial_pressure = state.pressure();
                solver_.step(timer.currentStepLength(), state, well_state);

                // Stop timer and report.
                solver_timer.stop();
                const double st = solver_timer.secsSinceStart();
                std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;

                stime += st;
                sreport.pressure_time = st;

                // Optionally, check if well controls are satisfied.
                if (check_well_controls_) {
                    std::cout << "Checking well conditions." << std::endl;
                    // For testing we set surface := reservoir
                    well_control_passed = wells_manager_.conditionsMet(well_state.bhp(), well_state.wellRates(), well_state.wellRates());
                    ++well_control_iteration;
                    if (!well_control_passed && well_control_iteration > max_well_control_iterations_) {
                        OPM_THROW(std::runtime_error, "Could not satisfy well conditions in " << max_well_control_iterations_ << " tries.");
                    }
                    if (!well_control_passed) {
                        std::cout << "Well controls not passed, solving again." << std::endl;
                    } else {
                        std::cout << "Well conditions met." << std::endl;
                    }
                }
            } while (!well_control_passed);

            // Update pore volumes if rock is compressible.
            if (rock_comp_props_ && rock_comp_props_->isActive()) {
                initial_porevol = porevol;
                computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
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
            ++timer;
            break; // this is a temporary measure
        }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }


} // namespace Opm
