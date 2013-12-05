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

#include <opm/autodiff/polymer/SimulatorFullyImplicitTwophase.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/polymer/FullyImplicitTwoPhaseSolver.hpp>
#include <opm/autodiff/polymer/IncompPropsAdInterface.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>


#include <opm/core/utility/miscUtilities.hpp>

#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/core/simulator/TwophaseState.hpp>

#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>
#include <iostream>
#include <Eigen/Eigen>
namespace Opm
{

    class SimulatorFullyImplicitTwophase::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             const IncompPropsAdInterface& props,
             LinearSolverInterface& linsolver,
             std::vector<double>& src);

        SimulatorReport run(SimulatorTimer& timer,
                            TwophaseState& state,
                            std::vector<double>& src);

    private:

        // Parameters for output.
        bool output_;
        bool output_vtk_;
        std::string output_dir_;
        int output_interval_;
        // Parameters for well control
        // Observed objects.
        const UnstructuredGrid& grid_;
        const IncompPropsAdInterface& props_;
        // Solvers
        FullyImplicitTwoPhaseSolver solver_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorFullyImplicitTwophase::SimulatorFullyImplicitTwophase(const parameter::ParameterGroup& param,
                                                                   const UnstructuredGrid& grid,
                                                                   const IncompPropsAdInterface& props,
                                                                   LinearSolverInterface& linsolver,
                                                                   std::vector<double>& src)
    {
        pimpl_.reset(new Impl(param, grid, props, linsolver, src));
    }





    SimulatorReport SimulatorFullyImplicitTwophase::run(SimulatorTimer& timer,
                                                        TwophaseState& state,
                                                        std::vector<double>& src)
    {
        return pimpl_->run(timer, state, src);
    }



    static void outputStateVtk(const UnstructuredGrid& grid,
                               const Opm::TwophaseState& state,
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
/*
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
    
   */ 
    
    
    
    
    SimulatorFullyImplicitTwophase::Impl::Impl(const parameter::ParameterGroup& param,
                                               const UnstructuredGrid& grid,
                                               const IncompPropsAdInterface& props,
                                               LinearSolverInterface& linsolver,
                                               std::vector<double>& src)
        : grid_(grid),
          props_(props),
          solver_(grid_, props_, linsolver)
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
        const int num_cells = grid.number_of_cells;
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
    }




    SimulatorReport SimulatorFullyImplicitTwophase::Impl::run(SimulatorTimer& timer,
                                                              TwophaseState& state,
                                                              std::vector<double>& src)
    {
        // Initialisation.
        std::vector<double> porevol;
        computePorevolume(grid_, props_.porosity(), porevol);
        // const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        std::vector<double> initial_porevol = porevol;

        // Main simulation loop.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        // These must be changed for three-phase.
        std::vector<double> fractional_flows;
        std::fstream tstep_os;
        if (output_) {
            std::string filename = output_dir_ + "/step_timing.param";
            tstep_os.open(filename.c_str(), std::fstream::out | std::fstream::app);
        }
        for (; !timer.done(); ++timer) {
            // Report timestep and (optionally) write state to disk.
            step_timer.start();
            timer.report(std::cout);
            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
            }
   
            SimulatorReport sreport;

            // Run solver.
            solver_timer.start();
            std::vector<double> initial_pressure = state.pressure();
            solver_.step(timer.currentStepLength(), state, src);

            // Stop timer and report.
            solver_timer.stop();
            const double st = solver_timer.secsSinceStart();
            std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;
            stime += st;
            sreport.pressure_time = st;

            sreport.total_time =  step_timer.secsSinceStart();
            if (output_) {
                sreport.reportParam(tstep_os);
            }
        }

        if (output_) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
     //       outputWaterCut(watercut, output_dir_);
            tstep_os.close();
        }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }


} // namespace Opm
