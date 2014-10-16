/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/polymer/fullyimplicit/SimulatorFullyImplicitCompressiblePolymer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/polymer/fullyimplicit/FullyImplicitCompressiblePolymerSolver.hpp>
#include <opm/polymer/fullyimplicit/utilities.hpp>
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/wells/WellsManager.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>
#include <iostream>


namespace Opm
{



    namespace 
    {
        static void outputStateVtk(const UnstructuredGrid& grid,
                                   const Opm::PolymerBlackoilState& state,
                                   const int step,
                                   const std::string& output_dir);
        static void outputStateMatlab(const UnstructuredGrid& grid,
                                      const Opm::PolymerBlackoilState& state,
                                      const int step,
                                      const std::string& output_dir);
        static  void outputWaterCut(const Opm::Watercut& watercut,
                	            const std::string& output_dir);
    } // anonymous namespace

    class SimulatorFullyImplicitCompressiblePolymer::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             const DerivedGeology& geo,
             const BlackoilPropsAdInterface& props,
             const PolymerPropsAd&          polymer_props,
             const RockCompressibility* rock_comp_props,
             WellsManager& wells_manager,
             PolymerInflowInterface& polymer_inflow,
             NewtonIterationBlackoilInterface& linsolver,
             const double* gravity);

        SimulatorReport run(SimulatorTimer& timer,
                            PolymerBlackoilState& state,
                            WellState& well_state);

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
        const BlackoilPropsAdInterface& props_;
        const PolymerPropsAd&       polymer_props_;
        const RockCompressibility* rock_comp_props_;
        WellsManager& wells_manager_;
        const Wells* wells_;
        PolymerInflowInterface& polymer_inflow_;
        const double* gravity_;
        // Solvers
        DerivedGeology geo_;
        FullyImplicitCompressiblePolymerSolver solver_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorFullyImplicitCompressiblePolymer::
    SimulatorFullyImplicitCompressiblePolymer(const parameter::ParameterGroup& param,
                                                                   const UnstructuredGrid& grid,
                                                          const DerivedGeology& geo,
                                                                   const BlackoilPropsAdInterface& props,
                                                                   const PolymerPropsAd&    polymer_props,
                                                                   const RockCompressibility* rock_comp_props,
                                                                   WellsManager& wells_manager,
                                                                   PolymerInflowInterface&  polymer_inflow,
                                                                   NewtonIterationBlackoilInterface& linsolver,
                                                                   const double* gravity)

    {
        pimpl_.reset(new Impl(param, grid, geo, props, polymer_props, rock_comp_props, wells_manager, polymer_inflow, linsolver, gravity));
    }





    SimulatorReport SimulatorFullyImplicitCompressiblePolymer::run(SimulatorTimer& timer,
                                                        PolymerBlackoilState& state,
                                                        WellState& well_state)
    {
        return pimpl_->run(timer, state, well_state);
    }





    // \TODO: Treat bcs.
    SimulatorFullyImplicitCompressiblePolymer::Impl::Impl(const parameter::ParameterGroup& param,
                   			                              const UnstructuredGrid& grid,
                                                          const DerivedGeology& geo,
                              			                  const BlackoilPropsAdInterface& props,
                                         			      const PolymerPropsAd&    polymer_props,
                                  			              const RockCompressibility* rock_comp_props,
                                  		                  WellsManager& wells_manager,
                                 			              PolymerInflowInterface&  polymer_inflow,
                         			                      NewtonIterationBlackoilInterface& linsolver,
                                    		              const double* gravity)
        : grid_(grid),
          props_(props),
          polymer_props_(polymer_props),
          rock_comp_props_(rock_comp_props),
          wells_manager_(wells_manager),
          wells_(wells_manager.c_wells()),
          polymer_inflow_(polymer_inflow),
          gravity_(gravity),
          geo_(geo),
          solver_(grid_, props_, geo_, rock_comp_props, polymer_props, *wells_manager.c_wells(), linsolver)

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




    SimulatorReport SimulatorFullyImplicitCompressiblePolymer::Impl::run(SimulatorTimer& timer,
                                                              PolymerBlackoilState& state,
                                                              WellState& well_state)
    {

        // Initialisation.
        std::vector<double> porevol;
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
        } else {
            computePorevolume(grid_, props_.porosity(), porevol);
        }
        const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        std::vector<double> initial_porevol = porevol;

        std::vector<double> polymer_inflow_c(grid_.number_of_cells);
		std::vector<double> transport_src(grid_.number_of_cells);
        // Main simulation loop.
        Opm::time::StopWatch solver_timer;
        double stime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        double tot_injected[2] = { 0.0 };
        double tot_produced[2] = { 0.0 };
        Opm::Watercut watercut;
        watercut.push(0.0, 0.0, 0.0);
#if 0
        // These must be changed for three-phase.
        double init_surfvol[2] = { 0.0 };
        double inplace_surfvol[2] = { 0.0 };
        Opm::computeSaturatedVol(porevol, state.surfacevol(), init_surfvol);
        Opm::WellReport wellreport;
#endif
        std::vector<double> fractional_flows;
        std::vector<double> well_resflows_phase;
        if (wells_) {
            well_resflows_phase.resize((wells_->number_of_phases)*(wells_->number_of_wells), 0.0);
#if 0
            wellreport.push(props_, *wells_,
                            state.pressure(), state.surfacevol(), state.saturation(),
                            0.0, well_state.bhp(), well_state.perfRates());
#endif
        }
        std::fstream tstep_os;
        if (output_) {
            std::string filename = output_dir_ + "/step_timing.param";
            tstep_os.open(filename.c_str(), std::fstream::out | std::fstream::app);
        }
//        while (!timer.done()) {
            // Report timestep and (optionally) write state to disk.
            step_timer.start();
            timer.report(std::cout);
            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
     //           outputWellStateMatlab(well_state,timer.currentStepNum(), output_dir_);

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
            	// Process transport sources (to include bdy terms and well flows).
//            	Opm::computeTransportSource(props_, wells_, well_state, transport_src);
                // Run solver.
                const double current_time = timer.simulationTimeElapsed();
                double stepsize = timer.currentStepLength();
                polymer_inflow_.getInflowValues(current_time, current_time + stepsize, polymer_inflow_c);
                solver_timer.start();
                std::vector<double> initial_pressure = state.pressure();
            	solver_.step(timer.currentStepLength(), state, well_state, polymer_inflow_c, transport_src);
                // Stop timer and report.
                solver_timer.stop();
                const double st = solver_timer.secsSinceStart();
                std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;

                stime += st;
                sreport.pressure_time = st;

                // Optionally, check if well controls are satisfied.
                if (check_well_controls_) {
                    Opm::computePhaseFlowRatesPerWell(*wells_,
                                                      well_state.perfRates(),
                                                      fractional_flows,
                                                      well_resflows_phase);
                    std::cout << "Checking well conditions." << std::endl;
                    // For testing we set surface := reservoir
                    well_control_passed = wells_manager_.conditionsMet(well_state.bhp(), well_resflows_phase, well_resflows_phase);
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

         	double injected[2] = { 0.0 };
          	double produced[2] = { 0.0 };
			double polyinj = 0;
			double polyprod = 0;

            Opm::computeInjectedProduced(props_, polymer_props_,
                                         state,
                                         transport_src, polymer_inflow_c, timer.currentStepLength(),
                                         injected, produced,
                                         polyinj, polyprod);
            tot_injected[0] += injected[0];
            tot_injected[1] += injected[1];
            tot_produced[0] += produced[0];
            tot_produced[1] += produced[1];
         	watercut.push(timer.simulationTimeElapsed() + timer.currentStepLength(),
                          	  produced[0]/(produced[0] + produced[1]),
                          	  tot_produced[0]/tot_porevol_init);
            std::cout.precision(5);
            const int width = 18;
            std::cout << "\nMass balance report.\n";
            std::cout << "    Injected reservoir volumes:      "
                      << std::setw(width) << injected[0]
                      << std::setw(width) << injected[1] << std::endl;
            std::cout << "    Produced reservoir volumes:      "
                      << std::setw(width) << produced[0]
                      << std::setw(width) << produced[1] << std::endl;
            std::cout << "    Total inj reservoir volumes:     "
                      << std::setw(width) << tot_injected[0]
                      << std::setw(width) << tot_injected[1] << std::endl;
            std::cout << "    Total prod reservoir volumes:    "
                      << std::setw(width) << tot_produced[0]
                      << std::setw(width) << tot_produced[1] << std::endl;
            sreport.total_time =  step_timer.secsSinceStart();
            if (output_) {
                sreport.reportParam(tstep_os);

                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
                outputWaterCut(watercut, output_dir_);
                tstep_os.close();
            }

            // advance to next timestep before reporting at this location
      //      ++timer;

    //    }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = stime;
        report.transport_time = 0.0;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }



    namespace
    {

        static void outputStateVtk(const UnstructuredGrid& grid,
                                   const Opm::PolymerBlackoilState& state,
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
            dm["cmax"] = &state.maxconcentration();
            dm["concentration"] = &state.concentration();
            std::vector<double> cell_velocity;
            Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
            dm["velocity"] = &cell_velocity;
            Opm::writeVtkData(grid, dm, vtkfile);
        }


        static void outputStateMatlab(const UnstructuredGrid& grid,
                                      const Opm::PolymerBlackoilState& state,
                                      const int step,
                                      const std::string& output_dir)
        {
            Opm::DataMap dm;
            dm["saturation"] = &state.saturation();
            dm["pressure"] = &state.pressure();
            dm["cmax"] = &state.maxconcentration();
            dm["concentration"] = &state.concentration();
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
    }

} // namespace Opm
