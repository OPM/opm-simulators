/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/simulator/SimulatorCompressibleTwophase.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/core/pressure/CompressibleTpfa.hpp>

#include <opm/core/grid.h>
#include <opm/core/newwells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/wells/WellsManager.hpp>

#include <opm/core/fluid/BlackoilPropertiesInterface.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/utility/ColumnExtract.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>


namespace Opm
{

    class SimulatorCompressibleTwophase::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             const BlackoilPropertiesInterface& props,
             const RockCompressibility* rock_comp_props,
             WellsManager& wells_manager,
             const std::vector<double>& src,
             const FlowBoundaryConditions* bcs,
             LinearSolverInterface& linsolver,
             const double* gravity);

        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state,
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
        // Parameters for transport solver.
        int num_transport_substeps_;
        bool use_segregation_split_;
        // Observed objects.
        const UnstructuredGrid& grid_;
        const BlackoilPropertiesInterface& props_;
        const RockCompressibility* rock_comp_props_;
        WellsManager& wells_manager_;
        const Wells* wells_;
        const std::vector<double>& src_;
        const FlowBoundaryConditions* bcs_;
        const double* gravity_;
        // Solvers
        CompressibleTpfa psolver_;
        TransportSolverCompressibleTwophaseReorder tsolver_;
        // Needed by column-based gravity segregation solver.
        std::vector< std::vector<int> > columns_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorCompressibleTwophase::SimulatorCompressibleTwophase(const parameter::ParameterGroup& param,
                                                                 const UnstructuredGrid& grid,
                                                                 const BlackoilPropertiesInterface& props,
                                                                 const RockCompressibility* rock_comp_props,
                                                                 WellsManager& wells_manager,
                                                                 const std::vector<double>& src,
                                                                 const FlowBoundaryConditions* bcs,
                                                                 LinearSolverInterface& linsolver,
                                                                 const double* gravity)
    {
        pimpl_.reset(new Impl(param, grid, props, rock_comp_props, wells_manager, src, bcs, linsolver, gravity));
    }





    SimulatorReport SimulatorCompressibleTwophase::run(SimulatorTimer& timer,
                                                       BlackoilState& state,
                                                       WellState& well_state)
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
          THROW("Creating directories failed: " << fpath);
        }
        vtkfilename << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        if (!vtkfile) {
            THROW("Failed to open " << vtkfilename.str());
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
              THROW("Creating directories failed: " << fpath);
            }
            fname << "/" << std::setw(3) << std::setfill('0') << step << ".txt";
            std::ofstream file(fname.str().c_str());
            if (!file) {
                THROW("Failed to open " << fname.str());
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
            THROW("Failed to open " << fname);
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
            THROW("Failed to open " << fname);
        }
        wellreport.write(os);
    }



    // \TODO: make CompressibleTpfa take src and bcs.
    SimulatorCompressibleTwophase::Impl::Impl(const parameter::ParameterGroup& param,
                                              const UnstructuredGrid& grid,
                                              const BlackoilPropertiesInterface& props,
                                              const RockCompressibility* rock_comp_props,
                                              WellsManager& wells_manager,
                                              const std::vector<double>& src,
                                              const FlowBoundaryConditions* bcs,
                                              LinearSolverInterface& linsolver,
                                              const double* gravity)
        : grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          wells_manager_(wells_manager),
          wells_(wells_manager.c_wells()),
          src_(src),
          bcs_(bcs),
          gravity_(gravity),
          psolver_(grid, props, rock_comp_props, linsolver,
                   param.getDefault("nl_pressure_residual_tolerance", 0.0),
                   param.getDefault("nl_pressure_change_tolerance", 1.0),
                   param.getDefault("nl_pressure_maxiter", 10),
                   gravity, wells_manager.c_wells() /*, src, bcs*/),
          tsolver_(grid, props,
                   param.getDefault("nl_tolerance", 1e-9),
                   param.getDefault("nl_maxiter", 30))
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
                THROW("Creating directories failed: " << fpath);
            }
            output_interval_ = param.getDefault("output_interval", 1);
        }

        // Well control related init.
        check_well_controls_ = param.getDefault("check_well_controls", false);
        max_well_control_iterations_ = param.getDefault("max_well_control_iterations", 10);

        // Transport related init.
        num_transport_substeps_ = param.getDefault("num_transport_substeps", 1);
        use_segregation_split_ = param.getDefault("use_segregation_split", false);
        if (gravity != 0 && use_segregation_split_){
            tsolver_.initGravity(gravity);
            extractColumn(grid_, columns_);
        }

        // Misc init.
        const int num_cells = grid.number_of_cells;
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
    }




    SimulatorReport SimulatorCompressibleTwophase::Impl::run(SimulatorTimer& timer,
                                                             BlackoilState& state,
                                                             WellState& well_state)
    {
        std::vector<double> transport_src;

        // Initialisation.
        std::vector<double> porevol;
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            computePorevolume(grid_, props_.porosity(), *rock_comp_props_, state.pressure(), porevol);
        } else {
            computePorevolume(grid_, props_.porosity(), porevol);
        }
        const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        std::vector<double> initial_porevol = porevol;

        // Main simulation loop.
        Opm::time::StopWatch pressure_timer;
        double ptime = 0.0;
        Opm::time::StopWatch transport_timer;
        double ttime = 0.0;
        Opm::time::StopWatch step_timer;
        Opm::time::StopWatch total_timer;
        total_timer.start();
        double init_surfvol[2] = { 0.0 };
        double inplace_surfvol[2] = { 0.0 };
        double tot_injected[2] = { 0.0 };
        double tot_produced[2] = { 0.0 };
        Opm::computeSaturatedVol(porevol, state.surfacevol(), init_surfvol);
        Opm::Watercut watercut;
        watercut.push(0.0, 0.0, 0.0);
        Opm::WellReport wellreport;
        std::vector<double> fractional_flows;
        std::vector<double> well_resflows_phase;
        if (wells_) {
            well_resflows_phase.resize((wells_->number_of_phases)*(wells_->number_of_wells), 0.0);
            wellreport.push(props_, *wells_,
                            state.pressure(), state.surfacevol(), state.saturation(),
                            0.0, well_state.bhp(), well_state.perfRates());
        }
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
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
            }

            SimulatorReport sreport;

            // Solve pressure equation.
            if (check_well_controls_) {
                computeFractionalFlow(props_, allcells_,
                                      state.pressure(), state.surfacevol(), state.saturation(),
                                      fractional_flows);
                wells_manager_.applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
            }
            bool well_control_passed = !check_well_controls_;
            int well_control_iteration = 0;
            do {
                // Run solver.
                pressure_timer.start();
                std::vector<double> initial_pressure = state.pressure();
                psolver_.solve(timer.currentStepLength(), state, well_state);

                // Renormalize pressure if both fluids and rock are
                // incompressible, and there are no pressure
                // conditions (bcs or wells).  It is deemed sufficient
                // for now to renormalize using geometric volume
                // instead of pore volume.
                if (psolver_.singularPressure()) {
                    // Compute average pressures of previous and last
                    // step, and total volume.
                    double av_prev_press = 0.0;
                    double av_press = 0.0;
                    double tot_vol = 0.0;
                    const int num_cells = grid_.number_of_cells;
                    for (int cell = 0; cell < num_cells; ++cell) {
                        av_prev_press += initial_pressure[cell]*grid_.cell_volumes[cell];
                        av_press      += state.pressure()[cell]*grid_.cell_volumes[cell];
                        tot_vol       += grid_.cell_volumes[cell];
                    }
                    // Renormalization constant
                    const double ren_const = (av_prev_press - av_press)/tot_vol;
                    for (int cell = 0; cell < num_cells; ++cell) {
                        state.pressure()[cell] += ren_const;
                    }
                    const int num_wells = (wells_ == NULL) ? 0 : wells_->number_of_wells;
                    for (int well = 0; well < num_wells; ++well) {
                        well_state.bhp()[well] += ren_const;
                    }
                }

                // Stop timer and report.
                pressure_timer.stop();
                double pt = pressure_timer.secsSinceStart();
                std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
                ptime += pt;
                sreport.pressure_time = pt;

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
                        THROW("Could not satisfy well conditions in " << max_well_control_iterations_ << " tries.");
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

            // Process transport sources from well flows.
            Opm::computeTransportSource(props_, wells_, well_state, transport_src);

            // Solve transport.
            transport_timer.start();
            double stepsize = timer.currentStepLength();
            if (num_transport_substeps_ != 1) {
                stepsize /= double(num_transport_substeps_);
                std::cout << "Making " << num_transport_substeps_ << " transport substeps." << std::endl;
            }
            double injected[2] = { 0.0 };
            double produced[2] = { 0.0 };
            for (int tr_substep = 0; tr_substep < num_transport_substeps_; ++tr_substep) {
                tsolver_.solve(&state.faceflux()[0], &state.pressure()[0],
                               &initial_porevol[0], &porevol[0], &transport_src[0], stepsize,
                               state.saturation(), state.surfacevol());
                double substep_injected[2] = { 0.0 };
                double substep_produced[2] = { 0.0 };
                Opm::computeInjectedProduced(props_, state, transport_src, stepsize,
                                             substep_injected, substep_produced);
                injected[0] += substep_injected[0];
                injected[1] += substep_injected[1];
                produced[0] += substep_produced[0];
                produced[1] += substep_produced[1];
                if (gravity_ != 0 && use_segregation_split_) {
                    tsolver_.solveGravity(columns_, stepsize, state.saturation(), state.surfacevol());
                }
            }
            transport_timer.stop();
            double tt = transport_timer.secsSinceStart();
            sreport.transport_time = tt;
            std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
            ttime += tt;
            // Report volume balances.
            Opm::computeSaturatedVol(porevol, state.surfacevol(), inplace_surfvol);
            tot_injected[0] += injected[0];
            tot_injected[1] += injected[1];
            tot_produced[0] += produced[0];
            tot_produced[1] += produced[1];
            std::cout.precision(5);
            const int width = 18;
            std::cout << "\nMass balance report.\n";
            std::cout << "    Injected surface volumes:      "
                      << std::setw(width) << injected[0]
                      << std::setw(width) << injected[1] << std::endl;
            std::cout << "    Produced surface volumes:      "
                      << std::setw(width) << produced[0]
                      << std::setw(width) << produced[1] << std::endl;
            std::cout << "    Total inj surface volumes:     "
                      << std::setw(width) << tot_injected[0]
                      << std::setw(width) << tot_injected[1] << std::endl;
            std::cout << "    Total prod surface volumes:    "
                      << std::setw(width) << tot_produced[0]
                      << std::setw(width) << tot_produced[1] << std::endl;
            const double balance[2] = { init_surfvol[0] - inplace_surfvol[0] - tot_produced[0] + tot_injected[0],
                                        init_surfvol[1] - inplace_surfvol[1] - tot_produced[1] + tot_injected[1] };
            std::cout << "    Initial - inplace + inj - prod: "
                      << std::setw(width) << balance[0]
                      << std::setw(width) << balance[1]
                      << std::endl;
            std::cout << "    Relative mass error:            "
                      << std::setw(width) << balance[0]/(init_surfvol[0] + tot_injected[0])
                      << std::setw(width) << balance[1]/(init_surfvol[1] + tot_injected[1])
                      << std::endl;
            std::cout.precision(8);

            watercut.push(timer.currentTime() + timer.currentStepLength(),
                          produced[0]/(produced[0] + produced[1]),
                          tot_produced[0]/tot_porevol_init);
            if (wells_) {
                wellreport.push(props_, *wells_,
                                state.pressure(), state.surfacevol(), state.saturation(),
                                timer.currentTime() + timer.currentStepLength(),
                                well_state.bhp(), well_state.perfRates());
            }
            sreport.total_time =  step_timer.secsSinceStart();
            if (output_) {
                sreport.reportParam(tstep_os);
            }
        }

        if (output_) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
            outputWaterCut(watercut, output_dir_);
            if (wells_) {
                outputWellReport(wellreport, output_dir_);
            }
            tstep_os.close();
        }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = ptime;
        report.transport_time = ttime;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }


} // namespace Opm
