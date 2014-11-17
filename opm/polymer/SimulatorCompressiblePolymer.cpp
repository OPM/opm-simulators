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

#include <opm/polymer/SimulatorCompressiblePolymer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/polymer/CompressibleTpfaPolymer.hpp>

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

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/polymer/TransportSolverTwophaseCompressiblePolymer.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/polymerUtilities.hpp>

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
        void outputStateVtk(const UnstructuredGrid& grid,
                            const Opm::PolymerBlackoilState& state,
                            const int step,
                            const std::string& output_dir);
        void outputStateMatlab(const UnstructuredGrid& grid,
                               const Opm::PolymerBlackoilState& state,
                               const int step,
                               const std::string& output_dir);
        void outputWaterCut(const Opm::Watercut& watercut,
                            const std::string& output_dir);
        void outputWellReport(const Opm::WellReport& wellreport,
                              const std::string& output_dir);

    } // anonymous namespace



    class SimulatorCompressiblePolymer::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             const BlackoilPropertiesInterface& props,
             const PolymerProperties& poly_props,
             const RockCompressibility* rock_comp_props,
             WellsManager& wells_manager,
             const PolymerInflowInterface& polymer_inflow,
             const std::vector<double>& src,
             const FlowBoundaryConditions* bcs,
             LinearSolverInterface& linsolver,
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
        // Parameters for transport solver.
        int num_transport_substeps_;
        bool use_segregation_split_;
        // Observed objects.
        const UnstructuredGrid& grid_;
        const BlackoilPropertiesInterface& props_;
        const PolymerProperties& poly_props_;
        const RockCompressibility* rock_comp_props_;
        WellsManager& wells_manager_;
        const Wells* wells_;
        const PolymerInflowInterface& polymer_inflow_;
        const double* gravity_;
        // Solvers
        CompressibleTpfaPolymer psolver_;
        TransportSolverTwophaseCompressiblePolymer tsolver_;
        // Needed by column-based gravity segregation solver.
        std::vector< std::vector<int> > columns_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorCompressiblePolymer::SimulatorCompressiblePolymer(const parameter::ParameterGroup& param,
                                                               const UnstructuredGrid& grid,
                                                               const BlackoilPropertiesInterface& props,
                                                               const PolymerProperties& poly_props,
                                                               const RockCompressibility* rock_comp_props,
                                                               WellsManager& wells_manager,
                                                               const PolymerInflowInterface& polymer_inflow,
                                                               const std::vector<double>& src,
                                                               const FlowBoundaryConditions* bcs,
                                                               LinearSolverInterface& linsolver,
                                                               const double* gravity)
    {
        pimpl_.reset(new Impl(param, grid, props, poly_props, rock_comp_props,
                              wells_manager, polymer_inflow, src, bcs, linsolver, gravity));
    }




    SimulatorReport SimulatorCompressiblePolymer::run(SimulatorTimer& timer,
                                                      PolymerBlackoilState& state,
                                                      WellState& well_state)
    {
        return pimpl_->run(timer, state, well_state);
    }






    SimulatorCompressiblePolymer::Impl::Impl(const parameter::ParameterGroup& param,
                                             const UnstructuredGrid& grid,
                                             const BlackoilPropertiesInterface& props,
                                             const PolymerProperties& poly_props,
                                             const RockCompressibility* rock_comp_props,
                                             WellsManager& wells_manager,
                                             const PolymerInflowInterface& polymer_inflow,
                                             const std::vector<double>& src,
                                             const FlowBoundaryConditions* bcs,
                                             LinearSolverInterface& linsolver,
                                             const double* gravity)
        : grid_(grid),
          props_(props),
          poly_props_(poly_props),
          rock_comp_props_(rock_comp_props),
          wells_manager_(wells_manager),
          wells_(wells_manager.c_wells()),
          polymer_inflow_(polymer_inflow),
          gravity_(gravity),
          psolver_(grid, props, rock_comp_props, poly_props, linsolver,
                   param.getDefault("nl_pressure_residual_tolerance", 0.0),
                   param.getDefault("nl_pressure_change_tolerance", 1.0),
                   param.getDefault("nl_pressure_maxiter", 10),
                   gravity, wells_manager.c_wells()),
          tsolver_(grid, props, poly_props,
                   TransportSolverTwophaseCompressiblePolymer::Bracketing,
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
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }
            output_interval_ = param.getDefault("output_interval", 1);
        }

        // Well control related init.
        check_well_controls_ = param.getDefault("check_well_controls", false);
        max_well_control_iterations_ = param.getDefault("max_well_control_iterations", 10);

        // Transport related init.
        TransportSolverTwophaseCompressiblePolymer::SingleCellMethod method;
        std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
        if (method_string == "Bracketing") {
            method = Opm::TransportSolverTwophaseCompressiblePolymer::Bracketing;
        } else if (method_string == "Newton") {
            method = Opm::TransportSolverTwophaseCompressiblePolymer::Newton;
        } else {
            OPM_THROW(std::runtime_error, "Unknown method: " << method_string);
        }
        tsolver_.setPreferredMethod(method);
        num_transport_substeps_ = param.getDefault("num_transport_substeps", 1);
        use_segregation_split_ = param.getDefault("use_segregation_split", false);
        if (gravity != 0 && use_segregation_split_) {
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



    SimulatorReport SimulatorCompressiblePolymer::Impl::run(SimulatorTimer& timer,
                                                            PolymerBlackoilState& state,
                                                            WellState& well_state)
    {
        std::vector<double> transport_src(grid_.number_of_cells);
        std::vector<double> polymer_inflow_c(grid_.number_of_cells);

        // Initialisation.
        std::vector<double> initial_pressure;
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
        Opm::time::StopWatch total_timer;
        total_timer.start();
        double init_surfvol[2] = { 0.0 };
        double inplace_surfvol[2] = { 0.0 };
        double polymass = computePolymerMass(porevol, state.saturation(), state.concentration(), poly_props_.deadPoreVol());
        double polymass_adsorbed = computePolymerAdsorbed(grid_, props_, poly_props_, state, rock_comp_props_);
        double init_polymass = polymass + polymass_adsorbed;
        double tot_injected[2] = { 0.0 };
        double tot_produced[2] = { 0.0 };
        double tot_polyinj = 0.0;
        double tot_polyprod = 0.0;
        Opm::computeSaturatedVol(porevol, state.surfacevol(), init_surfvol);
        Opm::Watercut watercut;
        watercut.push(0.0, 0.0, 0.0);
        Opm::WellReport wellreport;
        std::vector<double> fractional_flows;
        std::vector<double> well_resflows_phase;
        if (wells_) {
            well_resflows_phase.resize((wells_->number_of_phases)*(wells_->number_of_wells), 0.0);
            wellreport.push(props_, *wells_, state.pressure(), state.surfacevol(),
                            state.saturation(), 0.0, well_state.bhp(), well_state.perfRates());
        }
        // Report timestep and (optionally) write state to disk.
        timer.report(std::cout);
        if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
        }

        initial_pressure = state.pressure();

        // Solve pressure equation.
        if (check_well_controls_) {
            computeFractionalFlow(props_, poly_props_, allcells_,
                                  state.pressure(), state.surfacevol(), state.saturation(),
                                  state.concentration(), state.maxconcentration(),
                                  fractional_flows);
            wells_manager_.applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
        }
        bool well_control_passed = !check_well_controls_;
        int well_control_iteration = 0;
        do {
            // Run solver
            pressure_timer.start();
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

            // Stop timer and report
            pressure_timer.stop();
            double pt = pressure_timer.secsSinceStart();
            std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
            ptime += pt;

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

        // Process transport sources (to include bdy terms and well flows).
        Opm::computeTransportSource(props_, wells_, well_state, transport_src);

        // Find inflow rate.
        const double current_time = timer.simulationTimeElapsed();
        double stepsize = timer.currentStepLength();
        polymer_inflow_.getInflowValues(current_time, current_time + stepsize, polymer_inflow_c);


        // Solve transport.
        transport_timer.start();
        if (num_transport_substeps_ != 1) {
            stepsize /= double(num_transport_substeps_);
            std::cout << "Making " << num_transport_substeps_ << " transport substeps." << std::endl;
        }
        double injected[2] = { 0.0 };
        double produced[2] = { 0.0 };
        double polyinj = 0.0;
        double polyprod = 0.0;
        for (int tr_substep = 0; tr_substep < num_transport_substeps_; ++tr_substep) {
            tsolver_.solve(&state.faceflux()[0], initial_pressure,
                           state.pressure(), &initial_porevol[0], &porevol[0],
                           &transport_src[0], &polymer_inflow_c[0], stepsize,
                           state.saturation(), state.surfacevol(),
                           state.concentration(), state.maxconcentration());
            double substep_injected[2] = { 0.0 };
            double substep_produced[2] = { 0.0 };
            double substep_polyinj = 0.0;
            double substep_polyprod = 0.0;
            Opm::computeInjectedProduced(props_, poly_props_,
                                         state,
                                         transport_src, polymer_inflow_c, stepsize,
                                         substep_injected, substep_produced,
                                         substep_polyinj, substep_polyprod);
            injected[0] += substep_injected[0];
            injected[1] += substep_injected[1];
            produced[0] += substep_produced[0];
            produced[1] += substep_produced[1];
            polyinj += substep_polyinj;
            polyprod += substep_polyprod;
            if (gravity_ != 0 && use_segregation_split_) {
                tsolver_.solveGravity(columns_, stepsize,
                                      state.saturation(), state.surfacevol(),
                                      state.concentration(), state.maxconcentration());
            }
        }
        transport_timer.stop();
        double tt = transport_timer.secsSinceStart();
        std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
        ttime += tt;

        // Report volume balances.
        Opm::computeSaturatedVol(porevol, state.surfacevol(), inplace_surfvol);
        polymass = Opm::computePolymerMass(porevol, state.saturation(), state.concentration(), poly_props_.deadPoreVol());
        polymass_adsorbed = Opm::computePolymerAdsorbed(grid_, props_, poly_props_,
                                                        state, rock_comp_props_);
        tot_injected[0] += injected[0];
        tot_injected[1] += injected[1];
        tot_produced[0] += produced[0];
        tot_produced[1] += produced[1];
        tot_polyinj += polyinj;
        tot_polyprod += polyprod;
        std::cout.precision(5);
        const int width = 18;
        std::cout << "\nMass balance:        "
            "                   water(surfvol)      oil(surfvol)       polymer(kg)\n";
        std::cout << "    In-place:                       "
                  << std::setw(width) << inplace_surfvol[0]
                  << std::setw(width) << inplace_surfvol[1]
                  << std::setw(width) << polymass << std::endl;
        std::cout << "    Adsorbed:                       "
                  << std::setw(width) << 0.0
                  << std::setw(width) << 0.0
                  << std::setw(width) << polymass_adsorbed << std::endl;
        std::cout << "    Injected:                       "
                  << std::setw(width) << injected[0]
                  << std::setw(width) << injected[1]
                  << std::setw(width) << polyinj << std::endl;
        std::cout << "    Produced:                       "
                  << std::setw(width) << produced[0]
                  << std::setw(width) << produced[1]
                  << std::setw(width) << polyprod << std::endl;
        std::cout << "    Total inj:                      "
                  << std::setw(width) << tot_injected[0]
                  << std::setw(width) << tot_injected[1]
                  << std::setw(width) << tot_polyinj << std::endl;
        std::cout << "    Total prod:                     "
                  << std::setw(width) << tot_produced[0]
                  << std::setw(width) << tot_produced[1]
                  << std::setw(width) << tot_polyprod << std::endl;
        const double balance[3] = { init_surfvol[0] - inplace_surfvol[0] - tot_produced[0] + tot_injected[0],
                                    init_surfvol[1] - inplace_surfvol[1] - tot_produced[1] + tot_injected[1],
                                    init_polymass - polymass - tot_polyprod + tot_polyinj - polymass_adsorbed };
        std::cout << "    Initial - inplace + inj - prod: "
                  << std::setw(width) << balance[0]
                  << std::setw(width) << balance[1]
                  << std::setw(width) << balance[2]
                  << std::endl;
        std::cout << "    Relative mass error:            "
                  << std::setw(width) << balance[0]/(init_surfvol[0] + tot_injected[0])
                  << std::setw(width) << balance[1]/(init_surfvol[1] + tot_injected[1])
                  << std::setw(width) << balance[2]/(init_polymass + tot_polyinj)
                  << std::endl;
        std::cout.precision(8);

        watercut.push(timer.simulationTimeElapsed() + timer.currentStepLength(),
                      produced[0]/(produced[0] + produced[1]),
                      tot_produced[0]/tot_porevol_init);
        if (wells_) {
        wellreport.push(props_, *wells_, state.pressure(), state.surfacevol(),
                        state.saturation(), timer.simulationTimeElapsed() + timer.currentStepLength(),
                        well_state.bhp(), well_state.perfRates());
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
        }

        total_timer.stop();

        SimulatorReport report;
        report.pressure_time = ptime;
        report.transport_time = ttime;
        report.total_time = total_timer.secsSinceStart();
        return report;
    }







    namespace
    {

        void outputStateVtk(const UnstructuredGrid& grid,
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
            vtkfilename << "/output-" << std::setw(5) << std::setfill('0') << step << ".vtu";
            std::ofstream vtkfile(vtkfilename.str().c_str());
            if (!vtkfile) {
                OPM_THROW(std::runtime_error, "Failed to open " << vtkfilename.str());
            }
            Opm::DataMap dm;
            dm["saturation"] = &state.saturation();
            dm["pressure"] = &state.pressure();
            dm["concentration"] = &state.concentration();
            dm["cmax"] = &state.maxconcentration();
            dm["surfvol"] = &state.surfacevol();
            std::vector<double> cell_velocity;
            Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
            dm["velocity"] = &cell_velocity;
            Opm::writeVtkData(grid, dm, vtkfile);
        }

        void outputStateMatlab(const UnstructuredGrid& grid,
                               const Opm::PolymerBlackoilState& state,
                               const int step,
                               const std::string& output_dir)
        {
            Opm::DataMap dm;
            dm["saturation"] = &state.saturation();
            dm["pressure"] = &state.pressure();
            dm["concentration"] = &state.concentration();
            dm["cmax"] = &state.maxconcentration();
            dm["surfvol"] = &state.surfacevol();
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
                fname << "/" << std::setw(5) << std::setfill('0') << step << ".txt";
                std::ofstream file(fname.str().c_str());
                if (!file) {
                    OPM_THROW(std::runtime_error, "Failed to open " << fname.str());
                }
                const std::vector<double>& d = *(it->second);
                std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
            }
        }

        void outputWaterCut(const Opm::Watercut& watercut,
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


        void outputWellReport(const Opm::WellReport& wellreport,
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


    } // anonymous namespace






} // namespace Opm
