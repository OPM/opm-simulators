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

#include <opm/core/simulator/SimulatorIncompTwophase.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/core/pressure/IncompTpfa.hpp>

#include <opm/core/grid.h>
#include <opm/core/newwells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>

#include <opm/core/wells/WellsManager.hpp>

#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/utility/ColumnExtract.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
//#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/transport/ImpliciteTwoPhaseTransportSolver.hpp>
#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>


namespace Opm
{

    class SimulatorIncompTwophase::Impl
    {
    public:
        Impl(const parameter::ParameterGroup& param,
             const UnstructuredGrid& grid,
             const IncompPropertiesInterface& props,
             const RockCompressibility* rock_comp_props,
             WellsManager& wells_manager,
             const std::vector<double>& src,
             const FlowBoundaryConditions* bcs,
             LinearSolverInterface& linsolver,
             const double* gravity);

        SimulatorReport run(SimulatorTimer& timer,
                            TwophaseState& state,
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
        const IncompPropertiesInterface& props_;
        const RockCompressibility* rock_comp_props_;
        WellsManager& wells_manager_;
        const Wells* wells_;
        const std::vector<double>& src_;
        const FlowBoundaryConditions* bcs_;
        // Solvers
        IncompTpfa psolver_;
        // this should maybe be packed in a separate file
        boost::scoped_ptr<TwoPhaseTransportSolver> tsolver_;
        //ImpliciteTwoPhaseTransportSolver tsolver_;
        // Needed by column-based gravity segregation solver.
        std::vector< std::vector<int> > columns_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorIncompTwophase::SimulatorIncompTwophase(const parameter::ParameterGroup& param,
                                                     const UnstructuredGrid& grid,
                                                     const IncompPropertiesInterface& props,
                                                     const RockCompressibility* rock_comp_props,
                                                     WellsManager& wells_manager,
                                                     const std::vector<double>& src,
                                                     const FlowBoundaryConditions* bcs,
                                                     LinearSolverInterface& linsolver,
                                                     const double* gravity)
    {
        pimpl_.reset(new Impl(param, grid, props, rock_comp_props, wells_manager, src, bcs, linsolver, gravity));
    }





    SimulatorReport SimulatorIncompTwophase::run(SimulatorTimer& timer,
                                                 TwophaseState& state,
                                                 WellState& well_state)
    {
        return pimpl_->run(timer, state, well_state);
    }

    static void reportVolumes(std::ostream &os, double satvol[2], double tot_porevol_init,
                              double tot_injected[2], double tot_produced[2],
                              double injected[2], double produced[2],
                              double init_satvol[2])
    {
        std::cout.precision(5);
        const int width = 18;
        os << "\nVolume balance report (all numbers relative to total pore volume).\n";
        os << "    Saturated volumes:     "
           << std::setw(width) << satvol[0]/tot_porevol_init
           << std::setw(width) << satvol[1]/tot_porevol_init << std::endl;
        os << "    Injected volumes:      "
           << std::setw(width) << injected[0]/tot_porevol_init
           << std::setw(width) << injected[1]/tot_porevol_init << std::endl;
        os << "    Produced volumes:      "
           << std::setw(width) << produced[0]/tot_porevol_init
           << std::setw(width) << produced[1]/tot_porevol_init << std::endl;
        os << "    Total inj volumes:     "
           << std::setw(width) << tot_injected[0]/tot_porevol_init
           << std::setw(width) << tot_injected[1]/tot_porevol_init << std::endl;
        os << "    Total prod volumes:    "
           << std::setw(width) << tot_produced[0]/tot_porevol_init
           << std::setw(width) << tot_produced[1]/tot_porevol_init << std::endl;
        os << "    In-place + prod - inj: "
           << std::setw(width) << (satvol[0] + tot_produced[0] - tot_injected[0])/tot_porevol_init
           << std::setw(width) << (satvol[1] + tot_produced[1] - tot_injected[1])/tot_porevol_init << std::endl;
        os << "    Init - now - pr + inj: "
           << std::setw(width) << (init_satvol[0] - satvol[0] - tot_produced[0] + tot_injected[0])/tot_porevol_init
           << std::setw(width) << (init_satvol[1] - satvol[1] - tot_produced[1] + tot_injected[1])/tot_porevol_init
           << std::endl;
        os.precision(8);
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

    static void outputVectorMatlab(const std::string& name,
                                   const std::vector<int>& vec,
                                   const int step,
                                   const std::string& output_dir)
    {
        std::ostringstream fname;
        fname << output_dir << "/" << name;
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
        std::copy(vec.begin(), vec.end(), std::ostream_iterator<double>(file, "\n"));
    }

    static void outputStateMatlab(const UnstructuredGrid& grid,
                                  const Opm::TwophaseState& state,
                                  const int step,
                                  const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
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


    static bool allNeumannBCs(const FlowBoundaryConditions* bcs)
    {
        if (bcs == NULL) {
            return true;
        } else {
            return std::find(bcs->type, bcs->type + bcs->nbc, BC_PRESSURE)
                == bcs->type + bcs->nbc;
        }
    }


    static bool allRateWells(const Wells* wells)
    {
        if (wells == NULL) {
            return true;
        }
        const int nw = wells->number_of_wells;
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells->ctrls[w];
            if (wc->current >= 0) {
                if (wc->type[wc->current] == BHP) {
                    return false;
                }
            }
        }
        return true;
    }





    SimulatorIncompTwophase::Impl::Impl(const parameter::ParameterGroup& param,
                                        const UnstructuredGrid& grid,
                                        const IncompPropertiesInterface& props,
                                        const RockCompressibility* rock_comp_props,
                                        WellsManager& wells_manager,
                                        const std::vector<double>& src,
                                        const FlowBoundaryConditions* bcs,
                                        LinearSolverInterface& linsolverp,
                                        const double* gravity)
        : grid_(grid),
          props_(props),
          rock_comp_props_(rock_comp_props),
          wells_manager_(wells_manager),
          wells_(wells_manager.c_wells()),
          src_(src),
          bcs_(bcs),
          psolver_(grid, props, rock_comp_props, linsolverp,
                   param.getDefault("nl_pressure_residual_tolerance", 0.0),
                   param.getDefault("nl_pressure_change_tolerance", 1.0),
                   param.getDefault("nl_pressure_maxiter", 10),
                   gravity, wells_manager.c_wells(), src, bcs)
    {
        const bool use_reorder = param.getDefault("use_reorder", true);
        if(use_reorder){
            /*
            tsolver_.reset(new Opm::TransportModelTwoPhase(grid, props, rock_comp_props, linsolver,
                                                           param.getDefault("nl_pressure_residual_tolerance", 0.0),
                                                           param.getDefault("nl_pressure_change_tolerance", 1.0),
                                                           param.getDefault("nl_pressure_maxiter", 10),
                                                           gravity, wells_manager.c_wells(), src, bcs));
                                                           */
        }else{
            //Opm::ImplicitTransportDetails::NRReport  rpt;
            Opm::ImplicitTransportDetails::NRControl ctrl;
            ctrl.max_it = param.getDefault("max_it", 20);
            ctrl.verbosity = param.getDefault("verbosity", 0);
            ctrl.max_it_ls = param.getDefault("max_it_ls", 5);
            const bool guess_old_solution = param.getDefault("guess_old_solution", false);
            Opm::SimpleFluid2pWrappingProps fluid(props);
            std::vector<double> porevol;
            //if (rock_comp->isActive()) {
            //    computePorevolume(grid, props->porosity(), rock_comp, state.pressure(), porevol);
            //} else {
                computePorevolume(grid, props.porosity(), porevol);
            //}
            SinglePointUpwindTwoPhase<Opm::SimpleFluid2pWrappingProps>
                    model(fluid, grid, porevol, gravity, guess_old_solution);
            tsolver_.reset(new Opm::ImpliciteTwoPhaseTransportSolver(
                               wells_manager,
                               *rock_comp_props,
                               ctrl,
                               model,
                               grid,
                               props,
                               param));

        }

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
        // Misc init.
        const int num_cells = grid.number_of_cells;
        allcells_.resize(num_cells);
        for (int cell = 0; cell < num_cells; ++cell) {
            allcells_[cell] = cell;
        }
    }




    SimulatorReport SimulatorIncompTwophase::Impl::run(SimulatorTimer& timer,
                                                       TwophaseState& state,
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
        double init_satvol[2] = { 0.0 };
        double satvol[2] = { 0.0 };
        double tot_injected[2] = { 0.0 };
        double tot_produced[2] = { 0.0 };
        Opm::computeSaturatedVol(porevol, state.saturation(), init_satvol);
        std::cout << "\nInitial saturations are    " << init_satvol[0]/tot_porevol_init
                  << "    " << init_satvol[1]/tot_porevol_init << std::endl;
        Opm::Watercut watercut;
        watercut.push(0.0, 0.0, 0.0);
        Opm::WellReport wellreport;
        std::vector<double> fractional_flows;
        std::vector<double> well_resflows_phase;
        if (wells_) {
            well_resflows_phase.resize((wells_->number_of_phases)*(wells_->number_of_wells), 0.0);
            wellreport.push(props_, *wells_, state.saturation(), 0.0, well_state.bhp(), well_state.perfRates());
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
                computeFractionalFlow(props_, allcells_, state.saturation(), fractional_flows);
                wells_manager_.applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
            }
            bool well_control_passed = !check_well_controls_;
            int well_control_iteration = 0;
            do {
                // Run solver.
                pressure_timer.start();
                std::vector<double> initial_pressure = state.pressure();
                psolver_.solve(timer.currentStepLength(), state, well_state);

                // Renormalize pressure if rock is incompressible, and
                // there are no pressure conditions (bcs or wells).
                // It is deemed sufficient for now to renormalize
                // using geometric volume instead of pore volume.
                if ((rock_comp_props_ == NULL || !rock_comp_props_->isActive())
                    && allNeumannBCs(bcs_) && allRateWells(wells_)) {
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

            // Process transport sources (to include bdy terms and well flows).
            Opm::computeTransportSource(grid_, src_, state.faceflux(), 1.0,
                                        wells_, well_state.perfRates(), transport_src);

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
                //tsolver_.solve(&state.faceflux()[0], &initial_porevol[0], &transport_src[0],
                //              stepsize, state.saturation());
                tsolver_->solve(&transport_src[0],
                               &porevol[0],
                               stepsize,
                               state,
                               well_state);

                double substep_injected[2] = { 0.0 };
                double substep_produced[2] = { 0.0 };
                Opm::computeInjectedProduced(props_, state.saturation(), transport_src, stepsize,
                                             substep_injected, substep_produced);
                injected[0] += substep_injected[0];
                injected[1] += substep_injected[1];
                produced[0] += substep_produced[0];
                produced[1] += substep_produced[1];


                watercut.push(timer.currentTime() + timer.currentStepLength(),
                              produced[0]/(produced[0] + produced[1]),
                              tot_produced[0]/tot_porevol_init);

                if (wells_) {
                    wellreport.push(props_, *wells_, state.saturation(),
                                    timer.currentTime() + timer.currentStepLength(),
                                    well_state.bhp(), well_state.perfRates());
                }

            }
            transport_timer.stop();
            double tt = transport_timer.secsSinceStart();
            sreport.transport_time = tt;
            std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
            ttime += tt;
            // Report volume balances.
            Opm::computeSaturatedVol(porevol, state.saturation(), satvol);
            tot_injected[0] += injected[0];
            tot_injected[1] += injected[1];
            tot_produced[0] += produced[0];
            tot_produced[1] += produced[1];
            reportVolumes(std::cout,satvol, tot_porevol_init,
                          tot_injected, tot_produced,
                          injected, produced,
                          init_satvol);
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
