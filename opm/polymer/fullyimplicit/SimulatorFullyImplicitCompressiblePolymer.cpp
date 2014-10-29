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
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

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

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
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
             std::shared_ptr<EclipseState> eclipse_state,
             EclipseWriter& output_writer,
             Opm::DeckConstPtr& deck,
             NewtonIterationBlackoilInterface& linsolver,
             const double* gravity);

        SimulatorReport run(SimulatorTimer& timer,
                            PolymerBlackoilState& state);

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
        std::shared_ptr<EclipseState> eclipse_state_;
        EclipseWriter& output_writer_;
        Opm::DeckConstPtr& deck_;
        NewtonIterationBlackoilInterface& linsolver_;
        const double* gravity_;
        // Solvers
        DerivedGeology geo_;
        // Misc. data
        std::vector<int> allcells_;
    };




    SimulatorFullyImplicitCompressiblePolymer::
    SimulatorFullyImplicitCompressiblePolymer(const parameter::ParameterGroup& param,
                                              const UnstructuredGrid& grid,
                                              const DerivedGeology& geo,
                                              const BlackoilPropsAdInterface& props,
                                              const PolymerPropsAd& polymer_props,
                                              const RockCompressibility* rock_comp_props,
                                              std::shared_ptr<EclipseState> eclipse_state,
                                              EclipseWriter& output_writer,
                                              Opm::DeckConstPtr& deck,
                                              NewtonIterationBlackoilInterface& linsolver,
                                              const double* gravity)

    {
        pimpl_.reset(new Impl(param, grid, geo, props, polymer_props, rock_comp_props, eclipse_state, output_writer, deck, linsolver, gravity));
    }





    SimulatorReport SimulatorFullyImplicitCompressiblePolymer::run(SimulatorTimer& timer,
                                                        PolymerBlackoilState& state)
    {
        return pimpl_->run(timer, state);
    }





    // \TODO: Treat bcs.
    SimulatorFullyImplicitCompressiblePolymer::Impl::Impl(const parameter::ParameterGroup& param,
                   			                              const UnstructuredGrid& grid,
                                                          const DerivedGeology& geo,
                              			                  const BlackoilPropsAdInterface& props,
                                         			      const PolymerPropsAd& polymer_props,
                                  			              const RockCompressibility* rock_comp_props,
                                                          std::shared_ptr<EclipseState> eclipse_state,
                                                          EclipseWriter& output_writer,
                                                          Opm::DeckConstPtr& deck,
                         			                      NewtonIterationBlackoilInterface& linsolver,
                                    		              const double* gravity)
        : grid_(grid),
          props_(props),
          polymer_props_(polymer_props),
          rock_comp_props_(rock_comp_props),
          eclipse_state_(eclipse_state),
          output_writer_(output_writer),
          deck_(deck),
          linsolver_(linsolver),
          gravity_(gravity),
          geo_(geo)
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
                                                                         PolymerBlackoilState& state)
    {
        WellStateFullyImplicitBlackoil prev_well_state;
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
        std::string tstep_filename = output_dir_ + "/step_timing.txt";
        std::ofstream tstep_os(tstep_filename.c_str());

        //Main simulation loop.
        while (!timer.done()) {
            double tot_injected[2] = { 0.0 };
            double tot_produced[2] = { 0.0 };
            Opm::Watercut watercut;
            watercut.push(0.0, 0.0, 0.0);
#if 0
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
#endif
            // Report timestep and (optionally) write state to disk.

            step_timer.start();
            timer.report(std::cout);

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
            well_state.init(wells, state.blackoilState());
            if (timer.currentStepNum() != 0) {
                // Transfer previous well state to current.
                well_state.partialCopy(prev_well_state, *wells, prev_well_state.numWells());
            }
            //Compute polymer inflow.
            std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
            if (deck_->hasKeyword("WPOLYMER")) {
                if (wells_manager.c_wells() == 0) {
                    OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
                }
                polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, *wells, Opm::UgGridHelpers::numCells(grid_)));
            } else {
                polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                                1.0*Opm::unit::day,
                                                                0.0));
            }
            std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(grid_));
            polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(),
                                                timer.simulationTimeElapsed() + timer.currentStepLength(),
                                                polymer_inflow_c);

            if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
                if (output_vtk_) {
                    outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
                }
                outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
            }
            if (output_) {
                if (timer.currentStepNum() == 0) {
                    output_writer_.writeInit(timer);
                }
                output_writer_.writeTimeStep(timer, state.blackoilState(), well_state.basicWellState());
            }
            // Run solver.
            solver_timer.start();
            FullyImplicitCompressiblePolymerSolver solver(grid_, props_, geo_, rock_comp_props_, polymer_props_, *wells_manager.c_wells(), linsolver_);
            solver.step(timer.currentStepLength(), state, well_state, polymer_inflow_c, transport_src);
            // Stop timer and report.
            solver_timer.stop();
            const double st = solver_timer.secsSinceStart();
            std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;

            stime += st;
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
            if (output_) {
                SimulatorReport step_report;
                step_report.pressure_time = st;
                step_report.total_time =  step_timer.secsSinceStart();
                step_report.reportParam(tstep_os);
                outputWaterCut(watercut, output_dir_);
            }
            ++timer;
            prev_well_state = well_state;
        }
        // Write final simulation state.
        if (output_) {
            if (output_vtk_) {
                outputStateVtk(grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(grid_, state, timer.currentStepNum(), output_dir_);
            output_writer_.writeTimeStep(timer, state.blackoilState(), prev_well_state.basicWellState());
        }

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
