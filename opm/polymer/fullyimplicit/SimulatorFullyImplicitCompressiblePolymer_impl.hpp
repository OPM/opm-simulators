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
#ifndef OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED

namespace Opm
{
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
#endif
}

/// Class collecting all necessary components for a two-phase simulation.
SimulatorFullyImplicitCompressiblePolymer::
SimulatorFullyImplicitCompressiblePolymer(const parameter::ParameterGroup& param,
                                          const UnstructuredGrid& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const PolymerPropsAd&    polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity)
: BaseType(param,
           grid,
           geo,
           props,
           rock_comp_props,
           linsolver,
           gravity,
           /*disgas=*/false,
           /*vapoil=*/false,
           eclipse_state,
           output_writer,
           /*threshold_pressures_by_face=*/std::vector<double>())
    , deck_(deck)
    , polymer_props_(polymer_props)

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
    BaseType::allcells_.resize(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        BaseType::allcells_[cell] = cell;
    }
}

SimulatorReport SimulatorFullyImplicitCompressiblePolymer::run(SimulatorTimer& timer,
                                                               typename BaseType::ReservoirState& state)
{
    WellStateFullyImplicitBlackoil prev_well_state;
    // Initialisation.
    std::vector<double> porevol;
    if (BaseType::rock_comp_props_ && BaseType::rock_comp_props_->isActive()) {
        computePorevolume(BaseType::grid_,
                          BaseType::props_.porosity(),
                          *BaseType::rock_comp_props_,
                          state.pressure(),
                          porevol);
    } else {
        computePorevolume(BaseType::grid_,
                          BaseType::props_.porosity(),
                          porevol);
    }
    std::vector<double> initial_porevol = porevol;

    std::vector<double> polymer_inflow_c(BaseType::grid_.number_of_cells);
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
#if 0
        double tot_injected[2] = { 0.0 };
        double tot_produced[2] = { 0.0 };
        Opm::Watercut watercut;
        watercut.push(0.0, 0.0, 0.0);
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

        WellsManager wells_manager(BaseType::eclipse_state_,
                                   timer.currentStepNum(),
                                   Opm::UgGridHelpers::numCells(BaseType::grid_),
                                   Opm::UgGridHelpers::globalCell(BaseType::grid_),
                                   Opm::UgGridHelpers::cartDims(BaseType::grid_),
                                   Opm::UgGridHelpers::dimensions(BaseType::grid_),
                                   Opm::UgGridHelpers::cell2Faces(BaseType::grid_),
                                   Opm::UgGridHelpers::beginFaceCentroids(BaseType::grid_),
                                   BaseType::props_.permeability());
        const Wells* wells = wells_manager.c_wells();
        WellStateFullyImplicitBlackoil well_state;
        well_state.init(wells, state, prev_well_state);
        //Compute polymer inflow.
        std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
        if (deck_->hasKeyword("WPOLYMER")) {
            if (wells_manager.c_wells() == 0) {
                OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
            }
            polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, BaseType::eclipse_state_, *wells, Opm::UgGridHelpers::numCells(BaseType::grid_), timer.currentStepNum()));
        } else {
            polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                            1.0*Opm::unit::day,
                                                            0.0));
        }
        std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(BaseType::grid_));
        polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(),
                                            timer.simulationTimeElapsed() + timer.currentStepLength(),
                                            polymer_inflow_c);

        if (output_ && (timer.currentStepNum() % output_interval_ == 0)) {
            if (output_vtk_) {
                outputStateVtk(BaseType::grid_, state, timer.currentStepNum(), output_dir_);
            }
            outputStateMatlab(BaseType::grid_, state, timer.currentStepNum(), output_dir_);
        }
        if (output_) {
            if (timer.currentStepNum() == 0) {
                BaseType::output_writer_.writeInit(timer);
            }
            BaseType::output_writer_.writeTimeStep(timer, state, well_state);
        }
        // Run solver.
        solver_timer.start();
        FullyImplicitCompressiblePolymerSolver solver(BaseType::grid_, BaseType::props_, BaseType::geo_, BaseType::rock_comp_props_, polymer_props_, *wells_manager.c_wells(), BaseType::solver_);
        solver.step(timer.currentStepLength(), state, well_state, polymer_inflow_c);
        // Stop timer and report.
        solver_timer.stop();
        const double st = solver_timer.secsSinceStart();
        std::cout << "Fully implicit solver took:  " << st << " seconds." << std::endl;

        stime += st;
        // Update pore volumes if rock is compressible.
        if (BaseType::rock_comp_props_ && BaseType::rock_comp_props_->isActive()) {
            initial_porevol = porevol;
            computePorevolume(BaseType::grid_, BaseType::props_.porosity(), *BaseType::rock_comp_props_, state.pressure(), porevol);
        }
/*
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
*/
        if (output_) {
            SimulatorReport step_report;
            step_report.pressure_time = st;
            step_report.total_time =  step_timer.secsSinceStart();
            step_report.reportParam(tstep_os);
        }
        ++timer;
        prev_well_state = well_state;
    }
    // Write final simulation state.
    if (output_) {
        if (output_vtk_) {
            outputStateVtk(BaseType::grid_, state, timer.currentStepNum(), output_dir_);
        }
        outputStateMatlab(BaseType::grid_, state, timer.currentStepNum(), output_dir_);
        BaseType::output_writer_.writeTimeStep(timer, state, prev_well_state);
    }

    total_timer.stop();
    SimulatorReport report;
    report.pressure_time = stime;
    report.transport_time = 0.0;
    report.total_time = total_timer.secsSinceStart();
    return report;
}

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED
