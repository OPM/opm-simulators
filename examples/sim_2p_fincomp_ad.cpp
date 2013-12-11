#include "config.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/polymer/fullyimplicit/FullyImplicitTwoPhaseSolver.hpp>
#include <opm/polymer/fullyimplicit/IncompPropsAdBasic.hpp>


#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

int main (int argc, char** argv)
try
{
    using namespace Opm;
    parameter::ParameterGroup param(argc, argv, false);
    int nx = param.getDefault("nx", 30);
    int ny = param.getDefault("ny", 1);
    int nz = 1;
    double dx = 10./nx;
    double dy = 1.0;
    double dz = 1.0;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    std::vector<double> viscosity(num_phases, 1.0*centi*Poise);
    viscosity[0] = 0.5 * centi * Poise;
    viscosity[1] = 5 * centi * Poise;
    double porosity = 0.35;
    double permeability = 10.0*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Linear;
    IncompPropsAdBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);
    std::vector<double> omega;
    std::vector<double> src(num_cells, 0.0);
    src[0] = 1. / day;
    src[num_cells-1] = -1. / day;

    FlowBCManager bcs;
    LinearSolverUmfpack linsolver;
    FullyImplicitTwoPhaseSolver solver(grid, props, linsolver);
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    const double dt = param.getDefault("dt", 10.) * day;
    const int num_time_steps = param.getDefault("nsteps", 10);
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }
    TwophaseState state;
    state.init(grid, 2);

    //initial sat
    for (int c = 0; c < num_cells; ++c) {
        state.saturation()[2*c] = 0.2;
        state.saturation()[2*c+1] = 0.8;
    }
    std::vector<double> p(num_cells, 100*Opm::unit::barsa);
    state.pressure() = p;
    std::ostringstream vtkfilename;

    // Write the initial state.
    vtkfilename.str("");
    vtkfilename << "sim_2p_fincomp_" << std::setw(3) << std::setfill('0') << 0 << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    Opm::DataMap dm;
    dm["saturation"] = &state.saturation();
    dm["pressure"] = &state.pressure();
    Opm::writeVtkData(grid, dm, vtkfile);
    for (int i = 0; i < num_time_steps; ++i) {
        solver.step(dt, state, src);
        vtkfilename.str("");
        vtkfilename << "sim_2p_fincomp_" << std::setw(3) << std::setfill('0') << i + 1 << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        Opm::writeVtkData(grid, dm, vtkfile);
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
