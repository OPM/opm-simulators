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
#include <opm/autodiff/polymer/FullyImplicitTwoPhaseSolver.hpp>
#include <opm/autodiff/polymer/IncompPropsAdBasic.hpp>


#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
int main ()
try
{
    int nx = 3;
    int ny = 3;
    int nz = 1;
    double dx = 10.0;
    double dy = 10.0;
    double dz = 10.0;
    using namespace Opm;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    std::vector<double> viscosity(num_phases, 1.0*centi*Poise);
    double porosity = 0.5;
    double permeability = 10.0*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Linear;
    IncompPropsAdBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);
    std::vector<double> omega;
    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;

    FlowBCManager bcs;
    LinearSolverUmfpack linsolver;
    FullyImplicitTwoPhaseSolver solver(grid, props, linsolver);
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
//    const double tolerance = 1e-9;
//    const int max_iterations = 30;
    const double dt = 0.1*day;
    const int num_time_steps = 20;
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }
    TwophaseState state;
    state.init(grid, 2);

    //initial sat
    std::vector<double> sw(num_cells, 0.2);
    state.saturation() = sw;
    //initial pressure
    std::vector<double> p(num_cells, 4000);
    state.pressure() = p;
//    state.setFirstSat(allcells, props, TwophaseState::MinSat);
    std::ostringstream vtkfilename;

    for (int i = 0; i < num_time_steps; ++i) {
        solver.step(dt, state, src);
        vtkfilename.str("");
        vtkfilename << "sim_2p_fincomp-" << std::setw(3) << std::setfill('0') << i << ".vtu";
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
