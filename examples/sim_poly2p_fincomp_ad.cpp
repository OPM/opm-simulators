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
#include <opm/polymer/fullyimplicit/FullyImplicitTwophasePolymerSolver.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/polymer/fullyimplicit/IncompPropsAdBasic.hpp>
#include <opm/polymer/PolymerState.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/polymer/PolymerProperties.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
int main (int argc, char** argv)
try
{
    using namespace Opm;
    parameter::ParameterGroup param(argc, argv, false);
    bool use_poly_deck = param.has("deck_filename");
    if (!use_poly_deck) {
        OPM_THROW(std::runtime_error, "Polymer Properties must be read from deck_filename\n");
    }
    std::string deck_filename = param.get<std::string>("deck_filename");
    EclipseGridParser deck = EclipseGridParser(deck_filename);
    int nx = param.getDefault("nx", 20);
    int ny = param.getDefault("ny", 20);
    int nz = 1;
    double dx = 2.0;
    double dy = 2.0;
    double dz = 0.5;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    int num_phases = 2;
    using namespace Opm::unit;
    using namespace Opm::prefix;
    std::vector<double> density(num_phases, 1000.0);
    std::vector<double> viscosity(num_phases, 1.0*centi*Poise);
    viscosity[0] = 0.5 * centi * Poise;
    viscosity[0] = 5 * centi * Poise;
    double porosity = 0.35;
    double permeability = 10.0*milli*darcy;
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Linear;
    IncompPropsAdBasic props(num_phases, rel_perm_func, density, viscosity,
                                porosity, permeability, grid.dimensions, num_cells);

    // Init polymer properties.
    // Setting defaults to provide a simple example case.
        PolymerProperties polymer_props(deck);
  #if 0
    if (use_poly_deck) {
    } else {
        double c_max = param.getDefault("c_max_limit", 5.0);
        double mix_param = param.getDefault("mix_param", 1.0);
        double rock_density = param.getDefault("rock_density", 1000.0);
        double dead_pore_vol = param.getDefault("dead_pore_vol", 0.15);
        double res_factor = param.getDefault("res_factor", 1.) ; // res_factor = 1 gives no change in permeability
        double c_max_ads = param.getDefault("c_max_ads", 1.);
        int ads_index = param.getDefault<int>("ads_index", Opm::PolymerProperties::NoDesorption);
        std::vector<double> c_vals_visc(2, -1e100);
        c_vals_visc[0] = 0.0;
        c_vals_visc[1] = 7.0;
        std::vector<double> visc_mult_vals(2, -1e100);
        visc_mult_vals[0] = 1.0;
        // poly_props.visc_mult_vals[1] = param.getDefault("c_max_viscmult", 30.0);
        visc_mult_vals[1] = 20.0;
        std::vector<double> c_vals_ads(3, -1e100);
        c_vals_ads[0] = 0.0;
        c_vals_ads[1] = 2.0;
        c_vals_ads[2] = 8.0;
        std::vector<double> ads_vals(3, -1e100);
        ads_vals[0] = 0.0;
        ads_vals[1] = 0.0015;
        ads_vals[2] = 0.0025;
        PolymerProperties polymer_props;
        polymer_props.set(c_max, mix_param, rock_density, dead_pore_vol, res_factor, c_max_ads,
                         static_cast<Opm::PolymerProperties::AdsorptionBehaviour>(ads_index),
                           c_vals_visc,  visc_mult_vals, c_vals_ads, ads_vals);
    }
    #endif
    PolymerPropsAd polymer_props_ad(polymer_props);
    std::vector<double> omega;
    std::vector<double> src(num_cells, 0.0);
    std::vector<double> src_polymer(num_cells);
    src[0] = 10. / day;
    src[num_cells-1] = -1. / day;

    PolymerInflowBasic polymer_inflow(param.getDefault("poly_start_days", 30.0)*Opm::unit::day,
                                      param.getDefault("poly_end_days", 80.0)*Opm::unit::day,
                                      param.getDefault("poly_amount", polymer_props.cMax()));
    FlowBCManager bcs;
    LinearSolverUmfpack linsolver;
    FullyImplicitTwophasePolymerSolver solver(grid, props,polymer_props_ad, linsolver);
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    const double dt = param.getDefault("dt", 10.) * day;
    const int num_time_steps = param.getDefault("nsteps", 10);
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }
    PolymerState state;
    state.init(grid, 2);
    //initial sat
    for (int c = 0; c < num_cells; ++c) {
        state.saturation()[2*c] = 0.2;
        state.saturation()[2*c+1] = 0.8;
    }
    std::vector<double> p(num_cells, 200*Opm::unit::barsa);
    state.pressure() = p;

    std::vector<double> c(num_cells, 0.0);
    state.concentration() = c;
    std::ostringstream vtkfilename;
    double currentime = 0;
    for (int i = 0; i < num_time_steps; ++i) {
        currentime += dt;
        polymer_inflow.getInflowValues(currentime, currentime+dt, src_polymer);
        solver.step(dt, state, src, src_polymer);
        vtkfilename.str("");
        vtkfilename << "sim_poly2p_fincomp_ad_" << std::setw(3) << std::setfill('0') << i << ".vtu";
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
