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

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/polymer/PolymerState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/polymer/IncompTpfaPolymer.hpp>
#include <opm/polymer/TransportSolverTwophasePolymer.hpp>
#include <opm/polymer/PolymerProperties.hpp>

#include <boost/scoped_ptr.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>




// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    // std::cout << "\n================    Test program for single-cell solves with polymer    ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    param.disableOutput();
    // std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    boost::scoped_ptr<GridManager> grid;
    boost::scoped_ptr<IncompPropertiesInterface> props;
    PolymerState state;
    Opm::PolymerProperties poly_props;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;

    // -------- Initialising section ----------

    // Grid init.
    grid.reset(new GridManager(2, 1, 1, 1.0, 1.0, 1.0));
    // Rock and fluid init.
    props.reset(new IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
    // Init state variables (saturation and pressure).
    initStateBasic(*grid->c_grid(), *props, param, 0.0, state);
    // Init Polymer state
    if (param.has("poly_init")) {
        double poly_init = param.getDefault("poly_init", 0.0);
        for (int cell = 0; cell < grid->c_grid()->number_of_cells; ++cell) {
            double smin[2], smax[2];
            props->satRange(1, &cell, smin, smax);
            if (state.saturation()[2*cell] > 0.5*(smin[0] + smax[0])) {
                state.concentration()[cell] = poly_init;
                state.maxconcentration()[cell] = poly_init;
            } else {
                state.saturation()[2*cell + 0] = 0.;
                state.saturation()[2*cell + 1] = 1.;
                state.concentration()[cell] = 0.;
                state.maxconcentration()[cell] = 0.;
            }
        }
    }
    // Init polymer properties.
    // Setting defaults to provide a simple example case.
    double c_max = param.getDefault("c_max_limit", 5.0);
    double mix_param = param.getDefault("mix_param", 1.0);
    double rock_density = param.getDefault("rock_density", 1000.0);
    double dead_pore_vol = param.getDefault("dead_pore_vol", 0.0); // Note that we default to no dps here!
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
    // ads_vals[1] = 0.0;
    // ads_vals[2] = 0.0;
    std::vector<double> water_vel_vals(2, -1e100);
    water_vel_vals[0] = 0.0;
    water_vel_vals[1] = 10.0;
    std::vector<double> shear_vrf_vals(2, -1e100);
    shear_vrf_vals[0] = 1.0;
    shear_vrf_vals[1] = 1.0;
    poly_props.set(c_max, mix_param, rock_density, dead_pore_vol, res_factor, c_max_ads,
                   static_cast<Opm::PolymerProperties::AdsorptionBehaviour>(ads_index),
                   c_vals_visc,  visc_mult_vals, c_vals_ads, ads_vals, water_vel_vals, shear_vrf_vals);

    // Initialising src
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> src(num_cells, 0.0);
    // Compute pore volumes, in order to enable specifying injection rate
    // terms of total pore volume.
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    const double default_injection = 1.0;
    const double flow_per_sec = param.getDefault<double>("injected_porevolumes_per_sec", default_injection)
        *porevol[0];
    src[0] = flow_per_sec;
    src[num_cells - 1] = -flow_per_sec;

    // Boundary conditions.
    FlowBCManager bcs;

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Reordering solver.
    const double nl_tolerance = param.getDefault("nl_tolerance", 1e-9);
    const int nl_maxiter = param.getDefault("nl_maxiter", 30);
    Opm::TransportSolverTwophasePolymer::SingleCellMethod method;
    std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
    if (method_string == "Bracketing") {
        method = Opm::TransportSolverTwophasePolymer::Bracketing;
    } else if (method_string == "Newton") {
        method = Opm::TransportSolverTwophasePolymer::Newton;
    } else if (method_string == "Gradient") {
        method = Opm::TransportSolverTwophasePolymer::Gradient;
    } else if (method_string == "NewtonSimpleSC") {
        method = Opm::TransportSolverTwophasePolymer::NewtonSimpleSC;
    } else if (method_string == "NewtonSimpleC") {
        method = Opm::TransportSolverTwophasePolymer::NewtonSimpleC;
    } else {
        OPM_THROW(std::runtime_error, "Unknown method: " << method_string);
    }
    Opm::TransportSolverTwophasePolymer reorder_model(*grid->c_grid(), *props, poly_props,
                                             method, nl_tolerance, nl_maxiter);

    // Warn if any parameters are unused.
    // if (param.anyUnused()) {
    //     std::cout << "--------------------   Unused parameters:   --------------------\n";
    //     param.displayUsage();
    //     std::cout << "----------------------------------------------------------------" << std::endl;
    // }

    // Write parameters to file for later reference.
    param.writeParam("test_singlecellsolves.param");

    // Setting up a number of input (s, c) pairs and solving.
    // HACK warning: we manipulate the source term,
    // but the compressibility term in the solver
    // assumes that all inflow is water inflow. Therefore
    // one must zero the compressibility term in
    // TransportSolverTwophasePolymer line 365 before compiling this program.
    // (To fix this we should add proper all-phase src terms.)
    std::vector<double> transport_src = src;
    const double dt = param.getDefault("dt", 1.0);
    const int num_sats = 501;
    const int num_concs = 501;
    // Find the face between cell 0 and 1...
    const UnstructuredGrid& ug = *grid->c_grid();
    int face01 = -1;
    for (int f = 0; f < ug.number_of_faces; ++f) {
        if (ug.face_cells[2*f] == 0 && ug.face_cells[2*f+1] == 1) {
            face01 = f;
            break;
        }
    }
    if (face01 == -1) {
        OPM_THROW(std::runtime_error, "Could not find face adjacent to cells [0 1]");
    }
    state.faceflux()[face01] = src[0];
    for (int sats = 0; sats < num_sats; ++sats) {
        const double s = double(sats)/double(num_sats - 1);
        const double ff = s; // Simplified a lot...
        for (int conc = 0; conc < num_concs; ++conc) {
            const double c = poly_props.cMax()*double(conc)/double(num_concs - 1);
            std::vector<double> polymer_inflow_c(num_cells, c);
            // std::cout << "(s, c) = (" << s << ", " << c << ")\n";
            transport_src[0] = src[0]*ff;
            // Resetting the state for next run.
            state.saturation()[0] = 0.0;
            state.saturation()[1] = 0.0;
            state.concentration()[0] = 0.0;
            state.concentration()[1] = 0.0;
            state.maxconcentration()[0] = 0.0;
            state.maxconcentration()[1] = 0.0;
            reorder_model.solve(&state.faceflux()[0],
                                &porevol[0],
                                &transport_src[0],
                                &polymer_inflow_c[0],
                                dt,
                                state.saturation(),
                                state.concentration(),
                                state.maxconcentration());

#ifdef PROFILING
            // Extract residual counts.
            typedef std::list<Opm::TransportSolverTwophasePolymer::Newton_Iter> ListRes;
            const ListRes& res_counts = reorder_model.res_counts;
            double counts[2] = { 0, 0 };
            for (ListRes::const_iterator it = res_counts.begin(); it != res_counts.end(); ++it) {
                if (it->cell == 0) {
                    ++counts[it->res_s];
                }
            }
            // std::cout << "c residual count: " << counts[0] << '\n';
            // std::cout << "s residual count: " << counts[1] << '\n';
            std::cout << counts[0] << ' ' << counts[1] << ' ' << s << ' ' << c << '\n';
#endif
        }
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

