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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/polymer/PolymerState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/polymer/SimulatorPolymer.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/polymer/PolymerProperties.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/simulators/ensureDirectoryExists.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>


namespace
{
    void warnIfUnusedParams(const Opm::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    OpmLog::setupSimpleDefaultLogging(false, true, 10);
    std::cout << "\n================    Test program for incompressible two-phase flow with polymer    ===============\n\n";
    ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    std::unique_ptr<Deck> deck;
    boost::scoped_ptr<GridManager> grid;
    boost::scoped_ptr<IncompPropertiesInterface> props;
    boost::scoped_ptr<RockCompressibility> rock_comp;
    std::shared_ptr< EclipseState > eclipseState;
    std::shared_ptr<Schedule> schedule;
    std::unique_ptr<PolymerState> state;
    Opm::PolymerProperties poly_props;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        Opm::ParseContext parseContext({{ ParseContext::PARSE_RANDOM_SLASH , InputError::IGNORE }});
        Parser parser;
        deck.reset(new Deck(parser.parseFile(deck_filename , parseContext)));

        eclipseState.reset(new Opm::EclipseState(*deck , parseContext));
        schedule.reset( new Opm::Schedule(*deck, eclipseState->getInputGrid(), eclipseState->get3DProperties(), eclipseState->runspec(), parseContext));
        // Grid init
        grid.reset(new GridManager(eclipseState->getInputGrid()));
        {
            const UnstructuredGrid& ug_grid = *(grid->c_grid());
            // Rock and fluid init
            props.reset(new IncompPropertiesFromDeck(*deck, *eclipseState, ug_grid ));
            // check_well_controls = param.getDefault("check_well_controls", false);
            // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
            state.reset( new PolymerState(  UgGridHelpers::numCells( ug_grid ) , UgGridHelpers::numFaces( ug_grid ), 2));

            // Rock compressibility.
            rock_comp.reset(new RockCompressibility(*eclipseState));
            // Gravity.
            gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;
            // Init state variables (saturation and pressure).
            if (param.has("init_saturation")) {
                initStateBasic(ug_grid, *props, param, gravity[2], *state);
            } else {
                initStateFromDeck(ug_grid, *props, *deck, gravity[2], *state);
            }
            // Init polymer properties.
            poly_props.readFromDeck(*deck, *eclipseState);
        }
    } else {
        // Grid init.
        const int nx = param.getDefault("nx", 100);
        const int ny = param.getDefault("ny", 100);
        const int nz = param.getDefault("nz", 1);
        const double dx = param.getDefault("dx", 1.0);
        const double dy = param.getDefault("dy", 1.0);
        const double dz = param.getDefault("dz", 1.0);
        grid.reset(new GridManager(nx, ny, nz, dx, dy, dz));
        {
            const UnstructuredGrid& ug_grid = *(grid->c_grid());

            // Rock and fluid init.
            props.reset(new IncompPropertiesBasic(param, ug_grid.dimensions, UgGridHelpers::numCells( ug_grid )));;
            state.reset( new PolymerState(  UgGridHelpers::numCells( ug_grid ) , UgGridHelpers::numFaces( ug_grid ) , 2));
            // Rock compressibility.
            rock_comp.reset(new RockCompressibility(param));
            // Gravity.
            gravity[2] = param.getDefault("gravity", 0.0);
            // Init state variables (saturation and pressure).
            initStateBasic(ug_grid, *props, param, gravity[2], *state);
            // Init Polymer state
            if (param.has("poly_init")) {
                double poly_init = param.getDefault("poly_init", 0.0);
                for (int cell = 0; cell < UgGridHelpers::numCells( ug_grid ); ++cell) {
                    double smin[2], smax[2];

                    auto& saturation = state->saturation();
                    auto& concentration = state->getCellData( state->CONCENTRATION );
                    auto& max_concentration = state->getCellData( state->CMAX );

                    props->satRange(1, &cell, smin, smax);
                    if (saturation[2*cell] > 0.5*(smin[0] + smax[0])) {
                        concentration[cell] = poly_init;
                        max_concentration[cell] = poly_init;
                    } else {
                        saturation[2*cell + 0] = 0.;
                        saturation[2*cell + 1] = 1.;
                        concentration[cell] = 0.;
                        max_concentration[cell] = 0.;
                    }
                }
            }
        }
        // Init polymer properties.
        // Setting defaults to provide a simple example case.
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
    }

    // Warn if gravity but no density difference.
    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    if (use_gravity) {
        if (props->density()[0] == props->density()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Initialising src
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> src(num_cells, 0.0);
    if (use_deck) {
        // Do nothing, wells will be the driving force, not source terms.
    } else {
        // Compute pore volumes, in order to enable specifying injection rate
        // terms of total pore volume.
        std::vector<double> porevol;
        if (rock_comp->isActive()) {
            computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state->pressure(), porevol);
        } else {
            computePorevolume(*grid->c_grid(), props->porosity(), porevol);
        }
        const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        const double default_injection = use_gravity ? 0.0 : 0.1;
        const double flow_per_sec = param.getDefault<double>("injected_porevolumes_per_day", default_injection)
            *tot_porevol_init/unit::day;
        src[0] = flow_per_sec;
        src[num_cells - 1] = -flow_per_sec;
    }

    // Boundary conditions.
    FlowBCManager bcs;
    if (param.getDefault("use_pside", false)) {
        int pside = param.get<int>("pside");
        double pside_pressure = param.get<double>("pside_pressure");
        bcs.pressureSide(*grid->c_grid(), FlowBCManager::Side(pside), pside_pressure);
    }

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    if (output) {
      std::string output_dir =
        param.getDefault("output_dir", std::string("output"));
      ensureDirectoryExists(output_dir);
      param.writeParam(output_dir + "/simulation.param");
    }


    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;

    SimulatorReport rep;
    if (!use_deck) {
        // Simple simulation without a deck.
        PolymerInflowBasic polymer_inflow(param.getDefault("poly_start_days", 300.0)*Opm::unit::day,
                                          param.getDefault("poly_end_days", 800.0)*Opm::unit::day,
                                          param.getDefault("poly_amount", poly_props.cMax()));
        WellsManager wells;
        SimulatorPolymer simulator(param,
                                   *grid->c_grid(),
                                   *props,
                                   poly_props,
                                   rock_comp->isActive() ? rock_comp.get() : 0,
                                   wells,
                                   polymer_inflow,
                                   src,
                                   bcs.c_bcs(),
                                   linsolver,
                                   grav);
        SimulatorTimer simtimer;
        simtimer.init(param);
        warnIfUnusedParams(param);
        WellState well_state;
        well_state.init(0, *state);
        rep = simulator.run(simtimer, *state, well_state);
    } else {
        // With a deck, we may have more epochs etc.

        WellState well_state;
        int step = 0;
        const auto& timeMap = schedule->getTimeMap();
        SimulatorTimer simtimer;
        simtimer.init(timeMap);
        // Check for WPOLYMER presence in last epoch to decide
        // polymer injection control type.
        const bool use_wpolymer = deck->hasKeyword("WPOLYMER");
        if (use_wpolymer) {
            if (param.has("poly_start_days")) {
                OPM_MESSAGE("Warning: Using WPOLYMER to control injection since it was found in deck. "
                        "You seem to be trying to control it via parameter poly_start_days (etc.) as well.");
            }
        }
        for (size_t reportStepIdx = 0; reportStepIdx < timeMap.numTimesteps(); ++reportStepIdx) {
            simtimer.setCurrentStepNum(reportStepIdx);

            // Report on start of report step.
            std::cout << "\n\n--------------    Starting report step " << reportStepIdx << "    --------------"
                      << "\n                  (number of remaining steps: "
                      << simtimer.numSteps() - step << ")\n\n" << std::flush;

            // Create new wells, polymer inflow controls.
            WellsManager wells(*eclipseState , *schedule, reportStepIdx , *grid->c_grid());
            boost::scoped_ptr<PolymerInflowInterface> polymer_inflow;
            if (use_wpolymer) {
                if (wells.c_wells() == 0) {
                    OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
                }
                polymer_inflow.reset(new PolymerInflowFromDeck(*schedule, *wells.c_wells(), props->numCells(), simtimer.currentStepNum()));
            } else {
                polymer_inflow.reset(new PolymerInflowBasic(param.getDefault("poly_start_days", 300.0)*Opm::unit::day,
                                                            param.getDefault("poly_end_days", 800.0)*Opm::unit::day,
                                                            param.getDefault("poly_amount", poly_props.cMax())));
            }

            // @@@ HACK: we should really make a new well state and
            // properly transfer old well state to it every report step,
            // since number of wells may change etc.
            if (reportStepIdx == 0) {
                well_state.init(wells.c_wells(), *state);
            }

            // Create and run simulator.
            SimulatorPolymer simulator(param,
                                       *grid->c_grid(),
                                       *props,
                                       poly_props,
                                       rock_comp->isActive() ? rock_comp.get() : 0,
                                       wells,
                                       *polymer_inflow,
                                       src,
                                       bcs.c_bcs(),
                                       linsolver,
                                       grav);
            if (reportStepIdx == 0) {
                warnIfUnusedParams(param);
            }
            SimulatorReport epoch_rep = simulator.run(simtimer, *state, well_state);

            // Update total timing report and remember step number.
            rep += epoch_rep;
            step = simtimer.currentStepNum();
        }
    }

    std::cout << "\n\n================    End of simulation     ===============\n\n";
    rep.report(std::cout);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
