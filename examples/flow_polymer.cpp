/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
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
#include "config.h"

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/initStateEquil.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/thresholdPressures.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>

#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/polymer/fullyimplicit/SimulatorFullyImplicitBlackoilPolymer.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>

#include <opm/parser/eclipse/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/OpmLog/StreamLog.hpp>
#include <opm/parser/eclipse/OpmLog/CounterLog.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>


namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
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

    std::cout << "\n================    Test program for fully implicit three-phase black-oil-polymer flow     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    if (!use_deck) {
        OPM_THROW(std::runtime_error, "This program must be run with an input deck. "
                  "Specify the deck with deck_filename=deckname.data (for example).");
    }
    std::shared_ptr<GridManager> grid;
    std::shared_ptr<BlackoilPropertiesInterface> props;
    std::shared_ptr<BlackoilPropsAdFromDeck> new_props;
    std::shared_ptr<RockCompressibility> rock_comp;
    PolymerBlackoilState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    std::string deck_filename = param.get<std::string>("deck_filename");

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::string output_dir;
    if (output) {
        // Create output directory if needed.
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        // Write simulation parameters.
        param.writeParam(output_dir + "/simulation.param");
    }

    std::string logFile = output_dir + "/LOGFILE.txt";
    Opm::ParserPtr parser(new Opm::Parser());
    {
        std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(logFile , Opm::Log::DefaultMessageTypes);
        std::shared_ptr<Opm::CounterLog> counterLog = std::make_shared<Opm::CounterLog>(Opm::Log::DefaultMessageTypes);

        Opm::OpmLog::addBackend( "STREAM" , streamLog );
        Opm::OpmLog::addBackend( "COUNTER" , counterLog );
    }


    Opm::DeckConstPtr deck;
    std::shared_ptr<EclipseState> eclipseState;
    try {
        deck = parser->parseFile(deck_filename);
        Opm::checkDeck(deck);
        eclipseState.reset(new Opm::EclipseState(deck));
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Failed to create valid ECLIPSESTATE object. See logfile: " << logFile << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Grid init
    std::vector<double> porv = eclipseState->getDoubleGridProperty("PORV")->getData();
    grid.reset(new GridManager(eclipseState->getEclipseGrid(), porv));
    auto &cGrid = *grid->c_grid();
    const PhaseUsage pu = Opm::phaseUsageFromDeck(deck);
    Opm::BlackoilOutputWriter outputWriter(cGrid,
                                           param,
                                           eclipseState,
                                           pu );

    // Rock and fluid init
    props.reset(new BlackoilPropertiesFromDeck(deck, eclipseState, *grid->c_grid(), param));
    new_props.reset(new BlackoilPropsAdFromDeck(deck, eclipseState, *grid->c_grid()));
    const bool polymer = deck->hasKeyword("POLYMER");
    const bool use_wpolymer = deck->hasKeyword("WPOLYMER");
    PolymerProperties polymer_props(deck, eclipseState);
    PolymerPropsAd polymer_props_ad(polymer_props);
    // check_well_controls = param.getDefault("check_well_controls", false);
    // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
    // Rock compressibility.
    rock_comp.reset(new RockCompressibility(deck, eclipseState));

    // Gravity.
    gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state.blackoilState());
        initBlackoilSurfvol(*grid->c_grid(), *props, state.blackoilState());
        enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
        if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
            const int np = props->numPhases();
            const int nc = grid->c_grid()->number_of_cells;
            for (int c = 0; c < nc; ++c) {
                state.gasoilratio()[c] = state.surfacevol()[c*np + pu.phase_pos[Gas]]
                    / state.surfacevol()[c*np + pu.phase_pos[Oil]];
            }
        }
    } else if (deck->hasKeyword("EQUIL") && props->numPhases() == 3) {
        state.init(*grid->c_grid(), props->numPhases());
        const double grav = param.getDefault("gravity", unit::gravity);
        initStateEquil(*grid->c_grid(), *props, deck, eclipseState, grav, state.blackoilState());
        state.faceflux().resize(grid->c_grid()->number_of_faces, 0.0);
    } else {
        initBlackoilStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], state.blackoilState());
    }

    // The capillary pressure is scaled in new_props to match the scaled capillary pressure in props.
    if (deck->hasKeyword("SWATINIT")) {
        const int nc = grid->c_grid()->number_of_cells;
        std::vector<int> cells(nc);
        for (int c = 0; c < nc; ++c) { cells[c] = c; }
        std::vector<double> pc = state.saturation();
        props->capPress(nc, state.saturation().data(), cells.data(), pc.data(),NULL);
        new_props->setSwatInitScaling(state.saturation(),pc);
    }

    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Solver for Newton iterations.
    std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver;
    if (param.getDefault("use_cpr", true)) {
        fis_solver.reset(new NewtonIterationBlackoilCPR(param));
    } else {
        fis_solver.reset(new NewtonIterationBlackoilSimple(param));
    }

    Opm::ScheduleConstPtr schedule = eclipseState->getSchedule();
    Opm::TimeMapConstPtr timeMap(schedule->getTimeMap());
    SimulatorTimer simtimer;

    // initialize variables
    simtimer.init(timeMap);
    if (polymer){
        if (!use_wpolymer) {
            OPM_MESSAGE("Warning: simulate polymer injection without WPOLYMER.");
        } else {
            if (param.has("polymer_start_days")) {
                OPM_MESSAGE("Warning: Using WPOLYMER to control injection since it was found in deck."
                            "You seem to be trying to control it via parameter poly_start_days (etc.) as well.");
            }
        }
    } else {
        if (use_wpolymer) {
            OPM_MESSAGE("Warning: use WPOLYMER in a non-polymer scenario.");
        }
    }    

    bool use_local_perm = param.getDefault("use_local_perm", true);
    Opm::DerivedGeology geology(*grid->c_grid(), *new_props, eclipseState, use_local_perm, grav);

    std::vector<double> threshold_pressures = thresholdPressures(eclipseState, *grid->c_grid());

    SimulatorFullyImplicitBlackoilPolymer<UnstructuredGrid> simulator(param,
                                             *grid->c_grid(),
                                             geology,
                                             *new_props,
                                             polymer_props_ad,
                                             rock_comp->isActive() ? rock_comp.get() : 0,
                                             *fis_solver,
                                             grav,
                                             deck->hasKeyword("DISGAS"),
                                             deck->hasKeyword("VAPOIL"),
                                             polymer,
                                             eclipseState,
                                             outputWriter,
                                             deck,
                                             threshold_pressures);

    if (!schedule->initOnly()){
        std::cout << "\n\n================ Starting main simulation loop ===============\n"
                  << std::flush;

        SimulatorReport fullReport = simulator.run(simtimer, state);

        std::cout << "\n\n================    End of simulation     ===============\n\n";
        fullReport.report(std::cout);

        if (output) {
            std::string filename = output_dir + "/walltime.txt";
            std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
            fullReport.reportParam(tot_os);
            warnIfUnusedParams(param);
        }
    } else {
        outputWriter.writeInit( simtimer );
        std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

