/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Statoil ASA.

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
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/createGlobalCellArray.hpp>

#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/polymer/fullyimplicit/SimulatorFullyImplicitCompressiblePolymer.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/polymer/PolymerState.hpp>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/GridHelpers.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>
#include <opm/common/OpmLog/CounterLog.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/simulators/ensureDirectoryExists.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>

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

    std::cout << "\n================    Test program for fully implicit three-phase black-oil flow     ===============\n\n";
    ParameterGroup param(argc, argv, false);
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
    std::unique_ptr<PolymerBlackoilState> state;
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
        ensureDirectoryExists(output_dir);
        // Write simulation parameters.
        param.writeParam(output_dir + "/simulation.param");
    }

    std::string logFile = output_dir + "/LOGFILE.txt";
    Opm::ParseContext parseContext({{ ParseContext::PARSE_RANDOM_SLASH , InputError::IGNORE }});
    Opm::Parser parser;
    {
        std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(logFile , Opm::Log::DefaultMessageTypes);
        std::shared_ptr<Opm::CounterLog> counterLog = std::make_shared<Opm::CounterLog>(Opm::Log::DefaultMessageTypes);

        Opm::OpmLog::addBackend( "STREAM" , streamLog );
        Opm::OpmLog::addBackend( "COUNTER" , counterLog );
    }

    Deck deck;
    std::shared_ptr<EclipseState> eclipseState;
    try {
        deck = parser.parseFile(deck_filename , parseContext);
        Opm::checkDeck(deck, parser);
        eclipseState.reset(new Opm::EclipseState(deck , parseContext));
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Failed to create valid ECLIPSESTATE object. See logfile: " << logFile << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Grid init

    if (eclipseState->get3DProperties().hasDeckDoubleGridProperty("PORV")) {
        const auto& porv = eclipseState->get3DProperties().getDoubleGridProperty("PORV").getData();
        grid.reset(new GridManager(eclipseState->getInputGrid(), porv));
    } else {
        grid.reset(new GridManager(eclipseState->getInputGrid()));
    }
    auto &cGrid = *grid->c_grid();
    const PhaseUsage pu = Opm::phaseUsageFromDeck(deck);

    // Rock and fluid init

    std::vector<int> compressedToCartesianIdx;
    Opm::createGlobalCellArray(*grid->c_grid(), compressedToCartesianIdx);

    typedef BlackoilPropsAdFromDeck::MaterialLawManager MaterialLawManager;
    auto materialLawManager = std::make_shared<MaterialLawManager>();
    materialLawManager->initFromDeck(deck, *eclipseState, compressedToCartesianIdx);

    props.reset(new BlackoilPropertiesFromDeck( deck, *eclipseState, materialLawManager,
                                                Opm::UgGridHelpers::numCells(cGrid),
                                                Opm::UgGridHelpers::globalCell(cGrid),
                                                Opm::UgGridHelpers::cartDims(cGrid),
                                                param));

    state.reset( new PolymerBlackoilState( Opm::UgGridHelpers::numCells(cGrid), Opm::UgGridHelpers::numFaces(cGrid), 2));
    new_props.reset(new BlackoilPropsAdFromDeck(deck, *eclipseState, materialLawManager, cGrid));
    PolymerProperties polymer_props(deck, *eclipseState);
    PolymerPropsAd polymer_props_ad(polymer_props);

    // Rock compressibility.
    rock_comp.reset(new RockCompressibility(*eclipseState));

    // Gravity.
    gravity[2] = deck.hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], *state);
        initBlackoilSurfvol(*grid->c_grid(), *props, *state);
    } else {
        initStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], *state);
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

    const auto timeMap = eclipseState->getSchedule().getTimeMap();
    SimulatorTimer simtimer;
    simtimer.init(timeMap);


    SimulatorReport rep;
    // With a deck, we may have more epochs etc.
    WellState well_state;
    // Check for WPOLYMER presence in last epoch to decide
    // polymer injection control type.
    const bool use_wpolymer = deck.hasKeyword("WPOLYMER");
    if (use_wpolymer) {
        if (param.has("poly_start_days")) {
            OPM_MESSAGE("Warning: Using WPOLYMER to control injection since it was found in deck. "
                        "You seem to be trying to control it via parameter poly_start_days (etc.) as well.");
        }
    }
    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;

    std::unique_ptr<Opm::EclipseIO>
        eclipseWriter(new Opm::EclipseIO(*eclipseState,
                                         UgGridHelpers
                                         ::createEclipseGrid( cGrid ,
                                                              eclipseState->getInputGrid())));
    Opm::BlackoilOutputWriter
        outputWriter(cGrid, param, *eclipseState, std::move(eclipseWriter), pu);

    SimulatorReport fullReport;
    // Create and run simulator.
    Opm::DerivedGeology geology(*grid->c_grid(), *new_props, *eclipseState, grav);
    SimulatorFullyImplicitCompressiblePolymer<UnstructuredGrid>
        simulator(param,
                  *grid->c_grid(),
                  geology,
                  *new_props,
                  polymer_props_ad,
                  rock_comp->isActive() ? rock_comp.get() : 0,
                  eclipseState,
                  outputWriter,
                  deck,
                  *fis_solver,
                  grav);
    fullReport= simulator.run(simtimer, *state);

    std::cout << "\n\n================    End of simulation     ===============\n\n";
    fullReport.report(std::cout);

    if (output) {
        std::string filename = output_dir + "/walltime.param";
        std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
        fullReport.reportParam(tot_os);
        warnIfUnusedParams(param);
    }

}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

