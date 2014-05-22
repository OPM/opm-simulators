/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/autodiff/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/core/utility/share_obj.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>

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

    std::cout << "\n================    Test program for fully implicit three-phase black-oil flow     ===============\n\n";
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
    std::shared_ptr<BlackoilPropsAdInterface> new_props;
    std::shared_ptr<RockCompressibility> rock_comp;
    BlackoilState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    std::string deck_filename = param.get<std::string>("deck_filename");

    Opm::ParserPtr newParser(new Opm::Parser() );
    Opm::DeckConstPtr deck = newParser->parseFile( deck_filename );
    std::shared_ptr<EclipseState> eclipseState(new EclipseState(deck));

    // Grid init
    grid.reset(new GridManager(eclipseState->getEclipseGrid()));
    Opm::EclipseWriter outputWriter(param, deck, share_obj(*grid->c_grid()));

    // Rock and fluid init
    props.reset(new BlackoilPropertiesFromDeck(deck, *grid->c_grid(), param));
    new_props.reset(new BlackoilPropsAdFromDeck(deck, *grid->c_grid()));

    // check_well_controls = param.getDefault("check_well_controls", false);
    // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
    // Rock compressibility.
    rock_comp.reset(new RockCompressibility(deck));

    // Gravity.
    gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
        initBlackoilSurfvol(*grid->c_grid(), *props, state);
        enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
        const PhaseUsage pu = props->phaseUsage();
        if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
            const int np = props->numPhases();
            const int nc = grid->c_grid()->number_of_cells;
            for (int c = 0; c < nc; ++c) {
                state.gasoilratio()[c] = state.surfacevol()[c*np + pu.phase_pos[Gas]]
                    / state.surfacevol()[c*np + pu.phase_pos[Oil]];
            }
        }
    } else {
        initBlackoilStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], state);
    }

    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Solver for Newton iterations.
    std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver;
    if (param.getDefault("use_cpr", false)) {
        fis_solver.reset(new NewtonIterationBlackoilCPR(param));
    } else {
        fis_solver.reset(new NewtonIterationBlackoilSimple(param));
    }

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::ofstream outStream;
    std::string output_dir;
    if (output) {
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        std::string filename = output_dir + "/timing.param";
        outStream.open(filename.c_str(), std::fstream::trunc | std::fstream::out);
        // open file to clean it. The file is appended to in SimulatorTwophase
        filename = output_dir + "/step_timing.param";
        std::fstream step_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
        step_os.close();
        param.writeParam(output_dir + "/simulation.param");
    }

    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << std::flush;

    WellStateFullyImplicitBlackoil well_state;
    Opm::TimeMapPtr timeMap(new Opm::TimeMap(deck));
    SimulatorTimer simtimer;

    // initialize variables
    simtimer.init(timeMap);

    SimulatorReport fullReport;
    for (size_t reportStepIdx = 0; reportStepIdx < timeMap->numTimesteps(); ++reportStepIdx) {
        // Report on start of a report step.
        std::cout << "\n"
                  << "---------------------------------------------------------------\n"
                  << "--------------    Starting report step " << reportStepIdx << "    --------------\n"
                  << "---------------------------------------------------------------\n"
                  << "\n";

        WellsManager wells(eclipseState,
                           reportStepIdx,
                           *grid->c_grid(),
                           props->permeability());

        if (reportStepIdx == 0) {
            // @@@ HACK: we should really make a new well state and
            // properly transfer old well state to it every epoch,
            // since number of wells may change etc.
            well_state.init(wells.c_wells(), state);
        }

        simtimer.setCurrentStepNum(reportStepIdx);

        if (reportStepIdx == 0) {
            outputWriter.writeInit(simtimer);
            outputWriter.writeTimeStep(simtimer, state, well_state.basicWellState());
        }

        // Create and run simulator.
        SimulatorFullyImplicitBlackoil<UnstructuredGrid> simulator(param,
                                                 *grid->c_grid(),
                                                 *new_props,
                                                 rock_comp->isActive() ? rock_comp.get() : 0,
                                                 wells,
                                                 *fis_solver,
                                                 grav);
        SimulatorReport episodeReport = simulator.run(simtimer, state, well_state);

        ++simtimer;

        outputWriter.writeTimeStep(simtimer, state, well_state.basicWellState());
        fullReport += episodeReport;
    }

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

