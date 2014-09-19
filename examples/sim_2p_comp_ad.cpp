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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/SimulatorCompressibleAd.hpp>
#include <opm/autodiff/GeoProps.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>


#include <boost/filesystem.hpp>

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

    std::cout << "\n================    Test program for weakly compressible two-phase flow     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    Opm::DeckConstPtr deck;
    std::unique_ptr<GridManager> grid;
    std::unique_ptr<BlackoilPropertiesInterface> props;
    std::unique_ptr<RockCompressibility> rock_comp;
    EclipseStateConstPtr eclipseState;
    BlackoilState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };

    std::string deck_filename = param.get<std::string>("deck_filename");
    Opm::ParserPtr parser(new Opm::Parser());
    deck = parser->parseFile( deck_filename );
    eclipseState.reset(new EclipseState(deck));

    // Grid init
    grid.reset(new GridManager(deck));
    // Rock and fluid init
    props.reset(new BlackoilPropertiesFromDeck(deck, eclipseState, *grid->c_grid(), param));
    // check_well_controls = param.getDefault("check_well_controls", false);
    // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
    // Rock compressibility.
    rock_comp.reset(new RockCompressibility(deck, eclipseState));
    // Gravity.
    gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;
    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
    } else {
        initStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], state);
    }
    initBlackoilSurfvol(*grid->c_grid(), *props, state);

    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Linear solver.
    LinearSolverFactory linsolver(param);

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::ofstream epoch_os;
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
        std::string filename = output_dir + "/epoch_timing.param";
        epoch_os.open(filename.c_str(), std::fstream::trunc | std::fstream::out);
        // open file to clean it. The file is appended to in SimulatorTwophase
        filename = output_dir + "/step_timing.param";
        std::fstream step_os(filename.c_str(), std::fstream::trunc | std::fstream::out);
        step_os.close();
        param.writeParam(output_dir + "/simulation.param");
    }

    std::cout << "\n\n================    Starting main simulation loop     ===============\n";

    SimulatorReport rep;
    // With a deck, we may have more epochs etc.
    WellState well_state;
    Opm::TimeMapConstPtr timeMap = eclipseState->getSchedule()->getTimeMap();
    Opm::DerivedGeology geology(*grid->c_grid(), *props, eclipseState);
    SimulatorTimer simtimer;

    for (size_t reportStepIdx = 0; reportStepIdx < timeMap->numTimesteps(); ++reportStepIdx) {
        // Report on start of report step.
        std::cout << "\n\n--------------    Starting report step " << reportStepIdx << "    --------------"
                  << "\n                  (number of steps left: "
                  << timeMap->numTimesteps() - reportStepIdx << ")\n\n" << std::flush;

        // Create new wells, well_state
        WellsManager wells(eclipseState, reportStepIdx, *grid->c_grid(), props->permeability());
        // @@@ HACK: we should really make a new well state and
        // properly transfer old well state to it every report step,
        // since number of wells may change etc.
        if (reportStepIdx == 0) {
            well_state.init(wells.c_wells(), state);
        }

        simtimer.setCurrentStepNum(reportStepIdx);

        // Create and run simulator.
        SimulatorCompressibleAd simulator(param,
                                          *grid->c_grid(),
                                          geology,
                                          *props,
                                          rock_comp->isActive() ? rock_comp.get() : 0,
                                          wells,
                                          linsolver,
                                          grav);
        if (reportStepIdx == 0) {
            warnIfUnusedParams(param);
        }
        SimulatorReport epoch_rep = simulator.run(simtimer, state, well_state);
        if (output) {
            epoch_rep.reportParam(epoch_os);
        }
        // Update total timing report and remember step number.
        rep += epoch_rep;
    }

    std::cout << "\n\n================    End of simulation     ===============\n\n";
    rep.report(std::cout);

    if (output) {
        std::string filename = output_dir + "/walltime.param";
        std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
        rep.reportParam(tot_os);
    }

}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

