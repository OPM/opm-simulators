/*
*/

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

#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/polymer/SimulatorFullyImplicitTwophase.hpp>
#include <opm/autodiff/polymer/IncompPropsAdInterface.hpp>
#include <opm/autodiff/polymer/IncompPropsAdBasic.hpp>
#include <opm/autodiff/polymer/IncompPropsAdFromDeck.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>

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

    std::cout << "\n================    Test program for incompressible two-phase flow     ===============\n\n";
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;


    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<EclipseGridParser> deck;
    boost::scoped_ptr<GridManager> grid;
    boost::scoped_ptr<IncompPropsAdInterface> props;
    TwophaseState state;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        deck.reset(new EclipseGridParser(deck_filename));
        // Grid init
        grid.reset(new GridManager(*deck));
        // Rock and fluid init
        props.reset(new IncompPropsAdFromDeck(*deck, *grid->c_grid()));
        // Gravity.
        gravity[2] = deck->hasField("NOGRAV") ? 0.0 : unit::gravity;
        // Init state variables (saturation and pressure).
        int num_cells = grid->c_grid()->number_of_cells;
        if (param.has("init_saturation")) {
            //initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
            const double init_saturation = param.get<double>("init_saturation");
            for (int c = 0; c < num_cells; ++c) {
                state.saturation()[2*c] = init_saturation;
                state.saturation()[2*c+1] = 1. - init_saturation;
            }
        } else {
            if (deck->hasField("PRESSURE")) {
                // Set saturations from SWAT/SGAS, pressure from PRESSURE.
                std::vector<double>& s = state.saturation();
                std::vector<double>& p = state.pressure();
                const std::vector<double>& p_deck = deck->getFloatingPointValue("PRESSURE");
                // water-oil or water-gas: we require SWAT
                if (!deck->hasField("SWAT")) {
                    OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SWAT keyword in 2-phase init");
                }
                const std::vector<double>& sw_deck = deck->getFloatingPointValue("SWAT");
                for (int c = 0; c < num_cells; ++c) {
                    int c_deck = (grid->c_grid()->global_cell == NULL) ? c : grid->c_grid()->global_cell[c];
                    s[2*c] = sw_deck[c_deck];
                    s[2*c + 1] = 1.0 - sw_deck[c_deck];
                    p[c] = p_deck[c_deck];
                }
            }
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
        // Rock and fluid init.
        props.reset(new IncompPropsAdBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Rock compressibility.
        // Gravity.
        gravity[2] = param.getDefault("gravity", 0.0);
        int num_cells = grid->c_grid()->number_of_cells;
    }

    // Warn if gravity but no density difference.
    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Initialising src
    std::vector<double> src(num_cells, 0.0);
    if (use_deck) {
        // Do nothing, wells will be the driving force, not source terms.
        if (deck->hasField("SRC")) {
            const std::vector<double>& src_deck = deck->getFloatingPointValue("SRC");
            for (int c = 0; c < num_cells; ++c) {
                int c_deck = (grid->c_grid()->global_cell == NULL) ? c : grid->c_grid()->global_cell[c];
                src[c] = src_deck[c_deck];
            }
        }
    } else {
        // Compute pore volumes, in order to enable specifying injection rate
        // terms of total pore volume.
        std::vector<double> porevol;
        computePorevolume(*grid->c_grid(), props->porosity(), porevol);
        const double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        const double default_injection = use_gravity ? 0.0 : 0.1;
        const double flow_per_sec = param.getDefault<double>("injected_porevolumes_per_day", default_injection)
            *tot_porevol_init/unit::day;
        src[0] = flow_per_sec;
        src[num_cells - 1] = -flow_per_sec;
    }

    in: num_cells = grid->c_grid()->number_of_cells;

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


    std::cout << "\n\n================    Starting main simulation loop     ===============\n"
              << "                        (number of epochs: "
              << (use_deck ? deck->numberOfEpochs() : 1) << ")\n\n" << std::flush;

    SimulatorReport rep;
    if (!use_deck) {
        // Simple simulation without a deck.
        SimulatorFullyImplicitTwophase simulator(param,
                                            *grid->c_grid(),
                                            *props,
                                            linsolver,
                                            src);
        SimulatorTimer simtimer;
        simtimer.init(param);
        warnIfUnusedParams(param);
        rep = simulator.run(simtimer, state, src);
    } else {
        // With a deck, we may have more epochs etc.
        int step = 0;
        SimulatorTimer simtimer;
        // Use timer for last epoch to obtain total time.
        deck->setCurrentEpoch(deck->numberOfEpochs() - 1);
        simtimer.init(*deck);
        const double total_time = simtimer.totalTime();
        for (int epoch = 0; epoch < deck->numberOfEpochs(); ++epoch) {
            // Set epoch index.
            deck->setCurrentEpoch(epoch);

            // Update the timer.
            if (deck->hasField("TSTEP")) {
                simtimer.init(*deck);
            } else {
                if (epoch != 0) {
                    OPM_THROW(std::runtime_error, "No TSTEP in deck for epoch " << epoch);
                }
                simtimer.init(param);
            }
            simtimer.setCurrentStepNum(step);
            simtimer.setTotalTime(total_time);

            // Report on start of epoch.
            std::cout << "\n\n--------------    Starting epoch " << epoch << "    --------------"
                      << "\n                  (number of steps: "
                      << simtimer.numSteps() - step << ")\n\n" << std::flush;


            // Create and run simulator.
            SimulatorFullyImplicitTwophase simulator(param,
                                                *grid->c_grid(),
                                                *props,
                                                linsolver,
                                                src);
            if (epoch == 0) {
                warnIfUnusedParams(param);
            }
            SimulatorReport epoch_rep = simulator.run(simtimer, state, src);
            if (output) {
                epoch_rep.reportParam(epoch_os);
            }
            // Update total timing report and remember step number.
            rep += epoch_rep;
            step = simtimer.currentStepNum();
        }
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

