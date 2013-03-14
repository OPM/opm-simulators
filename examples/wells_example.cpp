#include <iostream>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>


#include <opm/core/utility/initState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/GridManager.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/newwells.h>
#include <opm/core/grid.h>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

int main(int argc, char** argv)
{

    using namespace Opm::parameter;
    using namespace Opm;
    ParameterGroup parameters(argc, argv, false);
    std::string file_name = parameters.getDefault<std::string > ("inputdeck", "data.data");

    SimulatorTimer simtimer;
    simtimer.init(parameters);

    // Read input file
    EclipseGridParser parser(file_name);
    std::cout << "Done!" << std::endl;
    // Setup grid
    GridManager grid(parser);

    // Finally handle the wells
    WellsManager wells(parser, *grid.c_grid(), NULL);

    double gravity[3] = {0.0, 0.0, parameters.getDefault<double>("gravity", 0.0)};
    IncompPropertiesFromDeck incomp_properties(parser, *grid.c_grid());

    RockCompressibility rock_comp(parser);

    Opm::LinearSolverFactory linsolver(parameters);
    double nl_pressure_residual_tolerance = 1e-8;
    double nl_pressure_change_tolerance = 0.0;
    int nl_pressure_maxiter = 100;
    if (rock_comp.isActive()) {
        nl_pressure_residual_tolerance = parameters.getDefault("nl_pressure_residual_tolerance", 1e-8);
        nl_pressure_change_tolerance = parameters.getDefault("nl_pressure_change_tolerance", 1.0); // in Pascal
        nl_pressure_maxiter = parameters.getDefault("nl_pressure_maxiter", 10);
    }

    std::vector<double> src;
    Opm::FlowBCManager bcs;

    // EXPERIMENT_ISTL
    IncompTpfa pressure_solver(*grid.c_grid(), incomp_properties, &rock_comp, linsolver,
                               nl_pressure_residual_tolerance, nl_pressure_change_tolerance, nl_pressure_maxiter,
                               gravity, wells.c_wells(), src, bcs.c_bcs());


    std::vector<int> all_cells;
    for (int i = 0; i < grid.c_grid()->number_of_cells; i++) {
        all_cells.push_back(i);
    }

    Opm::TwophaseState state;

    initStateFromDeck(*grid.c_grid(), incomp_properties, parser, gravity[2], state);

    Opm::WellState well_state;
    well_state.init(wells.c_wells(), state);

    pressure_solver.solve(simtimer.currentStepLength(), state, well_state);

    const int np = incomp_properties.numPhases();
    std::vector<double> fractional_flows(grid.c_grid()->number_of_cells*np, 0.0);
    computeFractionalFlow(incomp_properties, all_cells, state.saturation(), fractional_flows);

    // This will be refactored into a separate function once done
    std::vector<double> well_resflows(wells.c_wells()->number_of_wells*np, 0.0);
    computePhaseFlowRatesPerWell(*wells.c_wells(), well_state.perfRates(), fractional_flows, well_resflows);
    // We approximate (for _testing_ that resflows = surfaceflows)
    for (int wc_iter = 0; wc_iter < 10 && !wells.conditionsMet(well_state.bhp(), well_resflows, well_resflows); ++wc_iter) {
        std::cout << "Conditions not met for well, trying again" << std::endl;
        pressure_solver.solve(simtimer.currentStepLength(), state, well_state);
        std::cout << "Solved" << std::endl;

        computePhaseFlowRatesPerWell(*wells.c_wells(), well_state.perfRates(), fractional_flows, well_resflows);
    }

#if 0
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), incomp_properties, porevol);



    TwophaseFluid fluid(incomp_properties);
    TransportModel model(fluid, *grid->c_grid(), porevol, gravity[2], true);

    TransportSolver tsolver(model);

    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[] = {1.0, 0.0};
    double ssink[] = {0.0, 1.0};
    double zdummy[] = {0.0, 0.0};

    {
        int well_cell_index = 0;
        for (int well = 0; well < wells.c_wells()->number_of_wells; ++well) {
            for (int cell = wells.c_wells()->well_connpos[well]; cell < wells.c_wells()->well_connpos[well + 1]; ++cell) {
                if (well_rate_per_cell[well_cell_index] > 0.0) {
                    append_transport_source(well_cell_index, 2, 0,
                            well_rate_per_cell[well_cell_index], ssrc, zdummy, tsrc);
                } else if (well_rate_per_cell[well_cell_index] < 0.0) {
                    append_transport_source(well_cell_index, 2, 0,
                            well_rate_per_cell[well_cell_index], ssink, zdummy, tsrc);
                }
            }
        }
    }

    tsolver.solve(*grid->c_grid(), tsrc, stepsize, ctrl, state, linsolve, rpt);

    Opm::computeInjectedProduced(*props, state.saturation(), src, stepsize, injected, produced);
#endif
    return 0;
}

