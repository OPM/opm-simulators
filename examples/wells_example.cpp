#include <iostream>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>


#include "opm/core/utility/initState.hpp"
#include "opm/core/utility/SimulatorTimer.hpp"
#include <opm/core/WellsManager.hpp>
#include <opm/core/GridManager.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/newwells.h>
#include <opm/core/grid.h>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/TwophaseState.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

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

    std::vector<int> global_cells(grid.c_grid()->global_cell, grid.c_grid()->global_cell + grid.c_grid()->number_of_cells);

    double gravity[3] = {0.0, 0.0, parameters.getDefault<double>("gravity", 0.0)};
    IncompPropertiesFromDeck incomp_properties(parser, global_cells);

    RockCompressibility rock_comp(parser);

    Opm::LinearSolverFactory linsolver(parameters);

    // EXPERIMENT_ISTL
    IncompTpfa pressure_solver(*grid.c_grid(), incomp_properties.permeability(),
            gravity, linsolver, wells.c_wells());


    std::vector<int> all_cells;
    for (int i = 0; i < grid.c_grid()->number_of_cells; i++) {
        all_cells.push_back(i);
    }

    Opm::TwophaseState state;

    initStateFromDeck(*grid.c_grid(), incomp_properties, parser, gravity[2], state);

    // Compute phase mobilities
    std::vector<double> phase_mob;
    computePhaseMobilities(incomp_properties, all_cells, state.saturation(), phase_mob);
    // Compute total mobility and omega
    std::vector<double> totmob;
    std::vector<double> omega;
    computeTotalMobilityOmega(incomp_properties, all_cells, state.saturation(), totmob, omega);

    std::vector<double> wdp;
    computeWDP(*wells.c_wells(), *grid.c_grid(), state.saturation(), incomp_properties.density(), gravity[2], true, wdp);

    std::vector<double> src;
    Opm::FlowBCManager bcs;

    std::vector<double> pressure;
    std::vector<double> face_flux;

    std::vector<double> well_bhp;
    std::vector<double> well_rate_per_cell;
    std::vector<double> rc;
    rc.resize(grid.c_grid()->number_of_cells);
    
    int nl_pressure_maxiter = 100;
    double nl_pressure_tolerance = 0.0;
    if (rock_comp.isActive()) {
        nl_pressure_maxiter = parameters.getDefault("nl_pressure_maxiter", 10);
        nl_pressure_tolerance = parameters.getDefault("nl_pressure_tolerance", 1.0); // in Pascal
    }

    const int num_cells = grid.c_grid()->number_of_cells;
    std::vector<double> porevol;
    if (rock_comp.isActive()) {
        computePorevolume(*grid.c_grid(), incomp_properties.porosity(), rock_comp, state.pressure(), porevol);
    } else {
        computePorevolume(*grid.c_grid(), incomp_properties.porosity(), porevol);
    }
    if (rock_comp.isActive()) {
        std::vector<double> initial_pressure = state.pressure();
        std::vector<double> prev_pressure;
        for (int iter = 0; iter < nl_pressure_maxiter; ++iter) {
            prev_pressure = state.pressure();
            for (int cell = 0; cell < num_cells; ++cell) {
                rc[cell] = rock_comp.rockComp(state.pressure()[cell]);
            }
            state.pressure() = initial_pressure;
            pressure_solver.solve(totmob, omega, src, wdp, bcs.c_bcs(), porevol, rc, simtimer.currentStepLength(),
                    state.pressure(), state.faceflux(), well_bhp, well_rate_per_cell);
            double max_change = 0.0;
            for (int cell = 0; cell < num_cells; ++cell) {
                max_change = std::max(max_change, std::fabs(state.pressure()[cell] - prev_pressure[cell]));
            }
            std::cout << "Pressure iter " << iter << "   max change = " << max_change << std::endl;
            if (max_change < nl_pressure_tolerance) {
                break;
            }
        }
        computePorevolume(*grid.c_grid(), incomp_properties.porosity(), rock_comp, state.pressure(), porevol);
    } else {
        pressure_solver.solve(totmob, omega, src, wdp, bcs.c_bcs(), state.pressure(), state.faceflux(),
                well_bhp, well_rate_per_cell);
    }

    const int np = incomp_properties.numPhases();
    std::vector<double> fractional_flows(grid.c_grid()->number_of_cells*np, 0.0);
    computeFractionalFlow(incomp_properties, all_cells, state.saturation(), fractional_flows);

    // This will be refactored into a separate function once done
    std::vector<double> well_resflows(wells.c_wells()->number_of_wells*np, 0.0);
    computePhaseFlowRatesPerWell(*wells.c_wells(), well_rate_per_cell, fractional_flows, well_resflows);
    // We approximate (for _testing_ that resflows = surfaceflows)
    for (int wc_iter = 0; wc_iter < 10 && !wells.conditionsMet(well_bhp, well_resflows, well_resflows); ++wc_iter) {
        std::cout << "Conditions not met for well, trying again" << std::endl;
        if (rock_comp.isActive()) {
            std::vector<double> initial_pressure = state.pressure();
            std::vector<double> prev_pressure;
            for (int iter = 0; iter < nl_pressure_maxiter; ++iter) {
                prev_pressure = state.pressure();
                for (int cell = 0; cell < num_cells; ++cell) {
                    rc[cell] = rock_comp.rockComp(state.pressure()[cell]);
                }
                state.pressure() = initial_pressure;
                pressure_solver.solve(totmob, omega, src, wdp, bcs.c_bcs(), porevol, rc, simtimer.currentStepLength(),
                        state.pressure(), state.faceflux(), well_bhp, well_rate_per_cell);
                double max_change = 0.0;
                for (int cell = 0; cell < num_cells; ++cell) {
                    max_change = std::max(max_change, std::fabs(state.pressure()[cell] - prev_pressure[cell]));
                }
                std::cout << "Pressure iter " << iter << "   max change = " << max_change << std::endl;
                if (max_change < nl_pressure_tolerance) {
                    break;
                }
            }
            computePorevolume(*grid.c_grid(), incomp_properties.porosity(), rock_comp, state.pressure(), porevol);
        } else {
            pressure_solver.solve(totmob, omega, src, wdp, bcs.c_bcs(), state.pressure(), state.faceflux(),
                    well_bhp, well_rate_per_cell);
        }
        std::cout << "Solved" << std::endl;

        computePhaseFlowRatesPerWell(*wells.c_wells(), well_rate_per_cell, fractional_flows, well_resflows);
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

