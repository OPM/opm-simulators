#include <iostream>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/WellsManager.hpp>
#include <opm/core/GridManager.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include "opm/core/newwells.h"
#include "opm/core/grid.h"

int main(int argc, char** argv) {

    using namespace Opm::parameter;
    using namespace Opm;
    ParameterGroup parameters( argc, argv, false );
    std::string file_name = parameters.getDefault<std::string>("inputdeck", "data.data");

    // Read input file
    EclipseGridParser parser(file_name);
    std::cout << "Done!" << std::endl;
    // Setup grid
    GridManager grid(parser);
    
    // Finally handle the wells
    WellsManager wells(parser, *grid.c_grid(), NULL);
    
    std::vector<int> global_cells(grid.c_grid()->global_cell, grid.c_grid()->global_cell + grid.c_grid()->number_of_cells);

    std::cout << "ahoi" << std::endl;
    double gravity[3] = {0.0, 0.0, parameters.getDefault<double>("gravity", 0.0)};
    IncompPropertiesFromDeck incomp_properties(parser, global_cells);
        std::cout << "there" << std::endl;

    LinearSolverUmfpack umfpack_solver;
    std::cout << "here" << std::endl;
    IncompTpfa pressure_solver(*grid.c_grid(), incomp_properties.permeability(), 
                               gravity, umfpack_solver,  wells.c_wells());
    
    return 0;
}

