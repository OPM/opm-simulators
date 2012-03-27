#include <iostream>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/WellsManager.hpp>
#include <opm/core/GridManager.hpp>

#include "opm/core/newwells.h"

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
    
    return 0;
}

