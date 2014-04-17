/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initStateEquil.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/simulator/BlackoilState.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <boost/filesystem.hpp>

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

    void outputData(const std::string& output_dir,
                    const std::string& name,
                    const std::vector<double>& data)
    {
        std::ostringstream fname;
        fname << output_dir << "/" << name;
        boost::filesystem::path fpath = fname.str();
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        fname << "/" << "initial.txt";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname.str());
        }
        std::copy(data.begin(), data.end(), std::ostream_iterator<double>(file, "\n"));
    }



} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    // Setup.
    parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;
    const std::string deck_filename = param.get<std::string>("deck_filename");
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile(deck_filename);
    const double grav = param.getDefault("gravity", unit::gravity);
    GridManager gm(deck);
    const UnstructuredGrid& grid = *gm.c_grid();
    BlackoilPropertiesFromDeck props(deck, grid, param);
    warnIfUnusedParams(param);

    // Initialisation.
    BlackoilState state;
    initStateEquil(grid, props, deck, grav, state);

    // Output.
    const std::string output_dir = param.getDefault<std::string>("output_dir", "output");
    outputData(output_dir, "pressure", state.pressure());
    outputData(output_dir, "saturation", state.saturation());
    outputData(output_dir, "rs", state.gasoilratio());
    outputData(output_dir, "rv", state.rv());
}
catch (const std::exception& e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
