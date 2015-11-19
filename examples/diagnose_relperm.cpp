/*
  Copyright 2015 Statoil ASA.

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
#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <boost/filesystem.hpp>
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>


namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Warning: unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "-------------------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;
    parameter::ParameterGroup param(argc, argv);
    // Read saturation tables.
    EclipseStateConstPtr eclState; 
    ParserPtr parser(new Opm::Parser);
    Opm::DeckConstPtr deck;
    //ParseMode parseMode;
    Opm::ParseMode parseMode({{ ParseMode::PARSE_RANDOM_SLASH , InputError::IGNORE }, 
                              { ParseMode::PARSE_UNKNOWN_KEYWORD, InputError::IGNORE},
                              { ParseMode::PARSE_RANDOM_TEXT, InputError::IGNORE}
                             });
    std::string deck_filename = param.get<std::string>("deck_filename");
    deck = parser->parseFile(deck_filename, parseMode);
    eclState.reset(new EclipseState(deck, parseMode));

    GridManager gm(deck);
    const UnstructuredGrid& grid = *gm.c_grid();
    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
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
        param.writeParam(output_dir + "/relperm.param");
    }

    // Issue a warning if any parameters were unused.
    warnIfUnusedParams(param);

    Opm::time::StopWatch timer;
    timer.start();
    RelpermDiagnostics diagnostic(eclState);
    diagnostic.diagnosis(eclState, deck);
    timer.stop();
    double tt = timer.secsSinceStart();
    std::cout << "relperm diagnostics: " << tt << " seconds." << std::endl;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
