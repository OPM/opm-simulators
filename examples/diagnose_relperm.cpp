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
#include <opm/core/props/satfunc/RelpermDiagnostics_impl.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>


void usage() {
    std::cout << std::endl << 
        "Usage: diagnose_relperm <eclipseFile>" << std::endl;
}


// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;
    if (argc <= 1) {
        usage();
        exit(1);
    } 
    const char* eclipseFilename = argv[1];
    EclipseStateConstPtr eclState; 
    ParserPtr parser(new Opm::Parser);
    Opm::ParseContext parseContext({{ ParseContext::PARSE_RANDOM_SLASH , InputError::IGNORE }, 
                              { ParseContext::PARSE_UNKNOWN_KEYWORD, InputError::IGNORE},
                              { ParseContext::PARSE_RANDOM_TEXT, InputError::IGNORE},
                              { ParseContext::UNSUPPORTED_SCHEDULE_GEO_MODIFIER, InputError::IGNORE},
                              { ParseContext::UNSUPPORTED_COMPORD_TYPE, InputError::IGNORE},
                              { ParseContext::UNSUPPORTED_INITIAL_THPRES, InputError::IGNORE},
                              { ParseContext::INTERNAL_ERROR_UNINITIALIZED_THPRES, InputError::IGNORE}
                             });
    Opm::DeckConstPtr deck(parser->parseFile(eclipseFilename, parseContext));
    eclState.reset(new EclipseState(deck, parseContext));

    GridManager gm(eclState->getInputGrid());
    const UnstructuredGrid& grid = *gm.c_grid();
    using boost::filesystem::path; 
    path fpath(eclipseFilename);
    std::string baseName;
    if (boost::to_upper_copy(path(fpath.extension()).string())== ".DATA") {
        baseName = path(fpath.stem()).string();
    } else {
        baseName = path(fpath.filename()).string();
    }

    std::string logFile = baseName + ".SATFUNCLOG";
    std::shared_ptr<EclipsePRTLog> prtLog = std::make_shared<EclipsePRTLog>(logFile, Log::DefaultMessageTypes);
    OpmLog::addBackend( "ECLIPSEPRTLOG" , prtLog );
    Opm::time::StopWatch timer;
    timer.start();
    RelpermDiagnostics diagnostic;
    diagnostic.diagnosis(eclState, deck, grid);
    timer.stop();
    double tt = timer.secsSinceStart();
    std::cout << "relperm diagnostics: " << tt << " seconds." << std::endl;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
