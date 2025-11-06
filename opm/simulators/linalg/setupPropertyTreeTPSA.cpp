// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include "config.h"

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTreeTPSA.hpp>
#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>

#include <filesystem>
#include <string>

#include <boost/version.hpp>
#include <fmt/format.h>


namespace Opm {

/*!
* \brief Setup linear solver property tree based on runtime/default parameters
*
* \param p Runtime/default parameters
*
* \note Same as Opm::setupPropertyTree() for Flow, but limit remove Flow specific linear solver variants like cpr_***.
* Moreover, TPSA specific presets should be done here!
*/
PropertyTree setupPropertyTreeTPSA(TpsaLinearSolverParameters p)
{
    // Get linear solver
    std::string conf = p.linsolver_;

    // Get configuration from JSON file
    if (conf.size() > 5 && conf.substr(conf.size() - 5, 5) == ".json") {
#if BOOST_VERSION / 100 % 1000 > 48
        if ( !std::filesystem::exists(conf) ) {
            std::string msg = fmt::format("TPSA: JSON file {} in --tpsa-linear-solver does not exist!", conf);
            OpmLog::error(msg);
            OPM_THROW(std::invalid_argument, msg);
        }
        try {
            return PropertyTree(conf);
        }
        catch (...) {
            std::string msg = fmt::format("TPSA: Failed reading linear solver configuration from JSON file  {}", conf);
            OpmLog::error(msg);
            OPM_THROW(std::invalid_argument, msg);
        }
#else
        std::string msg = fmt::format("TPSA: --tpsa-linear-solver with JSON file (={}) not supported with boost "
                                      "version <= 1.48!", conf);
        OpmLog::error(msg);
        OPM_THROW(std::invalid_argument, msg);
#endif
    }

    // We use lower case as the internal canonical representation of solver names
    std::transform(conf.begin(), conf.end(), conf.begin(), ::tolower);

    // Standard AMG
    if (conf == "amg") {
        return setupAMG(conf, p);
    }

    // ILU setups
    if (conf == "ilu0") {
        return setupILU(conf, p);
    }
    if (conf == "dilu") {
        return setupDILU(conf, p);
    }

    // UMFPACK direct solver
    if (conf == "umfpack") {
        return setupUMFPack(conf, p);
    }

    // At this point, the only separate ISAI implementation is with the OpenCL code, and
    // it will check this argument to see if it should be using ISAI. The parameter tree
    // will be ignored, so this is just a dummy configuration to avoid the throw below.
    // If we are using CPU dune-istl solvers, this will just make "isai" an alias of "ilu".
    if (conf == "isai") {
        return setupILU(conf, p);
    }

    // No valid configuration option found.
    std::string msg = fmt::format("No valid settings found for --tpsa-linear-solver={}! "
                                  "Valid preset options are: ilu0, dilu, amg, or umfpack.", conf);
    OpmLog::error(msg);
    OPM_THROW(std::invalid_argument, msg);
}

}  // namespace Opm