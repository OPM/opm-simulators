/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <boost/version.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Opm
{

/// Set up a property tree intended for FlexibleSolver by either reading
/// the tree from a JSON file or creating a tree giving the default solver
/// and preconditioner. If the latter, the parameters --linear-solver-reduction,
/// --linear-solver-maxiter and --linear-solver-verbosity are used, but if reading
/// from file the data in the JSON file will override any other options.
template<class TypeTag>
boost::property_tree::ptree
setupPropertyTree(FlowLinearSolverParameters p) // Note: copying the parameters to potentially override.
{
    std::string conf = p.linear_solver_configuration_;

    // Get configuration from file.
    if (conf == "file") {
#if BOOST_VERSION / 100 % 1000 > 48
        if (p.linear_solver_configuration_json_file_ == "none") {
            OPM_THROW(std::invalid_argument,
                      "--linear-solver-configuration=file requires that a filename "
                          << "is passed with "
                          << "--linear-solver-configuration-json-file=filename.");
        } else {
            boost::property_tree::ptree prm;
            boost::property_tree::read_json(p.linear_solver_configuration_json_file_, prm);
            return prm;
        }
#else
        OPM_THROW(std::invalid_argument,
                  "--linear-solver-configuration=file not supported with "
                      << "boost version. Needs version > 1.48.");
#endif
    }

    // Use CPR configuration.
    if ((conf == "cpr") || (conf == "cpr_trueimpes") || (conf == "cpr_quasiimpes")) {
        if (conf == "cpr") {
            // Treat "cpr" as short cut for the true IMPES variant.
            conf = "cpr_trueimpes";
        }
        if (!EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter)) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_maxiter_ = 20;
        }
        if (!EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter)) {
            // Use our own default unless it was explicitly overridden by user.
            p.cpr_max_ell_iter_ = 1;
        }
        return setupCPR(conf, p);
    }

    // Use ILU0 configuration.
    if (conf == "ilu0") {
        return setupILU(conf, p);
    }

    // No valid configuration option found.
    OPM_THROW(std::invalid_argument,
              conf << " is not a valid setting for --linear-solver-configuration."
              << " Please use ilu0, cpr, cpr_trueimpes, or cpr_quasiimpes");
}



} // namespace Opm
