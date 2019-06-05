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

#include <boost/property_tree/json_parser.hpp>

namespace Opm
{

/// Set up a property tree intended for FlexibleSolver by either reading
/// the tree from a JSON file or creating a tree giving the default solver
/// and preconditioner. If the latter, the parameters --linear-solver-reduction,
/// --linear-solver-maxiter and --linear-solver-verbosity are used, but if reading
/// from file the data in the JSON file will override any other options.
boost::property_tree::ptree
setupPropertyTree(const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    if (p.linear_solver_configuration_json_file_ != "none") {
        boost::property_tree::read_json(p.linear_solver_configuration_json_file_, prm);
    } else {
        prm.put("tol", p.linear_solver_reduction_);
        prm.put("maxiter", p.linear_solver_maxiter_);
        prm.put("verbosity", p.linear_solver_verbosity_);
        prm.put("solver", "bicgstab");
        prm.put("preconditioner.type", "ParOverILU0");
        prm.put("preconditioner.relaxation", 1.0);
    }
    return prm;
}

} // namespace Opm
