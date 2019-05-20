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

#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <boost/property_tree/json_parser.hpp>

namespace Opm
{

boost::property_tree::ptree setupPropertyTree(const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    if (p.linear_solver_configuration_json_file_ != "none") {
        boost::property_tree::read_json(p.linear_solver_configuration_json_file_, prm);
    } else {
        prm.put("tol", p.linear_solver_reduction_);
        prm.put("maxiter", p.linear_solver_maxiter_);
        prm.put("verbosity", p.linear_solver_verbosity_);
        prm.put("preconditioner", "ILU0");
        prm.put("solver", "bicgstab");
        prm.put("w", 1.0);
        prm.put("n", 1);
    }

    /*
    prm.put("tol", 0.5);
    prm.put("maxiter", 20);
    prm.put("preconditioner", "cpr");
    prm.put("w", 1);
    prm.put("n", 1);
    prm.put("cpr.finesmoother.preconditioner", "ILU0");
    prm.put("cpr.finesmoother.w", 1);
    prm.put("cpr.finesmoother.n", 1);
    prm.put("cpr.coarsesolver.tol", 0.5);
    prm.put("cpr.coarsesolver.maxiter", 20);
    prm.put("cpr.coarsesolver.preconditioner", "amg");
    prm.put("cpr.coarsesolver.amg.maxlevel", 5);
    prm.put("cpr.coarsesolver.amg.coarsenTarget", 1000);
    prm.put("cpr.coarsesolver.amg.smoother", "ILU0");
    prm.put("cpr.coarsesolver.amg.alpha", 0.2);
    prm.put("cpr.coarsesolver.amg.beta", 0.0001);
    prm.put("cpr.coarsesolver.amg.verbosity", 0);
    prm.put("cpr.coarsesolver.amg.n", 1);
    prm.put("cpr.coarsesolver.amg.w", 1);
    prm.put("cpr.coarsesolver.verbosity", 0);
    prm.put("cpr.coarsesolver.solver", "bicgstab");
    prm.put("cpr.verbosity", 11);
    prm.put("cpr.weights_filename" , "weight_cpr.txt");
    prm.put("cpr.pressure_var_index" , 1);
    prm.put("verbosity", 10);
    prm.put("solver", "bicgstab");
    prm.put("restart", 20);
    */

    return prm;
}

} // namespace Opm
/*
{
    "tol": "0.5",
    "maxiter": "20",
    "preconditioner": "cpr",
    "w": "1",
    "n": "1",
    "amg": "",
    "cpr": {
        "finesmoother": {
            "preconditioner": "ILU0",
            "w": "1",
            "n": "1"
        },
        "coarsesolver": {
            "tol": "0.5",
            "maxiter": "20",
            "preconditioner": "amg",
            "amg": {
                "maxlevel": "5",
                "coarsenTarget": "1000",
                "smoother": "ILU0",
                "alpha": "0.2",
                "beta": "0.0001",
                "verbosity": "0",
                "n": "1",
                "w": "1"
            },
            "verbosity": "0",
            "solver": "bicgstab"
        },
        "verbosity": "11",
        "weights_filename" : "weight_cpr.txt",
        "pressure_var_index" : "1"
    },
    "verbosity": "10",
    "solver": "bicgstab",
    "restart": "20"
}
*/
