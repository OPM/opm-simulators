/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.

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
#include <boost/version.hpp>

namespace Opm
{

boost::property_tree::ptree
setupCPR(const std::string& conf, const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "bicgstab");
    prm.put("preconditioner.type", "cpr");
    if (conf == "cpr_quasiimpes") {
        prm.put("preconditioner.weight_type", "quasiimpes");
    } else {
        prm.put("preconditioner.weight_type", "trueimpes");
    }
    prm.put("preconditioner.finesmoother.type", "ParOverILU0");
    prm.put("preconditioner.finesmoother.relaxation", 1.0);
    prm.put("preconditioner.pressure_var_index", 1);
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.maxiter", 1);
    prm.put("preconditioner.coarsesolver.tol", 1e-1);
    prm.put("preconditioner.coarsesolver.solver", "loopsolver");
    prm.put("preconditioner.coarsesolver.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.type", "amg");
    prm.put("preconditioner.coarsesolver.preconditioner.alpha", 0.333333333333);
    prm.put("preconditioner.coarsesolver.preconditioner.relaxation", 1.0);
    prm.put("preconditioner.coarsesolver.preconditioner.iterations", p.cpr_max_ell_iter_);
    prm.put("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200);
    prm.put("preconditioner.coarsesolver.preconditioner.pre_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.post_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.beta", 1e-5);
    prm.put("preconditioner.coarsesolver.preconditioner.smoother", "ILU0");
    prm.put("preconditioner.coarsesolver.preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.maxlevel", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.skip_isolated", 0);
    // We request to accumulate data to 1 process always as our matrix
    // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
    // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
    prm.put("preconditioner.coarsesolver.preconditioner.accumulate", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.prolongationdamping", 1.6);
    prm.put("preconditioner.coarsesolver.preconditioner.maxdistance", 2);
    prm.put("preconditioner.coarsesolver.preconditioner.maxconnectivity", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.maxaggsize", 6);
    prm.put("preconditioner.coarsesolver.preconditioner.minaggsize", 4);
    return prm;
}


boost::property_tree::ptree
setupAMG([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "bicgstab");
    prm.put("preconditioner.type", "amg");
    prm.put("preconditioner.alpha", 0.333333333333);
    prm.put("preconditioner.relaxation", 1.0);
    prm.put("preconditioner.iterations", 20);
    prm.put("preconditioner.coarsenTarget", 1200);
    prm.put("preconditioner.pre_smooth", 1);
    prm.put("preconditioner.post_smooth", 1);
    prm.put("preconditioner.beta", 1e-5);
    prm.put("preconditioner.smoother", "ILU0");
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.maxlevel", 15);
    prm.put("preconditioner.skip_isolated", 0);
    // We request to accumulate data to 1 process always as our matrix
    // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
    // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
    prm.put("preconditioner.accumulate", 1);
    prm.put("preconditioner.prolongationdamping", 1.6);
    prm.put("preconditioner.maxdistance", 2);
    prm.put("preconditioner.maxconnectivity", 15);
    prm.put("preconditioner.maxaggsize", 6);
    prm.put("preconditioner.minaggsize", 4);
    return prm;
}


boost::property_tree::ptree
setupILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "bicgstab");
    prm.put("preconditioner.type", "ParOverILU0");
    prm.put("preconditioner.relaxation", p.ilu_relaxation_);
    prm.put("preconditioner.ilulevel", p.ilu_fillin_level_);
    return prm;
}



} // namespace Opm
