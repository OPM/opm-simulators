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

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include <filesystem>
#include <boost/version.hpp>

namespace Opm
{

/// Set up a property tree intended for FlexibleSolver by either reading
/// the tree from a JSON file or creating a tree giving the default solver
/// and preconditioner. If the latter, the parameters --linear-solver-reduction,
/// --linear-solver-maxiter and --linear-solver-verbosity are used, but if reading
/// from file the data in the JSON file will override any other options.
PropertyTree
setupPropertyTree(FlowLinearSolverParameters p, // Note: copying the parameters to potentially override.
                  bool linearSolverMaxIterSet,
                  bool linearSolverReductionSet)
{
    std::string conf = p.linsolver_;

    // Get configuration from file.
    if (conf.size() > 5 && conf.substr(conf.size() - 5, 5) == ".json") { // the ends_with() method is not available until C++20
#if BOOST_VERSION / 100 % 1000 > 48
        if ( !std::filesystem::exists(conf) ) {
            OPM_THROW(std::invalid_argument, "JSON file " + conf + " does not exist.");
        }
        try {
            return PropertyTree(conf);
        }
        catch (...) {
            OPM_THROW(std::invalid_argument, "Failed reading linear solver configuration from JSON file " + conf);
        }
#else
        OPM_THROW(std::invalid_argument,
                  "--linear-solver-configuration=file.json not supported with "
                  "boost version. Needs version > 1.48.");
#endif
    }

    // Use CPR configuration.
    if ((conf == "cpr_trueimpes") || (conf == "cpr_quasiimpes") || (conf == "cpr_trueimpesanalytic")) {
        if (!linearSolverMaxIterSet) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_maxiter_ = 20;
        }
        if (!linearSolverReductionSet) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_reduction_ = 0.005;
        }
        return setupCPR(conf, p);
    }

    if ((conf == "cpr") || (conf == "cprw")) {
        if (!linearSolverMaxIterSet) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_maxiter_ = 20;
        }
        if (!linearSolverReductionSet) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_reduction_ = 0.005;
        }
        return setupCPRW(conf, p);
    }

    if (conf == "amg") {
        return setupAMG(conf, p);
    }

    // Use ILU0 configuration.
    if (conf == "ilu0") {
        return setupILU(conf, p);
    }

    if (conf == "dilu") {
        return setupDILU(conf, p);
    }

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
    OPM_THROW(std::invalid_argument,
              conf + " is not a valid setting for --linear-solver-configuration."
              " Please use ilu0, dilu, cpr, cprw, cpr_trueimpes, cpr_quasiimpes, cpr_trueimpesanalytic or isai");
}

std::string getSolverString(const FlowLinearSolverParameters& p)
{
    if (p.newton_use_gmres_)
    {
        return {"gmres"};
    }
    else
    {
        return {"bicgstab"};
    }
}

PropertyTree
setupCPRW(const std::string& /*conf*/, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));
    prm.put("preconditioner.type", "cprw"s);
    prm.put("preconditioner.use_well_weights", "false"s);
    prm.put("preconditioner.add_wells", "true"s);
    prm.put("preconditioner.weight_type", "trueimpes"s);
    prm.put("preconditioner.finesmoother.type", "ParOverILU0"s);
    prm.put("preconditioner.finesmoother.relaxation", 1.0);
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.maxiter", 1);
    prm.put("preconditioner.coarsesolver.tol", 1e-1);
    prm.put("preconditioner.coarsesolver.solver", "loopsolver"s);
    prm.put("preconditioner.coarsesolver.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.type", "amg"s);
    prm.put("preconditioner.coarsesolver.preconditioner.alpha", 0.333333333333);
    prm.put("preconditioner.coarsesolver.preconditioner.relaxation", 1.0);
    prm.put("preconditioner.coarsesolver.preconditioner.iterations", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200);
    prm.put("preconditioner.coarsesolver.preconditioner.pre_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.post_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.beta", 0.0);
    prm.put("preconditioner.coarsesolver.preconditioner.smoother", "ILU0"s);
    prm.put("preconditioner.coarsesolver.preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.maxlevel", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.skip_isolated", 0);
    // We request to accumulate data to 1 process always as our matrix
    // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
    // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
    prm.put("preconditioner.coarsesolver.preconditioner.accumulate", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.prolongationdamping", 1.0);
    prm.put("preconditioner.coarsesolver.preconditioner.maxdistance", 2);
    prm.put("preconditioner.coarsesolver.preconditioner.maxconnectivity", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.maxaggsize", 6);
    prm.put("preconditioner.coarsesolver.preconditioner.minaggsize", 4);
    return prm;
}

PropertyTree
setupCPR(const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));
    prm.put("preconditioner.type", "cpr"s);
    if (conf == "cpr_quasiimpes") {
        prm.put("preconditioner.weight_type", "quasiimpes"s);
    } else if (conf == "cpr_trueimpes") {
        prm.put("preconditioner.weight_type", "trueimpes"s);
    } else {
        prm.put("preconditioner.weight_type", "trueimpesanalytic"s);
    }
    prm.put("preconditioner.finesmoother.type", "ParOverILU0"s);
    prm.put("preconditioner.finesmoother.relaxation", 1.0);
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.maxiter", 1);
    prm.put("preconditioner.coarsesolver.tol", 1e-1);
    prm.put("preconditioner.coarsesolver.solver", "loopsolver"s);
    prm.put("preconditioner.coarsesolver.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.type", "amg"s);
    prm.put("preconditioner.coarsesolver.preconditioner.alpha", 0.333333333333);
    prm.put("preconditioner.coarsesolver.preconditioner.relaxation", 1.0);
    prm.put("preconditioner.coarsesolver.preconditioner.iterations", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200);
    prm.put("preconditioner.coarsesolver.preconditioner.pre_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.post_smooth", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.beta", 0.0);
    prm.put("preconditioner.coarsesolver.preconditioner.smoother", "ILU0"s);
    prm.put("preconditioner.coarsesolver.preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.maxlevel", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.skip_isolated", 0);
    // We request to accumulate data to 1 process always as our matrix
    // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
    // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
    prm.put("preconditioner.coarsesolver.preconditioner.accumulate", 1);
    prm.put("preconditioner.coarsesolver.preconditioner.prolongationdamping", 1.0);
    prm.put("preconditioner.coarsesolver.preconditioner.maxdistance", 2);
    prm.put("preconditioner.coarsesolver.preconditioner.maxconnectivity", 15);
    prm.put("preconditioner.coarsesolver.preconditioner.maxaggsize", 6);
    prm.put("preconditioner.coarsesolver.preconditioner.minaggsize", 4);
    return prm;
}


PropertyTree
setupAMG([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));
    prm.put("preconditioner.type", "amg"s);
    prm.put("preconditioner.alpha", 0.333333333333);
    prm.put("preconditioner.relaxation", 1.0);
    prm.put("preconditioner.iterations", 20);
    prm.put("preconditioner.coarsenTarget", 1200);
    prm.put("preconditioner.pre_smooth", 1);
    prm.put("preconditioner.post_smooth", 1);
    prm.put("preconditioner.beta", 1e-5);
    prm.put("preconditioner.smoother", "ILU0"s);
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


PropertyTree
setupILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));
    prm.put("preconditioner.type", "ParOverILU0"s);
    prm.put("preconditioner.relaxation", p.ilu_relaxation_);
    prm.put("preconditioner.ilulevel", p.ilu_fillin_level_);
    return prm;
}

PropertyTree
setupDILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));
    prm.put("preconditioner.type", "DILU"s);
    return prm;
}


PropertyTree
setupUMFPack([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "umfpack"s);
    return prm;
}


} // namespace Opm
