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
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include <filesystem>
#include <boost/version.hpp>

namespace Opm
{

namespace
{
    // Populate DUNE AMG parameters under the given root prefix.
    inline void setupDuneAMG(PropertyTree& prm, const std::string& root)
    {
        using namespace std::string_literals;
        prm.put(root + "alpha", prm.get<double>(root + "alpha", 0.333333333333));
        prm.put(root + "relaxation", prm.get<double>(root + "relaxation", 1.0));
        prm.put(root + "iterations", prm.get<int>(root + "iterations", 1));
        prm.put(root + "coarsenTarget", prm.get<int>(root + "coarsenTarget", 1200));
        prm.put(root + "pre_smooth", prm.get<int>(root + "pre_smooth", 1));
        prm.put(root + "post_smooth", prm.get<int>(root + "post_smooth", 1));
        prm.put(root + "beta", prm.get<double>(root + "beta", 0.0));
        prm.put(root + "smoother", prm.get<std::string>(root + "smoother", "ilu0"s));
        prm.put(root + "verbosity", prm.get<int>(root + "verbosity", 0));
        prm.put(root + "maxlevel", prm.get<int>(root + "maxlevel", 15));
        prm.put(root + "skip_isolated", prm.get<int>(root + "skip_isolated", 0));
        // We request to accumulate data to 1 process always as our matrix
        // graph might be unsymmetric and hence not supported by the PTScotch/ParMetis
        // calls in DUNE. Accumulating to 1 skips PTScotch/ParMetis
        prm.put(root + "accumulate", prm.get<int>(root + "accumulate", 1));
        prm.put(root + "prolongationdamping", prm.get<double>(root + "prolongationdamping", 1.0));
        prm.put(root + "maxdistance", prm.get<int>(root + "maxdistance", 2));
        prm.put(root + "maxconnectivity", prm.get<int>(root + "maxconnectivity", 15));
        prm.put(root + "maxaggsize", prm.get<int>(root + "maxaggsize", 6));
        prm.put(root + "minaggsize", prm.get<int>(root + "minaggsize", 4));
    }

    // Populate Hypre BoomerAMG defaults under the given root prefix.
    // These defaults mirror HyprePreconditioner.hpp
#if HAVE_HYPRE && (HYPRE_USING_CUDA || HYPRE_USING_HIP)
    inline void setupHypreAMG(PropertyTree& prm, const std::string& root, bool useGpu)
    {
        prm.put(root + "print_level", prm.get<int>(root + "print_level", 0));
        prm.put(root + "max_iter", prm.get<int>(root + "max_iter", 1));
        prm.put(root + "strong_threshold", prm.get<double>(root + "strong_threshold", 0.5));
        prm.put(root + "agg_trunc_factor", prm.get<double>(root + "agg_trunc_factor", 0.3));
        prm.put(root + "interp_type", prm.get<int>(root + "interp_type", 6));
        prm.put(root + "max_levels", prm.get<int>(root + "max_levels", 15));
        prm.put(root + "tolerance", prm.get<double>(root + "tolerance", 0.0));

        prm.put(root + "use_gpu", useGpu);
        if (useGpu) {
            prm.put(root + "relax_type", prm.get<int>(root + "relax_type", 16));
            prm.put(root + "coarsen_type", prm.get<int>(root + "coarsen_type", 8));
            prm.put(root + "agg_num_levels", prm.get<int>(root + "agg_num_levels", 0));
            prm.put(root + "agg_interp_type", prm.get<int>(root + "agg_interp_type", 6));
        } else {
            prm.put(root + "relax_type", prm.get<int>(root + "relax_type", 13));
            prm.put(root + "coarsen_type", prm.get<int>(root + "coarsen_type", 10));
            prm.put(root + "agg_num_levels", prm.get<int>(root + "agg_num_levels", 1));
            prm.put(root + "agg_interp_type", prm.get<int>(root + "agg_interp_type", 4));
        }
    }
#endif

#if HAVE_AMGX
    // Populate AMGX AMG defaults under the given root prefix.
    // These defaults mirror AmgxPreconditioner.hpp (AmgxConfig)
    inline void setupAMGXAMG(PropertyTree& prm, const std::string& root)
    {
        using namespace std::string_literals;
        prm.put(root + "determinism_flag", prm.get<int>(root + "determinism_flag", 0));
        prm.put(root + "print_grid_stats", prm.get<int>(root + "print_grid_stats", 0));
        prm.put(root + "print_solve_stats", prm.get<int>(root + "print_solve_stats", 0));
        prm.put(root + "solver", prm.get<std::string>(root + "solver", "AMG"s));
        prm.put(root + "algorithm", prm.get<std::string>(root + "algorithm", "CLASSICAL"s));
        prm.put(root + "interpolator", prm.get<std::string>(root + "interpolator", "D2"s));
        prm.put(root + "selector", prm.get<std::string>(root + "selector", "PMIS"s));
        prm.put(root + "smoother", prm.get<std::string>(root + "smoother", "BLOCK_JACOBI"s));
        prm.put(root + "presweeps", prm.get<int>(root + "presweeps", 3));
        prm.put(root + "postsweeps", prm.get<int>(root + "postsweeps", 3));
        prm.put(root + "strength_threshold", prm.get<double>(root + "strength_threshold", 0.5));
        prm.put(root + "max_iters", prm.get<int>(root + "max_iters", 1));
        prm.put(root + "setup_frequency",
                prm.get<int>(root + "setup_frequency",
                             Opm::Parameters::Get<Opm::Parameters::CprReuseInterval>()));
    }
#endif

    // Decide and configure the GPU AMG backend. Throws if none available.
    inline void setupGpuAmgBackend(PropertyTree& prm,
                                   const std::string& typeKey,
                                   [[maybe_unused]] const std::string& root)
    {
        using namespace std::string_literals;

        // Respect explicit request if user already set a backend
        std::string requested = prm.get(typeKey, "amg"s);
        std::transform(requested.begin(), requested.end(), requested.begin(), ::tolower);

        if (requested == "amgx") {
#if HAVE_AMGX
            OpmLog::info("\nUsing AMGX for GPU AMG backend\n");
            setupAMGXAMG(prm, root);
            return;
#else
            OPM_THROW(std::invalid_argument, "Requested AMGX, but AMGX support is not available in this build.");
#endif
        }

        if (requested == "hypre") {
#if HAVE_HYPRE
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
            OpmLog::info("\nUsing Hypre for GPU AMG backend\n");
            setupHypreAMG(prm, root, true);
            return;
#else
            OPM_THROW(std::invalid_argument, "Requested Hypre on GPU, but Hypre was built without CUDA/HIP support.");
#endif
#else
            OPM_THROW(std::invalid_argument, "Requested Hypre, but Hypre support is not available in this build.");
#endif
        }

#if HAVE_AMGX
        OpmLog::info("\nAuto-selecting AMGX for GPU AMG backend\n");
        prm.put(typeKey, "amgx"s);
        setupAMGXAMG(prm, root);
        return;
#elif HAVE_HYPRE
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        OpmLog::info("\nAuto-selecting Hypre for GPU AMG backend\n");
        prm.put(typeKey, "hypre"s);
        setupHypreAMG(prm, root, true);
        return;
#else
        OPM_THROW(std::invalid_argument,
                  "GPU accelerator selected, but Hypre is built without GPU support. "
                  "Enable AMGX or build Hypre with CUDA/HIP, or use a different preconditioner (e.g., DILU by setting "
                  "--linear-solver=dilu).");
#endif
#else
        OPM_THROW(std::invalid_argument,
                  "GPU accelerator selected, but no GPU AMG backend is available. "
                  "Enable AMGX or build Hypre with CUDA/HIP, or use a different preconditioner (e.g., DILU by setting "
                  "--linear-solver=dilu).");
#endif
    }
} // anonymous namespace

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

    // We use lower case as the internal canonical representation of solver names.
    std::transform(conf.begin(), conf.end(), conf.begin(), ::tolower);

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
    prm.put("preconditioner.pre_smooth", 0);
    prm.put("preconditioner.post_smooth", 1);
    prm.put("preconditioner.finesmoother.type", "paroverilu0"s);
    prm.put("preconditioner.finesmoother.relaxation", 1.0);
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.maxiter", 1);
    prm.put("preconditioner.coarsesolver.tol", 1e-1);
    prm.put("preconditioner.coarsesolver.solver", "loopsolver"s);
    prm.put("preconditioner.coarsesolver.verbosity", 0);
    prm.put("preconditioner.coarsesolver.preconditioner.type", "amg"s);
    setupDuneAMG(prm, "preconditioner.coarsesolver.preconditioner.");
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
    prm.put("preconditioner.pre_smooth", 0);
    prm.put("preconditioner.post_smooth", 1);
    // Choose finesmoother based on accelerator backend
    if (p.linear_solver_accelerator_ == Parameters::LinearSolverAcceleratorType::GPU) {
        // TODO: Set this to opmilu0 to match CPU setup once ILU0 performance matches DILU
        prm.put("preconditioner.finesmoother.type", "dilu"s);
    } else {
        prm.put("preconditioner.finesmoother.type", "paroverilu0"s);
    }
    prm.put("preconditioner.finesmoother.relaxation", 1.0);
    prm.put("preconditioner.verbosity", 0);
    prm.put("preconditioner.coarsesolver.maxiter", 1);
    prm.put("preconditioner.coarsesolver.tol", 1e-1);
    prm.put("preconditioner.coarsesolver.solver", "loopsolver"s);
    prm.put("preconditioner.coarsesolver.verbosity", 0);
    // Choose coarsesolver AMG backend based on accelerator backend and available AMG backends
    if (p.linear_solver_accelerator_ == Parameters::LinearSolverAcceleratorType::GPU) {
        setupGpuAmgBackend(
            prm, "preconditioner.coarsesolver.preconditioner.type", "preconditioner.coarsesolver.preconditioner.");
        return prm;
    } else {
        prm.put("preconditioner.coarsesolver.preconditioner.type", "amg"s);
        setupDuneAMG(prm, "preconditioner.coarsesolver.preconditioner.");
    }
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

    // Choose AMG backend based on accelerator backend and available AMG backends
    if (p.linear_solver_accelerator_ == Parameters::LinearSolverAcceleratorType::GPU
        || prm.get("preconditioner.type", "amg"s) == "amgx"s
        || prm.get("preconditioner.type", "amg"s) == "hypre"s) {
        setupGpuAmgBackend(prm, "preconditioner.type", "preconditioner.");
    } else {
        prm.put("preconditioner.type", "amg"s);
        setupDuneAMG(prm, "preconditioner.");
        // Override fields that differ for standalone AMG (compared to CPR coarse solver)
        prm.put("preconditioner.beta", 1e-5);
        prm.put("preconditioner.prolongationdamping", 1.6);
    }
    prm.put("preconditioner.iterations", 20);
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
    if (p.linear_solver_accelerator_ == Parameters::LinearSolverAcceleratorType::GPU && !p.is_nldd_local_solver_) {
        // TODO: We could add ParOverILU0 as an alias in the GPU path to simplify this.
        prm.put("preconditioner.type", "opmilu0"s);
    } else {
        prm.put("preconditioner.type", "paroverilu0"s);
    }
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
    prm.put("preconditioner.type", "dilu"s);
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
