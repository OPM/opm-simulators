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
#include <fmt/format.h>

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
        std::ranges::transform(requested, requested.begin(), ::tolower);

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
                  bool linearSolverReductionSet,
                  bool tpsaSetup)
{
    std::string conf = p.linsolver_;

    // Get configuration from file.
    if (conf.size() > 5 && conf.substr(conf.size() - 5, 5) == ".json") { // the ends_with() method is not available until C++20
        if ( !std::filesystem::exists(conf) ) {
            OPM_THROW(std::invalid_argument, "JSON file " + conf + " does not exist.");
        }
        PropertyTree tree;
        try {
            tree = PropertyTree(conf);
        }
        catch (...) {
            OPM_THROW(std::invalid_argument, "Failed reading linear solver configuration from JSON file " + conf);
        }
        validateSystemCPRTree(tree);
        return tree;
    }

    // We use lower case as the internal canonical representation of solver names.
    std::ranges::transform(conf, conf.begin(), ::tolower);

    // Use CPR configuration.
    if (!tpsaSetup) {
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
        
        // System CPR configuration (coupled reservoir-well system solver).
        // system_cprw differs only in that its pressure stage carries the well
        // unknowns, exactly as cprw does relative to cpr.
        if ((conf == "system_cpr") || (conf == "system_cprw")) {
            if (!linearSolverMaxIterSet) {
                p.linear_solver_maxiter_ = 20;
            }
            if (!linearSolverReductionSet) {
                p.linear_solver_reduction_ = 0.005;
            }
            return setupSystemCPR(conf, p);
        }

        // The same algorithm composed from named parts rather than a fixed
        // three-stage sequence; identical results with these settings.
        if ((conf == "general_system_cpr") || (conf == "general_system_cprw")) {
            if (!linearSolverMaxIterSet) {
                p.linear_solver_maxiter_ = 20;
            }
            if (!linearSolverReductionSet) {
                p.linear_solver_reduction_ = 0.005;
            }
            return setupGeneralSystemCPR(conf, p);
        }
    }

    if (conf == "amg") {
        return setupAMG(conf, p);
    }

    // Use ILU0 configuration.
    if (conf == "ilu0") {
        return setupILU(conf, p);
    }

    // mixed-precision ILU0
    if (conf == "mixed-ilu0") {
        return setupMixedILU(conf, p);
    }

    // mixed-precision ILU0
    if (conf == "mixed-dilu") {
        return setupMixedDILU(conf, p);
    }

    // mixed-precision ILU0
    if (conf == "legacy-mixed-ilu0") {
        return setupLegacyMixedILU(conf, p);
    }

    // mixed-precision ILU0
    if (conf == "legacy-mixed-dilu") {
        return setupLegacyMixedDILU(conf, p);
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
    if (tpsaSetup) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("No valid settings found for --tpsa-linear-solver={}! "
                              "Valid preset options are: ilu0, dilu, amg, or umfpack.", conf));
    }
    else {
        OPM_THROW(std::invalid_argument,
                conf + " is not a valid setting for --linear-solver-configuration."
                " Please use ilu0, dilu, isai, cpr, cprw, cpr_trueimpes, cpr_quasiimpes, cpr_trueimpesanalytic,"
                " system_cpr, system_cprw, general_system_cpr or general_system_cprw");
    }


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
    // "auto" keeps the historical coarse well diagonal: contract D for standard
    // wells, minus the row sum for multisegment wells. "contract_d" contracts D
    // for both, which keeps the segment-to-segment coupling.
    prm.put("preconditioner.well_coarse_diagonal", "auto"s);
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
setupMixedILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "mixed-precision"s);
    prm.put("preconditioner.type", "mixed-ilu0"s);
    return prm;
}

PropertyTree
setupMixedDILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "mixed-precision"s);
    prm.put("preconditioner.type", "mixed-dilu"s);
    return prm;
}


PropertyTree
setupLegacyMixedILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "mixed-bicgstab"s);
    prm.put("preconditioner.type", "legacy-mixed-ilu0"s);
    return prm;
}

PropertyTree
setupLegacyMixedDILU([[maybe_unused]] const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    PropertyTree prm;
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", "mixed-bicgstab"s);
    prm.put("preconditioner.type", "legacy-mixed-dilu"s);
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


namespace {

// Join a node key with a leaf name, allowing the node to be the tree root.
std::string at_(const std::string& at, const std::string& leaf)
{
    return at.empty() ? leaf : at + "." + leaf;
}

// The sub-solvers of the system preconditioner, each written at a caller
// chosen key.  Both the fixed three-stage layout and the general one use
// these, so the trees they produce are identical and the two preconditioners
// can be compared setting for setting.
void setupSystemReservoirSmoother(PropertyTree& prm,
                                  const std::string& at,
                                  const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    prm.put(at_(at, "maxiter"), 1);
    prm.put(at_(at, "tol"), p.linear_solver_reduction_);
    prm.put(at_(at, "verbosity"), 0);
    // Apply the configured smoother once without loopsolver's extra defect
    // updates and convergence bookkeeping.
    prm.put(at_(at, "solver"), "preconditioner2inverseoperator"s);
    prm.put(at_(at, "preconditioner.type"), "paroverilu0"s);
    prm.put(at_(at, "preconditioner.relaxation"), 1.0);
}

void setupSystemReservoirSolver(PropertyTree& prm,
                                const std::string& at,
                                const FlowLinearSolverParameters& p,
                                const bool add_wells)
{
    using namespace std::string_literals;
    prm.put(at_(at, "maxiter"), 1);
    prm.put(at_(at, "tol"), p.linear_solver_reduction_);
    prm.put(at_(at, "verbosity"), 0);
    // Apply the configured reservoir preconditioner once without loopsolver's
    // extra defect updates and convergence bookkeeping.
    prm.put(at_(at, "solver"), "preconditioner2inverseoperator"s);
    prm.put(at_(at, "preconditioner.type"), "cpr"s);
    prm.put(at_(at, "preconditioner.relaxation"), 1.0);
    prm.put(at_(at, "preconditioner.use_well_weights"), "false"s);
    // add_wells promotes the pressure stage from reservoir-only CPR to CPRW
    // over the full (reservoir, well) system.
    prm.put(at_(at, "preconditioner.add_wells"), add_wells ? "true"s : "false"s);
    prm.put(at_(at, "preconditioner.weight_type"), "trueimpes"s);
    prm.put(at_(at, "preconditioner.pre_smooth"), 0);
    prm.put(at_(at, "preconditioner.post_smooth"), 0);
    // Set unused finesmoother to jac to avoid spending time setuping an ILU smoother that won't be used.
    prm.put(at_(at, "preconditioner.finesmoother.type"), "jac"s);
    prm.put(at_(at, "preconditioner.finesmoother.relaxation"), 1.0);
    prm.put(at_(at, "preconditioner.verbosity"), 0);
    prm.put(at_(at, "preconditioner.coarsesolver.maxiter"), 1);
    prm.put(at_(at, "preconditioner.coarsesolver.tol"), 1e-1);
    prm.put(at_(at, "preconditioner.coarsesolver.solver"), "loopsolver"s);
    prm.put(at_(at, "preconditioner.coarsesolver.verbosity"), 0);
    prm.put(at_(at, "preconditioner.coarsesolver.preconditioner.type"), "amg"s);
    setupDuneAMG(prm, at_(at, "preconditioner.coarsesolver.preconditioner."));
}

void setupSystemWellSolver(PropertyTree& prm,
                           const std::string& at,
                           const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    prm.put(at_(at, "maxiter"), 1);
    prm.put(at_(at, "tol"), p.linear_solver_reduction_);
    prm.put(at_(at, "verbosity"), 0);
    prm.put(at_(at, "solver"), "umfpack"s);
}

} // anonymous namespace

namespace {

// The knobs the CPRW pressure stage reads, shared by both layouts.
void setupSystemCPRWellOptions(PropertyTree& prm, const std::string& at = "preconditioner.")
{
    using namespace std::string_literals;
    // How the well equations are contracted to the one coarse unknown each
    // well carries in the CPRW pressure system. Only read when add_wells.
    //   cellavg      - average of the reservoir weights over the well's
    //                  perforated cells, on the conservation equations only.
    //                  This is what cprw does (use_well_weights = false).
    //   cellblockavg - the same average taken per block row instead of per
    //                  well. Identical to cellavg for a standard well, and
    //                  worse for multisegment wells (Norne per-connection:
    //                  2604 against 2569).
    //   quasiimpes   - inv(D)^T e_bhp, normalised. Catastrophic for
    //                  multisegment wells; do not make it the default again
    //                  without re-checking them.
    //   unit         - the pressure row as-is; a debugging baseline.
    prm.put(at + "well_weight_type", "cellavg"s);
    // How the well unknowns take part in the pressure-stage transfer:
    //   full            - restrict the well residual, prolong the bhp correction
    //   no_prolongation - restrict, but discard the bhp correction
    //   classic         - neither, i.e. the classic cprw formulation, so that
    //                     the only remaining difference is numerics
    // classic is the default: on full Norne with one segment per connection it
    // measures best (2558, against 2569 for no_prolongation and 2576 for full),
    // and on standard wells the three are within one iteration of each other.
    // The margin is thin. It is also taken with an exact well solve, which
    // nearly annihilates the well residual; restricting it may well pay once
    // the well solve is inexact.
    // Only read when add_wells.
    prm.put(at + "well_transfer", "classic"s);
    // Give a pressure-controlled well a trivial coarse equation, matching
    // StandardWellEquations::extractCPRPressureMatrix. Only read when add_wells.
    prm.put(at + "well_identity_on_pressure_control", "true"s);
    // How a well's coarse diagonal is formed:
    //   auto       - contract D for single-block wells, minus the row sum for
    //                multisegment ones, i.e. what classic cprw does
    //   contract_d - always contract D
    //   row_sum    - always minus the row sum
    // contract_d is the default. On full Norne with one segment per connection
    // it is worth ~0.4% over row_sum (2569 against 2580), and it is what makes
    // the classic cprw path 2646 rather than 2716 on the same case.
    prm.put(at + "well_coarse_diagonal", "contract_d"s);

}

} // anonymous namespace

PropertyTree
setupSystemCPR(const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    const bool add_wells = (conf == "system_cprw");
    PropertyTree prm;

    // Outer solver
    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));

    // Top-level preconditioner: system_cpr
    prm.put("preconditioner.type", "system_cpr"s);
    setupSystemCPRWellOptions(prm);

    setupSystemReservoirSmoother(prm, "preconditioner.reservoir_smoother"s, p);
    setupSystemReservoirSolver(prm, "preconditioner.reservoir_solver"s, p, add_wells);
    setupSystemWellSolver(prm, "preconditioner.well_solver"s, p);
    return prm;
}

// The same algorithm as setupSystemCPR, expressed as the parts the general
// preconditioner composes:
//
//   coarse_solver { reservoir_solver, well_solver }
//   smoother      { reservoir_smoother, well_solver }
//
// applied in that order.  The sub-trees are byte for byte the ones
// setupSystemCPR writes, and with this layout the general preconditioner
// reproduces the fixed three-stage one exactly.  It exists so that the parts
// can be left out or repeated, which the fixed sequence cannot express.
PropertyTree
setupGeneralSystemCPR(const std::string& conf, const FlowLinearSolverParameters& p)
{
    using namespace std::string_literals;
    const bool add_wells = (conf == "general_system_cprw");
    PropertyTree prm;

    prm.put("maxiter", p.linear_solver_maxiter_);
    prm.put("tol", p.linear_solver_reduction_);
    prm.put("verbosity", p.linear_solver_verbosity_);
    prm.put("solver", getSolverString(p));

    prm.put("preconditioner.type", "general_system_cpr"s);
    prm.put("preconditioner.verbosity", 0);
    // Sweeps before and after the coarse correction, as for cpr. 0/1 is the
    // composition the fixed three-stage system_cpr uses.
    prm.put("preconditioner.pre_smooth", 0);
    prm.put("preconditioner.post_smooth", 1);
    // What to do when a well change introduces something the initial build of
    // the structure did not have:
    //   rebuild - build every part again from this tree
    //   refresh - build only the parts the wells size (the well solves and the
    //             coarse system) and refresh the rest in place, which is what
    //             the fixed three-stage system_cpr does
    // These are different preconditioners, not two routes to the same one: a
    // CPR reservoir solve keeps its hierarchy when refreshed and re-aggregates
    // it when rebuilt.
    prm.put("preconditioner.well_structure_update", "rebuild"s);
    // The reservoir weighting, at the top of the preconditioner as it is for
    // cpr: every part that contracts the reservoir equations uses it.
    prm.put("preconditioner.weight_type", "trueimpes"s);

    if (add_wells) {
        // The coarse space, and the solver for it.  As for cpr, this node is
        // itself the solver spec; the keys beside it describe the space.
        prm.put("preconditioner.coarsesolver.type", "cprw_pressure"s);
        setupSystemCPRWellOptions(prm, "preconditioner.coarsesolver."s);
        prm.put("preconditioner.coarsesolver.maxiter", 1);
        prm.put("preconditioner.coarsesolver.tol", 1e-1);
        prm.put("preconditioner.coarsesolver.solver", "loopsolver"s);
        prm.put("preconditioner.coarsesolver.verbosity", 0);
        prm.put("preconditioner.coarsesolver.preconditioner.type", "amg"s);
        setupDuneAMG(prm, "preconditioner.coarsesolver.preconditioner."s);
    }

    // The sweep: a reservoir smoother followed by a well solve, which together
    // with the coarse correction is what system_cpr does.  Unlike cpr's single
    // finesmoother this is a list, because a coupled system has more than one
    // block to visit and the order matters.
    PropertyTree res;
    setupSystemReservoirSmoother(res, ""s, p);
    res.put("block", "reservoir"s);
    PropertyTree well;
    setupSystemWellSolver(well, ""s, p);
    well.put("block", "well"s);

    if (add_wells) {
        // coarse -> well -> reservoir smoother -> well, the sequence the fixed
        // three-stage system_cprw applies. The leading well solve cleans up
        // after the coarse correction.
        prm.put_child_list("preconditioner.finesmoother.steps", {well, res, well});
    } else {
        // No coarse space: the reservoir CPR solve is an ordinary block step
        // and leads the sweep instead.
        PropertyTree resSolver;
        setupSystemReservoirSolver(resSolver, ""s, p, /*add_wells=*/false);
        resSolver.put("block", "reservoir"s);
        prm.put_child_list("preconditioner.finesmoother.steps",
                           {resSolver, well, res, well});
    }
    return prm;
}

void validateGeneralSystemCPRTree(const PropertyTree& prm)
{
    // Only the composition is checked here: which parts exist is the point of
    // this layout, so nearly everything is optional. The parts themselves are
    // validated by whatever builds them.
    const auto coarse = prm.get_child_optional("preconditioner.coarsesolver");
    const auto smoother = prm.get_child_optional("preconditioner.finesmoother");
    if (!coarse.has_value() && !smoother.has_value()) {
        OPM_THROW(std::invalid_argument,
                  "general_system_cpr JSON configuration needs at least one of the "
                  "'preconditioner.coarsesolver' and 'preconditioner.finesmoother' sub-trees.");
    }
    if (coarse && !coarse->get_child_optional("preconditioner").has_value()) {
        OPM_THROW(std::invalid_argument,
                  "general_system_cpr's 'preconditioner.coarsesolver' is the solver for the "
                  "CPRW pressure system and needs its own 'preconditioner' sub-tree.");
    }
    if (smoother) {
        const auto steps = smoother->get_child_list("steps");
        if (!steps.has_value() || steps->empty()) {
            OPM_THROW(std::invalid_argument,
                      "general_system_cpr's 'preconditioner.finesmoother' needs a non-empty "
                      "'steps' array, each entry naming a \"block\" of 'reservoir' or 'well'.");
        }
        for (const auto& step : *steps) {
            const auto block = step.get("block", std::string{});
            if (block != "reservoir" && block != "well") {
                OPM_THROW(std::invalid_argument,
                          "Every 'preconditioner.finesmoother.steps' entry needs \"block\" set "
                          "to 'reservoir' or 'well'; got '" + block + "'.");
            }
        }
    }
}

void validateSystemCPRTree(const PropertyTree& prm)
{
    const auto type = prm.get("preconditioner.type", std::string{});
    if (type == "general_system_cpr") {
        validateGeneralSystemCPRTree(prm);
        return;
    }
    if (type != "system_cpr") {
        return;
    }
    for (const char* sub : {"reservoir_solver", "reservoir_smoother", "well_solver"}) {
        if (!prm.get_child_optional(std::string("preconditioner.") + sub).has_value()) {
            OPM_THROW(std::invalid_argument,
                      fmt::format("system_cpr JSON configuration is missing the required "
                                  "'preconditioner.{}' sub-tree.", sub));
        }
    }
    // Ensure reservoir_solver uses CPR preconditioner
    auto reservoir_solver = prm.get_child_optional("preconditioner.reservoir_solver");
    if (reservoir_solver) {
        std::string precond_type = reservoir_solver->get<std::string>("preconditioner.type", "");
        if (precond_type != "cpr") {
            OPM_THROW(std::invalid_argument,
                      "In system_cpr configuration, the reservoir_solver must use the CPR preconditioner "
                      "(preconditioner.reservoir_solver.preconditioner.type = 'cpr').");
        }
        // With add_wells the pressure stage is assembled and solved directly by
        // the system preconditioner, which takes its solver settings from the
        // coarsesolver sub-tree rather than from the reservoir_solver wrapper.
        if (reservoir_solver->get("preconditioner.add_wells", false)
            && !reservoir_solver->get_child_optional("preconditioner.coarsesolver").has_value()) {
            OPM_THROW(std::invalid_argument,
                      "In system_cpr configuration with "
                      "preconditioner.reservoir_solver.preconditioner.add_wells = true, the "
                      "'preconditioner.reservoir_solver.preconditioner.coarsesolver' sub-tree is "
                      "required: it configures the solver for the CPRW pressure system.");
        }
    }

    // A Krylov well solver stops on a tolerance, so it performs a different
    // number of inner iterations for each right-hand side. That makes the whole
    // system preconditioner non-stationary, which Krylov methods with short
    // recurrences (bicgstab, cg) and standard GMRES are not allowed to use: they
    // assume a fixed preconditioning operator. The outer solver has to be a
    // flexible one.
    auto well_solver = prm.get_child_optional("preconditioner.well_solver");
    if (well_solver) {
        const auto inner = well_solver->get<std::string>("solver", "bicgstab");
        const bool inner_is_krylov
            = (inner == "bicgstab") || (inner == "gmres") || (inner == "cg") || (inner == "flexgmres");
        const auto outer = prm.get<std::string>("solver", "bicgstab");
        const bool outer_is_flexible = (outer == "flexgmres");
        if (inner_is_krylov && !outer_is_flexible) {
            OPM_THROW(std::invalid_argument,
                      fmt::format("system_cpr is configured with an approximate (Krylov) well "
                                  "solver, 'preconditioner.well_solver.solver' = '{}', which makes "
                                  "the preconditioner vary between applications. The outer solver "
                                  "must then be flexible: set 'solver' to 'flexgmres' (it is "
                                  "currently '{}'), or use a stationary well solver such as "
                                  "'umfpack' or 'preconditioner2inverseoperator'.",
                                  inner, outer));
        }
    }
}


void checkSystemCPRMatrixAddWell(bool matrixAddWellContributions)
{
    if (matrixAddWellContributions) {
        OPM_THROW(std::invalid_argument,
                  "--matrix-add-well-contributions=true is incompatible with "
                  "--linear-solver=system_cpr because the system CPR implementation assumes that well contributions are not added to the matrix.");
    }
}


} // namespace Opm
