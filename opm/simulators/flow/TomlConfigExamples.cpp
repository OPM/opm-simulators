/*
  Copyright 2026 Equinor ASA.

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

#include <opm/simulators/flow/TomlConfigExamples.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/SetupPartitioningParams.hpp>
#endif

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

constexpr std::size_t commentWidth = 78;

/// Word-wrap text to fit in TOML comments.  Explicit line breaks are kept,
/// continuation lines of list items ("- ...") are indented, and indented
/// lines, e.g. command lines, are kept as they are.
std::string wrap(std::string_view text, const std::size_t width = commentWidth - 2)
{
    std::string out;
    std::size_t lineStart = 0;
    while (lineStart <= text.size()) {
        auto lineEnd = text.find('\n', lineStart);
        if (lineEnd == std::string_view::npos) {
            lineEnd = text.size();
        }
        const auto line = text.substr(lineStart, lineEnd - lineStart);
        lineStart = lineEnd + 1;

        if (line.starts_with(' ')) {
            out += std::string { line } + '\n';
            continue;
        }

        const auto indent = std::string(line.starts_with("- ") ? 2 : 0, ' ');
        std::string current;
        std::istringstream words { std::string{line} };
        for (std::string word; words >> word; ) {
            if (!current.empty() && (current.size() + 1 + word.size() > width)) {
                out += current + '\n';
                current = indent;
            }
            current += (current.empty() || current == indent) ? word : ' ' + word;
        }
        out += current + '\n';
    }

    // Remove trailing new line.
    if (!out.empty()) {
        out.pop_back();
    }
    return out;
}

std::string commentBlock(std::string_view text)
{
    std::string out;
    std::istringstream lines { wrap(text) };
    for (std::string line; std::getline(lines, line); ) {
        out += line.empty() ? "#\n" : "# " + line + '\n';
    }
    return out;
}

std::string toLower(std::string s)
{
    std::ranges::transform(s, s.begin(),
                           [](unsigned char c) { return std::tolower(c); });
    return s;
}

struct KeyDoc
{
    std::string_view description;
    std::string_view allowed{};
};

std::string formatDoc(const KeyDoc& doc)
{
    auto out = std::string { doc.description };
    if (! doc.allowed.empty()) {
        out += "\nAllowed values: ";
        out += doc.allowed;
    }
    return wrap(out);
}

// ===========================================================================
// Linear solver
// ===========================================================================

constexpr std::string_view solverTypeDoc =
    "Linear solver. Allowed values:\n"
    "- bicgstab: BiCGStab Krylov solver (default)\n"
    "- gmres: restarted GMRES, see 'restart'\n"
    "- flexgmres: flexible GMRES, for preconditioners that vary between "
    "iterations\n"
    "- cg: conjugate gradients (symmetric positive definite systems only)\n"
    "- loopsolver: plain preconditioned Richardson iteration\n"
    "- preconditioner2inverseoperator: apply the preconditioner exactly once\n"
    "- umfpack: sparse direct solver (serial runs only)\n"
    "- mixed-bicgstab: mixed-precision BiCGStab (serial, AVX2 builds, use with "
    "mixed-ilu0/mixed-dilu)\n"
    "- mixed-legacy: legacy mixed-precision solver (use with "
    "legacy-mixed-ilu0/legacy-mixed-dilu)";

constexpr std::string_view preconditionerTypeDoc =
    "Preconditioner. Allowed values:\n"
    "- paroverilu0: ILU(n), with n given by 'ilulevel', that handles the "
    "parallel overlap. ilu0 is ILU(0) and ilun is an alias of paroverilu0\n"
    "- duneilu: the DUNE ILU(n) implementation\n"
    "- dilu: diagonal ILU\n"
    "- jac, gs, sor, ssor: Jacobi, Gauss-Seidel, SOR and symmetric SOR, "
    "see 'relaxation' and 'repeats'\n"
    "- amg: DUNE aggregation based algebraic multigrid. kamg (Krylov cycle "
    "AMG) and famg (fast AMG) are available in serial runs only\n"
    "- cpr: two-stage Constrained Pressure Residual preconditioner, see "
    "'weight_type', 'finesmoother' and 'coarsesolver'. cprt uses the "
    "transposed weights\n"
    "- cprw: CPR with the wells included in the pressure system\n"
    "- system_cpr: CPR on the coupled reservoir-well system (CPU only), see "
    "the system_cpr/system_cprw presets of --linear-solver\n"
    "- hypre: BoomerAMG from Hypre (builds with Hypre only)\n"
    "- amgx: AMG from NVIDIA AMGX (builds with AMGX only)\n"
    "- mixed-ilu0, mixed-dilu: mixed-precision ILU0/DILU (serial, AVX2 builds, "
    "use with solver mixed-bicgstab)\n"
    "- legacy-mixed-ilu0, legacy-mixed-dilu: legacy mixed-precision variants "
    "(use with solver mixed-legacy)\n"
    "- tpsa_block: block preconditioner of the TPSA geomechanics solver\n"
    "With --linear-solver-accelerator=gpu, GPU preconditioners such as "
    "opmilu0, gpuilu0, gpudilu and gpujac are available as well.";

const std::map<std::string_view, KeyDoc>& solverDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "solver", { solverTypeDoc } },
        { "tol", { "Relative residual reduction at which the solver has "
                   "converged, i.e., the solver stops when ||r|| < tol * ||r0||. "
                   "Floating point number > 0." } },
        { "maxiter", { "Maximum number of iterations. Positive integer." } },
        { "verbosity", { "Amount of output from the solver: 0 = none, 1 = summary, "
                         "2 = output for each iteration." } },
        { "restart", { "Number of iterations after which GMRES restarts "
                       "(solver = gmres or flexgmres). Positive integer, default 15." } },
        { "preconditioner", { "Preconditioner of this solver." } },
    };
    return docs;
}

const std::map<std::string_view, KeyDoc>& preconditionerDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "type", { preconditionerTypeDoc } },
        { "relaxation", { "Relaxation (damping) factor of the ILU, Jacobi, "
                          "Gauss-Seidel or SOR preconditioner. Floating point "
                          "number, typically in (0, 1] for ILU and (0, 2) for SOR." } },
        { "ilulevel", { "Fill-in level n of ILU(n) for paroverilu0, ilun and duneilu. "
                        "Non-negative integer, default 0." } },
        { "milutype", { "Modified ILU variant for paroverilu0.",
                        "ilu (plain ILU, default), milu_1, milu_2, milu_3, milu_4 "
                        "(see --milu-variant)" } },
        { "redblack", { "Use red-black ordering in paroverilu0.", "true, false" } },
        { "reorder_spheres", { "Reorder the red-black ordered matrix in spheres "
                               "(paroverilu0 with redblack = true).", "true, false" } },
        { "resort", { "Resort the matrix in duneilu.", "true, false" } },
        { "repeats", { "Number of sweeps of the Jacobi, Gauss-Seidel, SOR or SSOR "
                       "preconditioner. Positive integer, default 1." } },
        { "verbosity", { "Amount of output from the preconditioner: 0 = none, "
                         "higher values give more output." } },
        { "weight_type", { "Weights used by CPR to form the pressure equation as a "
                           "combination of the conservation equations.",
                           "quasiimpes (from the diagonal of the Jacobian), "
                           "trueimpes (from the accumulation terms), "
                           "trueimpesanalytic (analytic true-IMPES weights, "
                           "not available for all models)" } },
        { "pre_smooth", { "Number of fine-scale smoother applications before the "
                          "pressure correction in CPR. Non-negative integer, Flow "
                          "uses --cpr-pre-smooth." } },
        { "post_smooth", { "Number of fine-scale smoother applications after the "
                           "pressure correction in CPR. Non-negative integer." } },
        { "add_wells", { "Include the wells in the CPR pressure system (cprw). "
                         "Should be true for cprw and false for cpr.", "true, false" } },
        { "use_well_weights", { "Use well-specific weights when forming the well "
                                "equations of the pressure system (cprw).",
                                "true, false" } },
        { "well_coarse_diagonal", { "How the diagonal of the pressure system is "
                                    "formed for multisegment wells (cprw) or wells "
                                    "(system_cpr).",
                                    "row_sum (minus the row sum), contract_d "
                                    "(contract the well diagonal block), and for "
                                    "system_cpr also auto (contract_d for single "
                                    "block wells, row_sum for multisegment wells)" } },
        { "weights_filename", { "File name used when writing the CPR weights for "
                                "debugging." } },
        { "finesmoother", { "Fine-scale smoother of CPR, applied to the full system. "
                            "Takes the same keys as a preconditioner." } },
        { "coarsesolver", { "Solver for the CPR pressure system. Takes the same "
                            "keys as the top-level solver." } },
        // System CPR
        { "well_weight_type", { "How the well equations are contracted to one "
                                "pressure equation per well (system_cpr with wells).",
                                "cellavg (mean reservoir weight of the perforated "
                                "cells, conservation equations only), cellblockavg "
                                "(the same mean per well block row), quasiimpes, "
                                "unit (the control or pressure equation as-is)" } },
        { "well_transfer", { "How the well unknowns take part in the pressure "
                             "stage transfer (system_cpr with wells).",
                             "full (restrict well residual, prolong bhp correction), "
                             "no_prolongation (restrict, discard bhp correction), "
                             "classic (neither, as in classic cprw)" } },
        { "well_identity_on_pressure_control", { "Give pressure controlled wells a "
                                                 "trivial pressure equation "
                                                 "(system_cpr).", "true, false" } },
        { "reservoir_smoother", { "Smoother for the reservoir part of the coupled "
                                  "system (system_cpr). Takes the same keys as a "
                                  "solver." } },
        { "reservoir_solver", { "Solver for the reservoir part of the coupled system "
                                "(system_cpr). Must use the cpr preconditioner." } },
        { "well_solver", { "Solver for the well part of the coupled system "
                           "(system_cpr). Takes the same keys as a solver." } },
        // TPSA
        { "disp_disp_solver", { "Solver for the displacement block (TPSA)." } },
        { "rot_rot_solver", { "Solver for the rotation block (TPSA)." } },
        { "spres_spres_solver", { "Solver for the solid pressure block (TPSA)." } },
    };
    return docs;
}

const std::map<std::string_view, KeyDoc>& duneAmgDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "alpha", { "Threshold for strong connections when aggregating: a "
                     "connection is strong if its weight exceeds alpha times the "
                     "largest weight. Floating point number in (0, 1)." } },
        { "beta", { "Connections with absolute value below beta are ignored when "
                    "aggregating. Floating point number >= 0." } },
        { "relaxation", { "Relaxation factor of the AMG smoother." } },
        { "iterations", { "Number of AMG cycles per preconditioner application. "
                          "Positive integer." } },
        { "coarsenTarget", { "Coarsening stops when a level has fewer unknowns than "
                             "this. Positive integer." } },
        { "maxlevel", { "Maximum number of levels in the AMG hierarchy. Positive "
                        "integer." } },
        { "pre_smooth", { "Number of pre-smoothing steps on each level. "
                          "Non-negative integer." } },
        { "post_smooth", { "Number of post-smoothing steps on each level. "
                           "Non-negative integer." } },
        { "smoother", { "Smoother used on each level of the AMG hierarchy.",
                        "ilu0, paroverilu0, ilun, dilu, jac, gs, sor, ssor" } },
        { "verbosity", { "Amount of output from the AMG setup: 0 = none, higher "
                         "values give more output." } },
        { "skip_isolated", { "Skip isolated nodes when aggregating.", "0, 1" } },
        { "accumulate", { "Redistribute the coarse levels to fewer processes in "
                          "parallel runs. Flow uses 1, since the pressure matrix "
                          "may be unsymmetric.",
                          "0 (no), 1 (to one process at once), "
                          "2 (successively)" } },
        { "prolongationdamping", { "Damping factor of the prolongation. Floating "
                                   "point number > 0." } },
        { "maxdistance", { "Maximum distance from the aggregate seed to any "
                           "aggregate member. Positive integer." } },
        { "maxconnectivity", { "Maximum connectivity of an aggregate. Positive "
                               "integer." } },
        { "maxaggsize", { "Maximum number of nodes in an aggregate. Positive "
                          "integer." } },
        { "minaggsize", { "Minimum number of nodes in an aggregate. Positive "
                          "integer." } },
    };
    return docs;
}

const std::map<std::string_view, KeyDoc>& hypreDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "print_level", { "BoomerAMG print level, 0 = no output." } },
        { "max_iter", { "Number of BoomerAMG cycles per preconditioner "
                        "application." } },
        { "strong_threshold", { "BoomerAMG strength threshold. Floating point "
                                "number in (0, 1)." } },
        { "agg_trunc_factor", { "Truncation factor of the aggressive coarsening "
                                "interpolation." } },
        { "interp_type", { "BoomerAMG interpolation type (integer, see the "
                           "HYPRE_BoomerAMGSetInterpType documentation)." } },
        { "max_levels", { "Maximum number of levels in the AMG hierarchy." } },
        { "tolerance", { "BoomerAMG convergence tolerance, 0 to always perform "
                         "max_iter cycles." } },
        { "use_gpu", { "Run BoomerAMG on the GPU.", "true, false" } },
        { "relax_type", { "BoomerAMG smoother (integer, see the "
                          "HYPRE_BoomerAMGSetRelaxType documentation)." } },
        { "coarsen_type", { "BoomerAMG coarsening algorithm (integer, see the "
                            "HYPRE_BoomerAMGSetCoarsenType documentation)." } },
        { "agg_num_levels", { "Number of levels with aggressive coarsening." } },
        { "agg_interp_type", { "Interpolation type on the aggressive coarsening "
                               "levels." } },
    };
    return docs;
}

const std::map<std::string_view, KeyDoc>& amgxDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "determinism_flag", { "Make AMGX results deterministic.", "0, 1" } },
        { "print_grid_stats", { "Print AMGX grid statistics.", "0, 1" } },
        { "print_solve_stats", { "Print AMGX solve statistics.", "0, 1" } },
        { "solver", { "AMGX solver (see the AMGX reference manual).",
                      "e.g., AMG" } },
        { "algorithm", { "AMGX AMG algorithm.", "CLASSICAL, AGGREGATION" } },
        { "interpolator", { "AMGX interpolation scheme (see the AMGX reference "
                            "manual).", "e.g., D1, D2" } },
        { "selector", { "AMGX coarse grid selector (see the AMGX reference "
                        "manual).", "e.g., PMIS, HMIS, SIZE_2" } },
        { "smoother", { "AMGX smoother (see the AMGX reference manual).",
                        "e.g., BLOCK_JACOBI, MULTICOLOR_GS, MULTICOLOR_DILU" } },
        { "presweeps", { "Number of pre-smoothing sweeps." } },
        { "postsweeps", { "Number of post-smoothing sweeps." } },
        { "strength_threshold", { "AMGX strength threshold." } },
        { "max_iters", { "Number of AMGX cycles per preconditioner application." } },
        { "setup_frequency", { "Number of preconditioner updates between full "
                               "AMGX setups." } },
    };
    return docs;
}

std::string parentOf(const std::string& key)
{
    const auto pos = key.rfind('.');
    return (pos == std::string::npos) ? std::string{} : key.substr(0, pos);
}

std::string leafOf(const std::string& key)
{
    const auto pos = key.rfind('.');
    return (pos == std::string::npos) ? key : key.substr(pos + 1);
}

std::string childKey(const std::string& parent, const std::string& leaf)
{
    return parent.empty() ? leaf : parent + '.' + leaf;
}

std::optional<KeyDoc> lookup(const std::map<std::string_view, KeyDoc>& docs,
                             const std::string& leaf)
{
    const auto pos = docs.find(leaf);
    if (pos == docs.end()) {
        return std::nullopt;
    }
    return pos->second;
}

/// Context dependent documentation of a linear solver configuration key.
///
/// The long lists of solvers and preconditioners are only included the
/// first time they are needed, as recorded in \p listed.
std::string linearSolverKeyComment(const Opm::PropertyTree& prm,
                                   const std::string& key,
                                   std::set<std::string>& listed)
{
    const auto parent = parentOf(key);
    const auto leaf = leafOf(key);
    const auto parentLeaf = leafOf(parent);
    const auto parentType = toLower(prm.get<std::string>(childKey(parent, "type"), ""));

    if (leaf == "type") {
        auto doc = std::string{};
        if (parentLeaf == "finesmoother") {
            doc = "Fine-scale smoother of CPR, typically paroverilu0, dilu or jac.";
        }
        else if ((parentLeaf == "preconditioner") &&
                 (leafOf(parentOf(parent)) == "coarsesolver"))
        {
            doc = "Preconditioner of the CPR pressure system, typically amg, "
                  "hypre or amgx.";
        }

        if (listed.insert("type").second) {
            return wrap(doc.empty() ? std::string { preconditionerTypeDoc }
                                    : doc + '\n' + std::string { preconditionerTypeDoc });
        }

        return wrap((doc.empty() ? std::string { "Preconditioner." } : doc) +
                    " See the first preconditioner 'type' above for all "
                    "allowed values.");
    }

    if (parentType == "hypre") {
        if (const auto doc = lookup(hypreDocs(), leaf)) {
            return formatDoc(*doc);
        }
    }
    else if (parentType == "amgx") {
        if (const auto doc = lookup(amgxDocs(), leaf)) {
            return formatDoc(*doc);
        }
    }
    else if ((parentType == "amg") || (parentType == "kamg") || (parentType == "famg")) {
        if (const auto doc = lookup(duneAmgDocs(), leaf)) {
            return formatDoc(*doc);
        }
    }

    // A table holding a 'solver' entry configures a solver, otherwise it
    // configures a preconditioner.
    const bool isSolverTable = prm.get_child_optional(childKey(parent, "solver")).has_value();
    const auto& primary = isSolverTable ? solverDocs() : preconditionerDocs();
    const auto& secondary = isSolverTable ? preconditionerDocs() : solverDocs();

    if (leaf == "solver") {
        if (listed.insert("solver").second) {
            return wrap(solverTypeDoc);
        }
        return wrap("Linear solver, e.g., bicgstab, gmres, loopsolver or umfpack. "
                    "See the first 'solver' above for all allowed values.");
    }

    if (parent.empty() && ((leaf == "tol") || (leaf == "maxiter"))) {
        const auto option = (leaf == "tol") ? "reduction" : "max-iter";
        return formatDoc(*lookup(solverDocs(), leaf)) + '\n' +
            wrap(fmt::format("Flow uses the value of --linear-solver-{} (or the "
                             "corresponding NLDD/TPSA option) for presets, "
                             "except that CPR presets default to {} unless the "
                             "option is given.", option,
                             (leaf == "tol") ? "0.005" : "20"));
    }

    if (const auto doc = lookup(primary, leaf)) {
        return formatDoc(*doc);
    }
    if (const auto doc = lookup(secondary, leaf)) {
        return formatDoc(*doc);
    }

    return {};
}

/// Keys holding floating point values.  The property tree stores, e.g.,
/// 1.0 as "1", which would otherwise be written as a TOML integer.
bool isFloatKey(const std::string& leaf)
{
    static const auto keys = std::set<std::string> {
        "tol", "relaxation", "alpha", "beta", "prolongationdamping",
        "strong_threshold", "agg_trunc_factor", "tolerance",
        "strength_threshold", "min_reduction",
    };
    return keys.count(leaf) > 0;
}

void writeFloatsAsFloats(Opm::PropertyTree& root,
                         const Opm::PropertyTree& node,
                         const std::string& path)
{
    for (const auto& key : node.get_child_keys()) {
        if (key.empty()) {
            continue;
        }

        const auto child = node.get_child(key);
        const auto fullKey = childKey(path, key);
        if (! child.get_child_keys().empty()) {
            writeFloatsAsFloats(root, child, fullKey);
            continue;
        }

        const auto value = child.get<std::string>("");
        const auto digits = std::string_view { value }.substr(value.starts_with('-') ? 1 : 0);
        if (isFloatKey(key) && !digits.empty() &&
            std::ranges::all_of(digits, [](unsigned char c) { return std::isdigit(c); }))
        {
            root.put(fullKey, value + ".0");
        }
    }
}

std::string linearSolverHeader(const std::string& optionName,
                               const std::string& origin)
{
    const auto text = fmt::format(
        "Linear solver configuration for OPM Flow, in TOML format "
        "(https://toml.io).\n"
        "\n"
        "This file was generated by 'flow --write-toml-example' and holds the "
        "configuration used by {0}. Use it with\n"
        "\n"
        "    flow {1}=<this file>.toml CASE.DATA\n"
        "\n"
        "When a configuration file is given, {1} only selects the file, and "
        "related options such as --linear-solver-reduction, "
        "--linear-solver-max-iter, --ilu-relaxation and --cpr-pre-smooth are "
        "ignored; the values in the file are used instead. Generate the "
        "configuration of another preset by combining --write-toml-example "
        "with {1}=<preset>, e.g., {1}=ilu0. Converting a JSON configuration "
        "to TOML works the same way, with {1}=<file>.json.\n"
        "\n"
        "Tables, e.g. [preconditioner], configure nested components. Keys that "
        "are left out take built-in defaults.",
        origin, optionName);

    return commentBlock(text);
}

// ===========================================================================
// Zoltan
// ===========================================================================

struct OptionalEntry
{
    std::string_view key;
    std::string_view exampleValue;
    KeyDoc doc;
};

const std::map<std::string_view, KeyDoc>& zoltanDocs()
{
    static const auto docs = std::map<std::string_view, KeyDoc> {
        { "LB_METHOD", { "Load-balancing method. Flow provides the grid as a graph, "
                         "so only the graph based methods are meaningful.",
                         "\"GRAPH\", \"HYPERGRAPH\"" } },
        { "GRAPH_PACKAGE", { "Graph partitioning package used when LB_METHOD = "
                             "\"GRAPH\". Scotch and ParMETIS require Zoltan to be "
                             "built with support for them.",
                             "\"PHG\", \"Scotch\", \"ParMETIS\"" } },
        { "HYPERGRAPH_PACKAGE", { "Hypergraph partitioning package used when "
                                  "LB_METHOD = \"HYPERGRAPH\".", "\"PHG\"" } },
        { "PHG_EDGE_SIZE_THRESHOLD", { "Hyperedges connecting more than this "
                                       "fraction of the vertices are ignored "
                                       "(0 removes all, 1 removes none). Flow "
                                       "uses --zoltan-phg-edge-size-threshold for "
                                       "the graph and hypergraph presets.",
                                       "floating point number in [0, 1]" } },
        { "IMBALANCE_TOL", { "Tolerated load imbalance, e.g., 1.1 allows 10% "
                             "imbalance. Flow uses --imbalance-tol unless set here.",
                             "floating point number >= 1" } },
        { "DEBUG_LEVEL", { "Amount of output from Zoltan.",
                           "integer from 0 (none, Flow's default) to 10" } },
        { "LB_APPROACH", { "Partitioning approach.",
                           "\"PARTITION\" (Flow's default), \"REPARTITION\", "
                           "\"REFINE\"" } },
        { "CHECK_GRAPH", { "Level of consistency checks of the graph.",
                           "0 (none), 1 (some), 2 (full)" } },
        { "PHG_COARSENING_METHOD", { "Coarsening method of the PHG partitioner.",
                                     "\"agg\", \"ipm\", \"l-ipm\", \"a-ipm\"" } },
        { "PHG_REPART_MULTIPLIER", { "Weight of migration cost relative to "
                                     "communication cost for LB_APPROACH = "
                                     "\"REPARTITION\".",
                                     "floating point number > 0" } },
        { "EnvelopeWellLayers", { "OPM specific parameter, only used with "
                                  "--partition-method=zoltanwell: number of "
                                  "layers of cells around each well that are "
                                  "kept on the same process as the well.",
                                  "non-negative integer, default 0" } },
    };
    return docs;
}

const std::vector<OptionalEntry>& zoltanOptionalEntries()
{
    static const auto entries = std::vector<OptionalEntry> {
        { "IMBALANCE_TOL", "", {} },
        { "DEBUG_LEVEL", "0", {} },
        { "LB_APPROACH", "\"PARTITION\"", {} },
        { "CHECK_GRAPH", "2", {} },
        { "PHG_COARSENING_METHOD", "\"agg\"", {} },
        { "EnvelopeWellLayers", "1", {} },
    };
    return entries;
}

// ===========================================================================
// METIS
// ===========================================================================

const std::vector<OptionalEntry>& metisEntries()
{
    static const auto entries = std::vector<OptionalEntry> {
        { "METIS_OPTION_PTYPE", "\"METIS_PTYPE_KWAY\"",
          { "Partitioning method. By default, Flow uses recursive bisection if "
            "the number of processes is a power of two and less than 65, and "
            "k-way partitioning otherwise.",
            "\"METIS_PTYPE_RB\" (recursive bisection), \"METIS_PTYPE_KWAY\" "
            "(k-way)" } },
        { "METIS_OPTION_OBJTYPE", "\"METIS_OBJTYPE_CUT\"",
          { "Objective of the k-way partitioning.",
            "\"METIS_OBJTYPE_CUT\" (minimise edge cut, default), "
            "\"METIS_OBJTYPE_VOL\" (minimise communication volume)" } },
        { "METIS_OPTION_CTYPE", "\"METIS_CTYPE_SHEM\"",
          { "Matching scheme used when coarsening.",
            "\"METIS_CTYPE_RM\" (random matching), \"METIS_CTYPE_SHEM\" "
            "(sorted heavy-edge matching, default)" } },
        { "METIS_OPTION_IPTYPE", "\"METIS_IPTYPE_GROW\"",
          { "Algorithm for the initial partitioning.",
            "\"METIS_IPTYPE_GROW\", \"METIS_IPTYPE_RANDOM\", "
            "\"METIS_IPTYPE_EDGE\", \"METIS_IPTYPE_NODE\"" } },
        { "METIS_OPTION_RTYPE", "\"METIS_RTYPE_FM\"",
          { "Algorithm for refinement.",
            "\"METIS_RTYPE_FM\", \"METIS_RTYPE_GREEDY\", "
            "\"METIS_RTYPE_SEP2SIDED\", \"METIS_RTYPE_SEP1SIDED\"" } },
        { "METIS_OPTION_NCUTS", "1",
          { "Number of different partitionings to compute; the best one is "
            "used.", "positive integer, default 1" } },
        { "METIS_OPTION_NITER", "10",
          { "Number of refinement iterations on each level.",
            "positive integer, default 10" } },
        { "METIS_OPTION_SEED", "0",
          { "Seed of the random number generator.", "integer" } },
        { "METIS_OPTION_MINCONN", "0",
          { "Minimise the maximum connectivity between partitions (k-way only).",
            "0, 1" } },
        { "METIS_OPTION_CONTIG", "0",
          { "Force contiguous partitions (k-way only).", "0, 1" } },
        { "METIS_OPTION_UFACTOR", "30",
          { "Maximum allowed load imbalance, as (1 + x)/1000. Note that Flow "
            "also passes --imbalance-tol to METIS.",
            "positive integer, METIS defaults are 1 (recursive bisection) and "
            "30 (k-way)" } },
        { "METIS_OPTION_NO2HOP", "0",
          { "Do not perform 2-hop matching when coarsening.", "0, 1" } },
        { "METIS_OPTION_DBGLVL", "\"METIS_DBG_INFO\"",
          { "Amount of diagnostic output from METIS.",
            "an integer bit mask or one of \"METIS_DBG_INFO\", "
            "\"METIS_DBG_TIME\", \"METIS_DBG_COARSEN\", \"METIS_DBG_REFINE\", "
            "\"METIS_DBG_IPART\", \"METIS_DBG_MOVEINFO\", \"METIS_DBG_SEPINFO\", "
            "\"METIS_DBG_CONNINFO\", \"METIS_DBG_CONTIGINFO\"" } },
    };
    return entries;
}

// ===========================================================================
// Common helpers
// ===========================================================================

Opm::PropertyTree flatTree(const std::map<std::string, std::string>& params)
{
    auto tree = Opm::PropertyTree{};
    for (const auto& [key, value] : params) {
        tree.put(key, value);
    }
    return tree;
}

std::string optionalEntries(const std::vector<OptionalEntry>& entries,
                            const std::map<std::string, std::string>& active,
                            const std::map<std::string_view, KeyDoc>& docs,
                            const std::map<std::string_view, std::string>& exampleOverrides = {})
{
    std::string out;
    for (const auto& entry : entries) {
        if (active.count(std::string{entry.key}) > 0) {
            continue;
        }

        const auto docPos = docs.find(entry.key);
        const auto& doc = (docPos != docs.end()) ? docPos->second : entry.doc;
        const auto custom = exampleOverrides.find(entry.key);
        const auto value = (custom != exampleOverrides.end())
            ? std::string_view { custom->second } : entry.exampleValue;

        out += '\n';
        out += commentBlock(formatDoc(doc));
        out += fmt::format("# {} = {}\n", entry.key, value);
    }
    return out;
}

// ===========================================================================
// Generation from run-time parameters
// ===========================================================================

std::string linearSolverExampleFromParameters(const std::string& name)
{
    std::string optionName;
    std::string conf;
    Opm::PropertyTree prm;

    if (name == "tpsa-linear-solver") {
        Opm::TpsaLinearSolverParameters p;
        p.init();
        optionName = "--tpsa-linear-solver";
        conf = p.linsolver_;
        prm = Opm::setupPropertyTree(p, false, false, /*tpsaSetup=*/true);
    }
    else {
        Opm::FlowLinearSolverParameters p;
        const bool nldd = (name == "nldd-local-linear-solver");
        p.is_nldd_local_solver_ = nldd;
        p.init(false);
        optionName = nldd ? "--nldd-local-linear-solver" : "--linear-solver";
        conf = p.linsolver_;

        if (conf == "hybrid") {
            throw std::invalid_argument {
                "The 'hybrid' linear solver uses two configurations and cannot be "
                "written as a single configuration file."
            };
        }

        const bool maxIterSet = nldd
            ? Opm::Parameters::IsSet<Opm::Parameters::NlddLocalLinearSolverMaxIter>()
            : Opm::Parameters::IsSet<Opm::Parameters::LinearSolverMaxIter>();
        const bool reductionSet = nldd
            ? Opm::Parameters::IsSet<Opm::Parameters::NlddLocalLinearSolverReduction>()
            : Opm::Parameters::IsSet<Opm::Parameters::LinearSolverReduction>();
        prm = Opm::setupPropertyTree(p, maxIterSet, reductionSet);
    }

    const auto origin = Opm::PropertyTree::isConfigFileName(conf)
        ? fmt::format("{}={} (read from file)", optionName, conf)
        : fmt::format("{}={}", optionName, conf);

    return Opm::linearSolverTomlExample(prm, optionName, origin);
}

#if HAVE_MPI
std::string zoltanExampleFromParameters()
{
    const auto conf = Opm::Parameters::Get<Opm::Parameters::ZoltanParams>();
    const auto threshold = Opm::Parameters::Get<Opm::Parameters::ZoltanPhgEdgeSizeThreshold>();
    const auto imbalanceTol =
        Opm::Parameters::IsSet<Opm::Parameters::ZoltanImbalanceTol<double>>(false)
        ? Opm::Parameters::Get<Opm::Parameters::ZoltanImbalanceTol<double>>(false)
        : Opm::Parameters::Get<Opm::Parameters::ImbalanceTol<double>>(false);

    return Opm::zoltanParamsTomlExample(Opm::setupZoltanParams(conf, threshold),
                                        fmt::format("--zoltan-params={}", conf),
                                        imbalanceTol);
}

std::string metisExampleFromParameters()
{
    const auto conf = Opm::Parameters::Get<Opm::Parameters::MetisParams>();
    return Opm::metisParamsTomlExample(Opm::setupMetisParams(conf),
                                       fmt::format("--metis-params={}", conf));
}
#endif // HAVE_MPI

std::string exampleFromParameters(const std::string& name)
{
    if ((name == "linear-solver") ||
        (name == "nldd-local-linear-solver") ||
        (name == "tpsa-linear-solver"))
    {
        return linearSolverExampleFromParameters(name);
    }

#if HAVE_MPI
    if (name == "zoltan-params") {
        return zoltanExampleFromParameters();
    }

    if (name == "metis-params") {
        return metisExampleFromParameters();
    }
#endif // HAVE_MPI

    if (name == "hybrid-newton") {
        return Opm::hybridNewtonTomlExample();
    }

    throw std::invalid_argument {
        fmt::format("Unknown example name '{}' in --write-toml-example", name)
    };
}

std::vector<std::string> parseExampleNames(const std::string& spec)
{
    const auto known = Opm::tomlExampleNames();

    std::vector<std::string> names;
    std::istringstream items { spec };
    for (std::string item; std::getline(items, item, ','); ) {
        const auto first = item.find_first_not_of(" \t");
        if (first == std::string::npos) {
            continue;
        }
        item = item.substr(first, item.find_last_not_of(" \t") - first + 1);

        if (item == "all") {
            for (const auto& name : known) {
                if (std::ranges::find(names, name) == names.end()) {
                    names.push_back(name);
                }
            }
            continue;
        }

        if (std::ranges::find(known, item) == known.end()) {
            std::string list;
            for (const auto& name : known) {
                list += name + ", ";
            }
            throw std::invalid_argument {
                fmt::format("Unknown example name '{}' in --write-toml-example. "
                            "Valid names are: {}or all.", item, list)
            };
        }

        if (std::ranges::find(names, item) == names.end()) {
            names.push_back(item);
        }
    }

    if (names.empty()) {
        throw std::invalid_argument { "No example names given in --write-toml-example" };
    }

    return names;
}

} // Anonymous namespace

namespace Opm {

void registerTomlExampleParameters()
{
    std::string names;
    for (const auto& name : tomlExampleNames()) {
        names += name + ", ";
    }

    Parameters::Register<Parameters::WriteTomlExample>
        (fmt::format("Write documented example TOML configuration files to the "
                     "current directory and exit. Comma-separated list of: {}or "
                     "all. Each example NAME is written to NAME.toml, and "
                     "existing files are not overwritten. The examples reflect the "
                     "values of the corresponding options, e.g., "
                     "--linear-solver=ilu0 --write-toml-example=linear-solver "
                     "writes the configuration of the ilu0 preset.", names).c_str());
}

std::vector<std::string> tomlExampleNames()
{
    return {
        "linear-solver",
        "nldd-local-linear-solver",
        "tpsa-linear-solver",
#if HAVE_MPI
        "zoltan-params",
        "metis-params",
#endif
        "hybrid-newton",
    };
}

std::vector<std::filesystem::path>
writeTomlExamples(const std::string& spec,
                  const std::filesystem::path& directory)
{
    // Generate all examples before writing anything, so that an error does
    // not leave a partial set of files behind.
    std::vector<std::pair<std::filesystem::path, std::string>> examples;
    for (const auto& name : parseExampleNames(spec)) {
        auto file = directory / (name + ".toml");
        if (std::filesystem::exists(file)) {
            throw std::invalid_argument {
                fmt::format("The file {} exists already. Remove or rename it to "
                            "write a new example.", file.generic_string())
            };
        }
        examples.emplace_back(std::move(file), exampleFromParameters(name));
    }

    std::vector<std::filesystem::path> written;
    for (const auto& [file, content] : examples) {
        std::ofstream os { file };
        os << content;
        // Close explicitly, so that failures to flush buffered output, e.g.,
        // because the disk is full, are detected.
        os.close();
        if (! os) {
            throw std::runtime_error {
                fmt::format("Failed to write {}", file.generic_string())
            };
        }
        written.push_back(file);
    }

    return written;
}

std::string linearSolverTomlExample(const PropertyTree& prm,
                                    const std::string& optionName,
                                    const std::string& origin)
{
    auto tree = prm;
    writeFloatsAsFloats(tree, prm, "");

    std::set<std::string> listed;
    std::ostringstream os;
    os << linearSolverHeader(optionName, origin) << '\n';
    tree.write_toml(os, [&tree, &listed](const std::string& key)
    {
        return linearSolverKeyComment(tree, key, listed);
    });
    return os.str();
}

std::string zoltanParamsTomlExample(const std::map<std::string, std::string>& params,
                                    const std::string& origin,
                                    const double imbalanceTol)
{
    const auto header = fmt::format(
        "Zoltan partitioning parameters for OPM Flow, in TOML format "
        "(https://toml.io).\n"
        "\n"
        "This file was generated by 'flow --write-toml-example' and holds the "
        "parameters used by {}. Use it with\n"
        "\n"
        "    flow --partition-method=zoltan --zoltan-params=<this file>.toml "
        "CASE.DATA\n"
        "\n"
        "The parameters are used with --partition-method=zoltan and "
        "--partition-method=zoltanwell (the default). Each key/value pair is "
        "passed to Zoltan_Set_Param(), after Flow's own defaults, so the values "
        "here take precedence. Only top-level keys are used. Values may be "
        "strings or numbers; use 0/1 rather than true/false. See "
        "https://sandialabs.github.io/Zoltan/ug_html/ug.html for all "
        "parameters.\n"
        "\n"
        "Commented out entries below show other commonly used parameters, with "
        "example values.",
        origin);

    auto tree = flatTree(params);

    std::ostringstream os;
    os << commentBlock(header) << '\n';
    tree.write_toml(os, [](const std::string& key)
    {
        const auto pos = zoltanDocs().find(key);
        return (pos != zoltanDocs().end()) ? formatDoc(pos->second) : std::string{};
    });

    os << optionalEntries(zoltanOptionalEntries(), params, zoltanDocs(),
                          { { "IMBALANCE_TOL", fmt::format("{}", imbalanceTol) } });

    return os.str();
}

std::string metisParamsTomlExample(const std::map<std::string, std::string>& params,
                                   const std::string& origin)
{
    const auto header = fmt::format(
        "METIS partitioning parameters for OPM Flow, in TOML format "
        "(https://toml.io).\n"
        "\n"
        "This file was generated by 'flow --write-toml-example' and holds the "
        "parameters used by {}. Use it with\n"
        "\n"
        "    flow --partition-method=metis --metis-params=<this file>.toml "
        "CASE.DATA\n"
        "\n"
        "Only top-level keys are used. Unset options take the METIS defaults. "
        "Options with symbolic values take either the name of the METIS "
        "constant as a string or an integer. The options are not available if "
        "METIS is provided by Scotch. See the METIS manual for details.\n"
        "\n"
        "Commented out entries below show all supported options, with example "
        "values.",
        origin);

    auto tree = flatTree(params);

    std::map<std::string_view, KeyDoc> docs;
    for (const auto& entry : metisEntries()) {
        docs.emplace(entry.key, entry.doc);
    }

    std::ostringstream os;
    os << commentBlock(header);
    if (! params.empty()) {
        os << '\n';
    }
    tree.write_toml(os, [&docs](const std::string& key)
    {
        const auto pos = docs.find(key);
        return (pos != docs.end()) ? formatDoc(pos->second) : std::string{};
    });

    os << optionalEntries(metisEntries(), params, docs);

    return os.str();
}

std::string hybridNewtonTomlExample()
{
    const auto header = std::string_view {
        "Hybrid Newton configuration for OPM Flow, in TOML format "
        "(https://toml.io).\n"
        "\n"
        "This file was generated by 'flow --write-toml-example'. Use it with\n"
        "\n"
        "    flow --use-hybrid-newton=true "
        "--hybrid-newton-config-file=<this file>.toml CASE.DATA\n"
        "\n"
        "Hybrid Newton uses trained neural network models to predict the "
        "initial guess of the nonlinear (Newton) solver at the start of "
        "selected time steps. Only three-phase black-oil models are "
        "supported.\n"
        "\n"
        "The file holds one table per model. The table names, here "
        "'model_1', are arbitrary. All file names and values below are "
        "examples and must be adapted to the trained model; there are no "
        "defaults."
    };

    const auto features = std::string_view {
        "Features. Input features are evaluated, transformed and scaled in the "
        "order given here, which must match the order expected by the model. "
        "Each feature is a table, named after the quantity, holding the "
        "optional keys\n"
        "- feature_engineering: transformation applied before scaling, one of "
        "\"none\" (default), \"log\", \"log10\", \"log1p\"\n"
        "- scaling_params: table with either 'mean' and 'std' (standard "
        "scaling, (x - mean)/std) or 'min' and 'max' (min-max scaling, "
        "(x - min)/(max - min)). No scaling if left out.\n"
        "\n"
        "Supported input features: PRESSURE, SWAT, SOIL, SGAS, RS, RV, PERMX "
        "(one value per cell) and TIMESTEP (current time step size in "
        "seconds, one value). Values are in the unit system of the deck.\n"
        "\n"
        "Supported output features: PRESSURE, SWAT, SOIL, SGAS, RS, RV. Prefix "
        "a name with DELTA_, e.g. DELTA_PRESSURE, if the model predicts the "
        "change from the current state rather than the new state. Missing "
        "saturations are computed from the others. RS and RV require models "
        "with dissolved gas/vaporised oil."
    };

    std::string out = commentBlock(header);
    out += '\n';
    out += "[model_1]\n";
    out += commentBlock("Neural network model file, in the OPM ML model format.");
    out += "model_path = \"models/model_1.model\"\n\n";
    out += commentBlock("Text file with the (0-based) active cell indices the model "
                        "applies to, one index per line. Empty lines and lines "
                        "starting with '#' are ignored.");
    out += "cell_indices_file = \"models/model_1_active_cells.txt\"\n\n";
    out += commentBlock("Simulation time(s), in seconds since the start of the "
                        "simulation, at which the model is applied. One value: "
                        "applied at the time step starting at that time. Two "
                        "values [start, end]: applied at every time step starting "
                        "in the interval.");
    out += "apply_times = [0.0, 8640000.0]\n\n";
    out += commentBlock(features);
    out += '\n';
    out += "[model_1.features.inputs.PRESSURE]\n"
           "feature_engineering = \"none\"\n"
           "scaling_params = { mean = 250.0, std = 50.0 }\n\n"
           "[model_1.features.inputs.SGAS]\n"
           "scaling_params = { min = 0.0, max = 1.0 }\n\n"
           "[model_1.features.inputs.PERMX]\n"
           "feature_engineering = \"log10\"\n"
           "scaling_params = { mean = 2.0, std = 0.5 }\n\n"
           "[model_1.features.inputs.TIMESTEP]\n"
           "feature_engineering = \"log1p\"\n\n"
           "[model_1.features.outputs.DELTA_PRESSURE]\n"
           "scaling_params = { mean = 0.0, std = 10.0 }\n\n"
           "[model_1.features.outputs.SGAS]\n"
           "scaling_params = { min = 0.0, max = 1.0 }\n";

    return out;
}

} // namespace Opm
