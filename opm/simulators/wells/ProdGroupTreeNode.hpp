/*
  Copyright 2026 Sintef Digital

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

#ifndef OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED
#define OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <array>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace Opm {

/// Represents the control status of a node in the balanced production group tree.
    enum class ProdNodeModeCategory {
    Group,          ///< Node is controlled by higher group
    Individual,     ///< Node is controlled by one of its own individual limits
    None,           ///< Group node with all children at individual or none category
    Transparent,    ///< Group node with no guide rate producing less than its limits
};

/// Node type: well vs. group
enum class ProdNodeType {
    Well,
    Group
};

/// One node in the production group tree used by the balancing algorithm.
/// The algorithm always uses 3-component [oil, water, gas] rate vectors;
/// missing phases are set to zero.
/// Rates follow the sign convention: negative = production.
template<class Scalar>
struct ProdGroupTreeNode {
    // ---- Identity --------------------------------------------------------
    std::string name;
    ProdNodeType type{ProdNodeType::Well};

    // ---- Tree topology ---------------------------------------------------
    std::string parent;                ///< parent node name; empty for FIELD
    std::vector<std::string> children; ///< child node names

    // ---- Group control status --------------------------------------------
    bool availableForGroupControl{true};      ///< well: isAvailableForGroupControl()
                                              ///< group: productionGroupControlAvailable()
    
    // A group/well that has no ancestor with production limits, needs to be balanced
    // individually if it has limits itself (it may still be flagged as available for 
    // group control).
    bool hasLimitedAncestor{true};

    ProdNodeModeCategory modeCategory{ProdNodeModeCategory::Group};

    /// Which (individual or group) control mode is currently active / most limiting.
    Well::ProducerCMode mode{Well::ProducerCMode::CMODE_UNDEFINED};

    /// Preferred control / control given in schedule for this node
    Group::ProductionCMode preferredMode{Group::ProductionCMode::NONE};

    // ---- Rates -----------------------------------------------------------
    /// [oil, water, gas] surface rates; negative = production.
    std::array<Scalar, 3> rates{};

    /// Snapshot of rates at tree-build time (beginning of timestep); never modified
    /// during balancing. Used for guide-rate lookups to match FractionCalculator's
    /// behaviour of always scaling guide rates against beginning-of-timestep rates.
    std::array<Scalar, 3> initialRates{};

    // ---- Individual limits -----------------------------------------------
    /// Maps each active limit type to the scalar value the algorithm compares
    /// against the corresponding projection of the node's rates.
    ///  ORAT  → oil limit (compared to -rates[0])
    ///  WRAT  → water limit (compared to -rates[1])
    ///  GRAT  → gas limit (compared to -rates[2])
    ///  LRAT  → liquid limit (compared to -(rates[0]+rates[1]))
    ///  RESV  → reservoir-volume limit (compared to sum_p -rates[p]*resv_coeff[p])
    ///  BHP/THP → equivalent total rate derived from IPR at the BHP limit
    ///          (for wells only; groups do not have a BHP limit)
    std::map<Well::ProducerCMode, Scalar> Limits;

    // ---- Group target (set by parent during target assignment) -----------
    struct GroupTarget {
        std::string groupName;
        Group::ProductionCMode ctrlMode{Group::ProductionCMode::NONE};
        Scalar value{0};
        Scalar guideRate{0};
        Scalar guideRateRatio{1};
    };
    GroupTarget groupTarget;

    // ---- RESV conversion coefficients (wells only) -----------------------
    /// convert_coeff[i] for each active phase so that
    ///   resv_rate = sum_p -rates[p] * resvCoeff[p].
    /// Stored in canonical [oil, water, gas] order; zero for inactive phases.
    std::array<Scalar, 3> resvCoeff{};

    // ---- Efficiency factor -----------------------------------------------
    /// Own efficiency factor: WEFAC * efficiency_scaling_factor for wells,
    /// GEFAC for groups. Rates flowing upward to the parent are multiplied
    /// by this factor; targets flowing downward from the parent are divided
    /// by this factor before being applied to this node.
    Scalar efficiencyFactor{1};

    // ---- Fields for new sorting-based balancing algorithm ----------------
    /// Whether this node participates in guide-rate balancing (false for transparent groups)
    bool hasGuideRate{false};

    /// Accumulated rates during balancing (sum of children's rates)
    std::array<Scalar, 3> rateSums{};

    /// Accumulated guide rate sums during balancing
    std::array<Scalar, 3> guideRateSums{};

    /// Flag for fallback mechanism (when mode fraction is too small)
    bool useFallback{false};

    /// Visited flag for free-path checking
    bool visited{false};

    /// Whether this node is a satellite group with fixed rates from GSATPROD
    bool isSatellite{false};

    /// Whether this subtree is balanced
    bool isBalanced{false};

    /// Iteration counter for this node (cumulative across all calls to balanceGroupTree)
    int balancingCount{0};

    /// Number of while-loop iterations in the most recent call to balanceGroupTree on this node
    int lastIterationCount{0};

    /// Fallback group target (for non-preferred mode)
    GroupTarget groupTargetFallback;
};

} // namespace Opm

#endif // OPM_PROD_GROUP_TREE_NODE_HEADER_INCLUDED
