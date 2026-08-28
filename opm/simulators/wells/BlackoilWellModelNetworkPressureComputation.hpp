/*
  Copyright 2020-2026 Equinor ASA

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

#ifndef OPM_BLACKOIL_WELL_MODEL_NETWORK_PRESSURE_COMPUTATION_HPP
#define OPM_BLACKOIL_WELL_MODEL_NETWORK_PRESSURE_COMPUTATION_HPP

#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/output/data/Groups.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/network/NetworkNodePressureUpdater.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <stack>
#include <string>
#include <vector>

namespace Opm {

/// Result of a single network branch VFP lookup.
template<typename Scalar>
struct NetworkBranchPressure
{
    Scalar pressure{0.0};
    // False when the table has no solution at this point (zero-filled cells give
    // bhp <= 1 atm); the pressure must then not be used as a node pressure.
    bool valid{true};
    // True when the flow rate or upstream pressure had to be clamped to the table axes.
    bool clamped{false};
};

namespace detail {
    /// Clamp the VFP lookup point to the table axes; the tables must not be extrapolated
    /// for network branches (a zero-filled tail extrapolates to negative pressures).
    /// Rates are scaled uniformly so that WFR/GFR fractions are preserved.
    template<typename Scalar, typename IndexTraits, typename Table>
    bool clampToTableAxes(const Table& table, std::vector<Scalar>& rates, Scalar& up_press)
    {
        bool clamped = false;
        const auto& thp_axis = table.getTHPAxis();
        const Scalar thp_lo = thp_axis.front();
        const Scalar thp_hi = thp_axis.back();
        if (up_press < thp_lo || up_press > thp_hi) {
            up_press = std::clamp(up_press, thp_lo, thp_hi);
            clamped = true;
        }
        const auto& flo_axis = table.getFloAxis();
        const Scalar flo = std::abs(getFlo(table,
                                           rates[IndexTraits::waterPhaseIdx],
                                           rates[IndexTraits::oilPhaseIdx],
                                           rates[IndexTraits::gasPhaseIdx]));
        const Scalar flo_hi = flo_axis.back();
        if (flo > flo_hi && flo > 0.0) {
            const Scalar s = flo_hi / flo;
            std::ranges::transform(rates, rates.begin(), [s](const auto r) { return s * r; });
            clamped = true;
        }
        return clamped;
    }
} // namespace detail

/// @brief  Helper class to insulate the NetworkPressureComputation class from
///         the differences between production and injection VFP tables.
template<typename Scalar, typename IndexTraits, typename VfpProperties>
struct NetworkVfpPressureCalculator;

// Production specialization.
template<typename Scalar, typename IndexTraits>
struct NetworkVfpPressureCalculator<Scalar, IndexTraits, VFPProdProperties<Scalar>>
{
    static void prepareRates(std::vector<Scalar>& rates)
    {
        // Network rates are positive, while production VFP expects negative rates.
        std::ranges::transform(rates, rates.begin(), [](const auto r) { return -r; });
    }

    template <class GroupState>
    static bool hasLeafNodeRate(const GroupState& group_state,
                                const std::string& node)
    {
        return group_state.has_network_leaf_node_production_rates(node);
    }

    template <class GroupState>
    static const std::vector<Scalar>
    leafNodeRate(const GroupState& group_state,
                 const std::string& node)
    {
        return group_state.network_leaf_node_production_rates(node);
    }

    template<typename Branch>
    static NetworkBranchPressure<Scalar> compute(const VFPProdProperties<Scalar>& vfp_props,
                                                 const int table_id,
                                                 std::vector<Scalar> rates,
                                                 Scalar up_press,
                                                 const Branch& upbranch,
                                                 const UnitSystem& unit_system)
    {
        // NB! ALQ in extended network is never implicitly the gas lift rate (GRAT), i.e., the
        //     gas lift rates only enters the network pressure calculations through the rates
        //     (e.g., in GOR calculations) unless a branch ALQ is set in BRANPROP.
        const auto& table = vfp_props.getTable(table_id);
        const auto alq_type = table.getALQType();
        const auto dimension = VFPProdTable::ALQDimension(alq_type, unit_system);
        const Scalar alq = upbranch.alq_value(dimension).value_or(0.0);

        NetworkBranchPressure<Scalar> result;
        result.clamped = detail::clampToTableAxes<Scalar, IndexTraits>(table, rates, up_press);
        result.pressure = vfp_props.bhp(table_id,
                                        rates[IndexTraits::waterPhaseIdx],
                                        rates[IndexTraits::oilPhaseIdx],
                                        rates[IndexTraits::gasPhaseIdx],
                                        up_press,
                                        alq,
                                        0.0, // explicit_wfr
                                        0.0, // explicit_gfr
                                        false); // use_expvfp we dont support explicit lookup
        result.valid = result.pressure > unit::atm;
        return result;
    }
};

// Injection specialization.
template<typename Scalar, typename IndexTraits>
struct NetworkVfpPressureCalculator<Scalar, IndexTraits, VFPInjProperties<Scalar>>
{
    static void prepareRates(std::vector<Scalar>&)
    {
    }

    template <class GroupState>
    static bool hasLeafNodeRate(const GroupState& group_state,
                                const std::string& node)
    {
        return group_state.has_network_leaf_node_injection_rates(node);
    }

    template <class GroupState>
    static const std::vector<Scalar>
    leafNodeRate(const GroupState& group_state,
                 const std::string& node)
    {
        return group_state.network_leaf_node_injection_rates(node);
    }

    template<typename Branch>
    static NetworkBranchPressure<Scalar> compute(const VFPInjProperties<Scalar>& vfp_props,
                                                 const int table_id,
                                                 std::vector<Scalar> rates,
                                                 Scalar up_press,
                                                 const Branch&,
                                                 const UnitSystem&)
    {
        NetworkBranchPressure<Scalar> result;
        result.clamped = detail::clampToTableAxes<Scalar, IndexTraits>(vfp_props.getTable(table_id), rates, up_press);
        result.pressure = vfp_props.bhp(table_id,
                                        rates[IndexTraits::waterPhaseIdx],
                                        rates[IndexTraits::oilPhaseIdx],
                                        rates[IndexTraits::gasPhaseIdx],
                                        up_press);
        result.valid = result.pressure > unit::atm;
        return result;
    }
};

/// @brief  Class to compute network pressures using VFP tables, given flow rates
///         for each group and fixed pressures at network roots.
///         Optionally, the ALQ values from wells can be included in the rates.
template<typename GenericWellModel, typename VfpProperties, typename Communication = Parallel::Communication>
class NetworkPressureComputation
{
public:
    NetworkPressureComputation(const GenericWellModel& well_model,
                               const Network::ExtNetwork& network,
                               const VfpProperties& vfp_props,
                               const UnitSystem& unit_system,
                               const int report_step_idx,
                               const Communication& comm)
        : well_model_(well_model)
        , network_(network)
        , vfp_props_(vfp_props)
        , unit_system_(unit_system)
        , report_step_idx_(report_step_idx)
        , comm_(comm)
    {
    }

    using Scalar = typename GenericWellModel::Scalar;
    using IndexTraits = GenericWellModel::IndexTraits;
    std::pair<std::map<std::string, Scalar>, std::map<std::string, data::BranchData>> run()
    {
        const auto roots = network_.roots();
        for (const auto& root : roots) {
            // Fixed pressure nodes of the network are the roots of trees.
            // Leaf nodes must correspond to groups in the group structure.
            // Let us first find all leaf nodes of the network. We also
            // create a vector of all nodes, ordered so that a child is
            // always after its parent.
            const auto [root_to_child_nodes, leaf_nodes] = collectTreeNodes(root.get().name());

            // Starting with the leaf nodes of the network, get the flow rates
            // from the corresponding groups.
            auto node_inflows = initializeLeafInflows(leaf_nodes);

            // Accumulate flow rates in the network, towards the roots.
            // Note that a root (i.e. fixed pressure node) can still be
            // contributing flow towards other nodes in the network, i.e.
            // a node can be the root of a subtree.
            accumulateInflows(root_to_child_nodes, node_inflows);

            // Going the other way (from roots to leafs), calculate the pressure
            // at each node using VFP tables and rates.
            computeNodePressures(root_to_child_nodes, node_inflows);

#ifdef OPM_NETWORK_PRESSURE_TRACE
            // Off unless the macro is defined: this builds a string per node per
            // sub-iteration per domain, whether or not the log discards it.
            OpmLog::debug("Network pressure computation completed for root " + root.get().name() + ". Node pressures:");
            for (const auto& [node, pressure] : node_pressures_) {
                OpmLog::debug("Network node " + node + " pressure: " + std::to_string(pressure/1e5) + " bar");
            }
            OpmLog::debug("Node inflows:");
            for (const auto& [node, inflows] : node_inflows) {
                OpmLog::debug("Network node " + node + " inflows: "
                    + std::to_string(inflows[0]*86400) + ", " + std::to_string(inflows[1]*86400) + ", " + std::to_string(inflows[2]*86400));
            }
#endif

        }

        return {node_pressures_, branch_data_};
    }

    /// Nodes whose pressure could not be computed from the VFP tables (and their
    /// descendants). Their entries in the pressure map are placeholders (the upstream
    /// pressure) and must not be used as node pressures.
    const std::set<std::string>& invalidNodes() const
    {
        return invalid_nodes_;
    }

private:
    std::pair<std::vector<std::string>, std::set<std::string>>
    collectTreeNodes(const std::string& root) const
    {
        std::stack<std::string> children;
        std::set<std::string> leaf_nodes;
        std::vector<std::string> root_to_child_nodes;
        children.push(root);
        while (!children.empty()) {
            const auto node = children.top();
            children.pop();
            root_to_child_nodes.push_back(node);
            auto branches = network_.downtree_branches(node);
            if (branches.empty()) {
                leaf_nodes.insert(node);
            }
            for (const auto& branch : branches) {
                children.push(branch.downtree_node());
            }
        }

        assert(children.empty());
        return {root_to_child_nodes, leaf_nodes};
    }

    std::map<std::string, std::vector<Scalar>>
    initializeLeafInflows(const std::set<std::string>& leaf_nodes) const
    {
        std::map<std::string, std::vector<Scalar>> node_inflows;
        const std::vector<Scalar> zero_rates(3, 0.0);

        for (const auto& node : leaf_nodes) {
            // Guard against empty leaf nodes (may not be present in GRUPTREE).
            // Use the domain-correct check so injection networks query the injection
            // rate map rather than the production rate map (which is always empty for
            // pure injection groups, causing zero-rate pressure calculations).
            using Calc = NetworkVfpPressureCalculator<Scalar, IndexTraits, VfpProperties>;
            if (!Calc::hasLeafNodeRate(well_model_.groupStateHelper().groupState(), node)) {
                node_inflows[node] = zero_rates;
                continue;
            }

            node_inflows[node] = Calc::leafNodeRate(well_model_.groupStateHelper().groupState(),
                                                    node);
            if (network_.node(node).add_gas_lift_gas()) {
                addGasLiftGas(node, node_inflows[node]);
            }
        }

        return node_inflows;
    }

    void addGasLiftGas(const std::string& node,
                       std::vector<Scalar>& rates) const
    {
        const auto& group = well_model_.schedule().getGroup(node, report_step_idx_);
        const auto& well_state = well_model_.groupStateHelper().wellState();
        Scalar alq = 0.0;
        // Add gas lift from all wells on this process
        for (const std::string& wellname : group.wells()) {
            const Well& well = well_model_.schedule().getWell(wellname, report_step_idx_);
            if (well.isInjector() || !well_state.isOpen(wellname)) {
                continue;
            }

            const Scalar efficiency = well.getEfficiencyFactor(/*network*/ true)
                * well_state.getGlobalEfficiencyScalingFactor(wellname);
            const auto& well_index = well_state.index(wellname);
            if (well_index.has_value() &&
                well_state.wellIsOwned(well_index.value(), wellname))
            {
                alq += well_state.well(wellname).alq_state.get() * efficiency;
            }
        }
        // Sum ALQ across all processes to get total ALQ for the node.
        // Note that communication is required here since each
        // process has different wells, and the loop above therefore
        // only considers local wells.
        // However, all processes have all groups and their rates available,
        // so we do not need to communicate those.
        alq = comm_.sum(alq);
        // Only add satellite production once for parallel runs
        // (i.e. add after communication)
        if (group.hasSatelliteProduction()) {
            const auto& gsat_prod = well_model_.schedule()[report_step_idx_].gsatprod().get(node, well_model_.summaryState());
            alq += gsat_prod.rate[GSatProd::GSatProdGroupProp::Rate::GLift];
        }

        rates[IndexTraits::gasPhaseIdx] += alq;
    }

    void accumulateInflows(const std::vector<std::string>& root_to_child_nodes,
                           std::map<std::string, std::vector<Scalar>>& node_inflows) const
    {
        const auto child_to_root_nodes = std::ranges::reverse_view(root_to_child_nodes);

        for (const auto& node : child_to_root_nodes) {
            const auto upbranch = network_.uptree_branch(node);
            if (!upbranch) {
                continue;
            }
            std::vector<Scalar>& up = node_inflows[(*upbranch).uptree_node()];
            const std::vector<Scalar>& down = node_inflows[node];
            // NEFAC support
            const Scalar efficiency = network_.node(node).efficiency();
            if (up.empty()) {
                up = std::vector<Scalar>(down.size(), 0.0);
            }
            assert(up.size() == down.size());
            for (std::size_t ii = 0; ii < up.size(); ++ii) {
                up[ii] += efficiency * down[ii];
            }
        }
    }

    void computeNodePressures(const std::vector<std::string>& root_to_child_nodes,
                              const std::map<std::string, std::vector<Scalar>>& node_inflows)
    {
        for (const auto& node : root_to_child_nodes) {
            // Do not traverse subtree more than once
            if (node_pressures_.find(node) != node_pressures_.end()) {
                continue;
            }

            const auto terminal_pressure = network_.node(node).terminal_pressure();
            const auto upbranch = network_.uptree_branch(node);
            assert(upbranch || terminal_pressure); // If not root, must have uptree branch, and if root, must have terminal pressure.
            using Calc = NetworkVfpPressureCalculator<Scalar, IndexTraits, VfpProperties>;

            if (terminal_pressure) {
                node_pressures_[node] = *terminal_pressure;
                if (upbranch) {
                    // If terminal pressure is specified on a non-root node, we still want to calculate the branch data for the uptree branch.
                    const Scalar up_press = node_pressures_[(*upbranch).uptree_node()];
                    auto rates = node_inflows.at(node);
                    branch_data_.try_emplace(node,
                                             *terminal_pressure - up_press,
                                             rates[IndexTraits::oilPhaseIdx],
                                             rates[IndexTraits::waterPhaseIdx],
                                             rates[IndexTraits::gasPhaseIdx]);
                } else {
                    // Root node with terminal pressure and no uptree branch, inserting a zero-valued placeholder.
                    branch_data_.emplace(node, data::BranchData{0.0, 0.0, 0.0, 0.0});
                }
                continue;
            }

            const std::string& up_node = (*upbranch).uptree_node();
            const Scalar up_press = node_pressures_[up_node];
            // Descendants of a node without a valid pressure have none either.
            if (invalid_nodes_.count(up_node) > 0) {
                invalid_nodes_.insert(node);
            }
            const auto vfp_table = (*upbranch).vfp_table();
            if (!vfp_table) {
                // Table number specified as 9999 in the deck, no pressure loss.
                if (network_.node(node).as_choke()) {
                    // Node pressure is set to the group THP.
                    node_pressures_[node] = well_model_.groupStateHelper().groupState().well_group_thp(node);
                } else {
                    node_pressures_[node] = up_press;
                }
                auto rates = node_inflows.at(node);
                branch_data_.try_emplace(node,
                                         node_pressures_[node] - up_press,
                                         rates[IndexTraits::oilPhaseIdx],
                                         rates[IndexTraits::waterPhaseIdx],
                                         rates[IndexTraits::gasPhaseIdx]);
                continue;
            }

            OPM_TIMEBLOCK(NetworkVfpCalculations);
            auto rates = node_inflows.at(node);
            assert(rates.size() == 3);
            Calc::prepareRates(rates);
            const auto branch = Calc::compute(vfp_props_, *vfp_table, rates, up_press, *upbranch, unit_system_);
            // An invalid lookup (zero-filled table cells) gets the upstream pressure as a
            // placeholder so downstream lookups stay in range; callers must consult invalidNodes().
            const Scalar node_pressure = branch.valid ? branch.pressure : up_press;
            if (!branch.valid) {
                invalid_nodes_.insert(node);
            }
            if (!branch.valid || branch.clamped) {
                OpmLog::debug(fmt::format("Network branch {} -> {}: VFP table {} {} at rates ({:.4g}, {:.4g}, {:.4g}) sm3/d, "
                                          "upstream pressure {:.2f} bar",
                                          up_node, node, *vfp_table,
                                          branch.valid ? "lookup clamped to the table axes" : "has no solution",
                                          rates[IndexTraits::waterPhaseIdx] * unit::day,
                                          rates[IndexTraits::oilPhaseIdx] * unit::day,
                                          rates[IndexTraits::gasPhaseIdx] * unit::day,
                                          up_press / unit::barsa));
            }
            node_pressures_[node] = node_pressure;
            // Prefer inserting after computing the pressure, hence negating rates
            branch_data_.try_emplace(node,
                                     node_pressure - up_press,
                                     -rates[IndexTraits::oilPhaseIdx],
                                     -rates[IndexTraits::waterPhaseIdx],
                                     -rates[IndexTraits::gasPhaseIdx]);
        }
    }

    const GenericWellModel& well_model_;
    const Network::ExtNetwork& network_;
    const VfpProperties& vfp_props_;
    const UnitSystem& unit_system_;
    const int report_step_idx_;
    const Communication& comm_;
    std::map<std::string, Scalar> node_pressures_;
    std::map<std::string, data::BranchData> branch_data_;
    std::set<std::string> invalid_nodes_;
};

} // namespace Opm

#endif // OPM_BLACKOIL_WELL_MODEL_NETWORK_PRESSURE_COMPUTATION_HPP
