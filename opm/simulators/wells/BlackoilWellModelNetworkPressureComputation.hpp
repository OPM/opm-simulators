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

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>

#include <algorithm>
#include <cassert>
#include <map>
#include <ranges>
#include <set>
#include <stack>
#include <string>
#include <vector>

namespace Opm {

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
    static const std::vector<Scalar>
    leafNodeRate(const GroupState& group_state, const std::string& node)
    {
        return group_state.network_leaf_node_production_rates(node);
    }

    template<typename Branch>
    static Scalar compute(const VFPProdProperties<Scalar>& vfp_props,
                          const int table_id,
                          const std::vector<Scalar>& rates,
                          const Scalar up_press,
                          const Branch& upbranch,
                          const UnitSystem& unit_system)
    {
        // NB! ALQ in extended network is never implicitly the gas lift rate (GRAT), i.e., the
        //     gas lift rates only enters the network pressure calculations through the rates
        //     (e.g., in GOR calculations) unless a branch ALQ is set in BRANPROP.
        const auto alq_type = vfp_props.getTable(table_id).getALQType();
        const auto dimension = VFPProdTable::ALQDimension(alq_type, unit_system);
        const Scalar alq = upbranch.alq_value(dimension).value_or(0.0);

        return vfp_props.bhp(table_id,
                             rates[IndexTraits::waterPhaseIdx],
                             rates[IndexTraits::oilPhaseIdx],
                             rates[IndexTraits::gasPhaseIdx],
                             up_press,
                             alq,
                             0.0, // explicit_wfr
                             0.0, // explicit_gfr
                             false); // use_expvfp we dont support explicit lookup
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
    static const std::vector<Scalar>
    leafNodeRate(const GroupState& group_state, const std::string& node)
    {
        return group_state.network_leaf_node_injection_rates(node);
    }

    template<typename Branch>
    static Scalar compute(const VFPInjProperties<Scalar>& vfp_props,
                          const int table_id,
                          const std::vector<Scalar>& rates,
                          const Scalar up_press,
                          const Branch&,
                          const UnitSystem&)
    {
        return vfp_props.bhp(table_id,
                             rates[IndexTraits::waterPhaseIdx],
                             rates[IndexTraits::oilPhaseIdx],
                             rates[IndexTraits::gasPhaseIdx],
                             up_press);
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

    std::map<std::string, Scalar> run()
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
        }

        return node_pressures_;
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
            // Guard against empty leaf nodes (may not be present in GRUPTREE)
            if (!well_model_.groupStateHelper().groupState().has_production_rates(node)) {
                node_inflows[node] = zero_rates;
                continue;
            }

            using Calc = NetworkVfpPressureCalculator<Scalar, IndexTraits, VfpProperties>;
            node_inflows[node] = Calc::leafNodeRate(well_model_.groupStateHelper().groupState(), node);
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
            const auto terminal_pressure = network_.node(node).terminal_pressure();
            if (terminal_pressure) {
                node_pressures_[node] = *terminal_pressure;
                continue;
            }

            const auto upbranch = network_.uptree_branch(node);
            assert(upbranch);
            const Scalar up_press = node_pressures_[(*upbranch).uptree_node()];
            const auto vfp_table = (*upbranch).vfp_table();
            if (!vfp_table) {
                // Table number specified as 9999 in the deck, no pressure loss.
                if (network_.node(node).as_choke()) {
                    // Node pressure is set to the group THP.
                    node_pressures_[node] = well_model_.groupStateHelper().groupState().well_group_thp(node);
                } else {
                    node_pressures_[node] = up_press;
                }
                continue;
            }

            OPM_TIMEBLOCK(NetworkVfpCalculations);
            auto rates = node_inflows.at(node);
            assert(rates.size() == 3);
            using Calc = NetworkVfpPressureCalculator<Scalar, IndexTraits, VfpProperties>;
            Calc::prepareRates(rates);
            node_pressures_[node]
                = Calc::compute(vfp_props_, *vfp_table, rates, up_press, *upbranch, unit_system_);
        }
    }

    const GenericWellModel& well_model_;
    const Network::ExtNetwork& network_;
    const VfpProperties& vfp_props_;
    const UnitSystem& unit_system_;
    const int report_step_idx_;
    const Communication& comm_;
    std::map<std::string, Scalar> node_pressures_;
};

} // namespace Opm

#endif // OPM_BLACKOIL_WELL_MODEL_NETWORK_PRESSURE_COMPUTATION_HPP
