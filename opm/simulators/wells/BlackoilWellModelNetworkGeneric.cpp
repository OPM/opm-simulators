/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#include <opm/simulators/wells/BlackoilWellModelNetworkGeneric.hpp>

#include <opm/common/TimingMacros.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/wells/BlackoilWellModelNetworkPressureComputation.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>

namespace Opm {

namespace details {
    /// Helper to check if any network (production, gas injection, water injection) is active at a given time step.
    bool anyNetworkActive(const Schedule& schedule, const int timeStepIdx)
    {
        const auto& sstate = schedule[timeStepIdx];
        return sstate.network().active()
            || (sstate.injectionNetwork.get_ptr(Phase::GAS) != nullptr
                && sstate.injectionNetwork.get_ptr(Phase::GAS)->active())
            || (sstate.injectionNetwork.get_ptr(Phase::WATER) != nullptr
                && sstate.injectionNetwork.get_ptr(Phase::WATER)->active());
    }

    /// Helper to get all active networks (production, gas injection, water injection) at a given time step.
    std::vector<ActiveNetworkDescriptor>
    activeNetworks(const Schedule& schedule, const int timeStepIdx)
    {
        std::vector<ActiveNetworkDescriptor> active_networks;
        const auto& sstate = schedule[timeStepIdx];
        if (sstate.network().active()) {
            active_networks.push_back({NetworkDomain::Production, std::cref(sstate.network())});
        }
        if (sstate.injectionNetwork.get_ptr(Phase::GAS) != nullptr
            && sstate.injectionNetwork.get_ptr(Phase::GAS)->active()) {
            active_networks.push_back({NetworkDomain::InjectionGas, std::cref(*sstate.injectionNetwork.get_ptr(Phase::GAS))});
        }
        if (sstate.injectionNetwork.get_ptr(Phase::WATER) != nullptr
            && sstate.injectionNetwork.get_ptr(Phase::WATER)->active()) {
            active_networks.push_back({NetworkDomain::InjectionWater, std::cref(*sstate.injectionNetwork.get_ptr(Phase::WATER))});
        }
        return active_networks;
    }

    std::optional<Phase> injectionPhaseForDomain(const NetworkDomain domain)
    {
        switch (domain) {
        case NetworkDomain::InjectionGas:
            return Phase::GAS;
        case NetworkDomain::InjectionWater:
            return Phase::WATER;
        case NetworkDomain::Production:
        case NetworkDomain::Count:
            return std::nullopt;
        }

        return std::nullopt;
    }
} // namespace details


template<typename Scalar, typename IndexTraits>
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
BlackoilWellModelNetworkGeneric(BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model)
    : well_model_(well_model)
{
    this->setFromRestart(well_model_.eclState().getRestartNetworkPressures());
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
setFromRestart(const std::optional<std::map<std::string, double>>& node_pressures)
{
    if (node_pressures.has_value()) {
        if constexpr (std::is_same_v<Scalar,double>) {
            this->node_pressures_ = node_pressures.value();
        } else {
            for (const auto& it : node_pressures.value()) {
                this->node_pressures_[it.first] = it.second;
            }
        }
        this->syncProductionDomainState_();
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateActiveState(const int report_step)
{
    this->active_ = false;
    for (const auto& network : details::activeNetworks(well_model_.schedule(), report_step)) {
        updateActiveStateImpl(network.network.get());
    }
    this->active_ = well_model_.comm().max(active_);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateActiveStateImpl(const Network::ExtNetwork& network)
{
    // Accumulates into active_ across the domains; an inactive network must not
    // clear what an earlier domain set.
    if (!network.active()) {
        return;
    }
    bool network_active = false;
    for (const auto& well : well_model_.genericWells()) {
        const bool is_partof_network = network.has_node(well->wellEcl().groupName());
        const bool prediction_mode = well->wellEcl().predictionMode();
        if (is_partof_network && prediction_mode) {
            network_active = true;
            break;
        }
    }
#ifdef RESERVOIR_COUPLING_ENABLED
    // A reservoir coupling master may have no local wells in the network
    // leaf groups, the leaf rates come from slave-reported
    // network_surface_rates instead.  Without this clause the master's
    // network solver short-circuits via active_=false and the leaf
    // pressures are never iterated.
    const auto& rescoup_proxy = well_model_.groupStateHelper().rescoup();
    if (!network_active && rescoup_proxy.isMaster()) {
        const auto& rescoup_master = rescoup_proxy.master();
        const auto num_slaves = rescoup_master.numSlaves();
        for (std::size_t s = 0; s < num_slaves && !network_active; ++s) {
            // Only an activated slave supplies leaf rates this step; a master
            // group whose slave is inactive must not keep the network active
            // (it would be solved against missing/stale slave rates). This
            // matches the slaveIsActivated gate used in the coupled-network
            // iteration.
            if (!rescoup_master.slaveIsActivated(s)) {
                continue;
            }
            for (const auto& master_group : rescoup_master.getMasterGroupNamesForSlave(s)) {
                if (network.has_node(master_group)) {
                    network_active = true;
                    break;
                }
            }
        }
    }
#endif
    this->active_ = this->active_ || network_active;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
needPreStepRebalance(const int report_step) const
{
    const auto active_networks = details::activeNetworks(well_model_.schedule(), report_step);
    bool network_rebalance_necessary = false;
    for (const auto& well : well_model_.genericWells()) {
        const bool is_partof_network = std::any_of(active_networks.begin(),
                                                   active_networks.end(),
                                                   [&](const auto& network)
                                                   {
                                                       return network.network.get().has_node(well->wellEcl().groupName());
                                                   });
        // TODO: we might find more relevant events to be included here (including network change events?)
        const auto& events = well_model_.wellState().well(well->indexOfWell()).events;
        if (is_partof_network && events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE)) {
            network_rebalance_necessary = true;
            break;
        }
    }
    network_rebalance_necessary = well_model_.comm().max(network_rebalance_necessary);
    return network_rebalance_necessary;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
shouldBalance(const int reportStepIdx) const
{
    // if network is not active, we do not need to balance the network
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return false;
    }

    const auto& balance = well_model_.schedule()[reportStepIdx].network_balance();
    const auto& iterCtx = well_model_.iterationContext();
    if (balance.mode() == Network::Balance::CalcMode::TimeStepStart) {
        return iterCtx.isFirstGlobalIteration();
    } else if (balance.mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = well_model_.schedule()[reportStepIdx].nupcol();
        return iterCtx.withinNupcol(nupcol);
    } else {
        // We do not support any other rebalancing modes,
        // i.e. TimeInterval based rebalancing is not available.
        // This should be warned about elsewhere, so we choose to
        // avoid spamming with a warning here.
        return false;
    }
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
willBalanceOnNextIteration(const int reportStepIdx) const
{
    // if network is not active, we do not need to balance the network
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return false;
    }

    const auto& schedule_state = well_model_.schedule()[reportStepIdx];
    if (schedule_state.network_balance().mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = schedule_state.nupcol();
        return well_model_.iterationContext().withinNupcol(nupcol - 1); // Note the -1 here!
    } else {
        // Any other rebalancing mode will only rebalance
        // at the start of the timestep.
        return false;
    }
}

template<typename Scalar, typename IndexTraits>
Scalar
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updatePressures(const int reportStepIdx,
                const Scalar damping_factor,
                const Scalar upper_update_bound)
{
    OPM_TIMEFUNCTION();
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return 0.0;
    }

    this->syncProductionDomainState_();
    const auto previous_node_pressures = this->domain_node_pressures_;

    for (const auto& network : details::activeNetworks(well_model_.schedule(), reportStepIdx)) {
        NetworkPressures result;
        if (network.domain == details::NetworkDomain::Production) {
            result = this->computePressures(network.network.get(),
                                            *well_model_.getVFPProperties().getProd(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
        } else {
            const auto injection_phase = details::injectionPhaseForDomain(network.domain);
            assert(injection_phase.has_value());
            result = this->computePressures(network.network.get(),
                                            *well_model_.getVFPProperties().getInj(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm(),
                                            *injection_phase);
        }
        this->nodePressures(network.domain) = std::move(result.node_pressures);
        this->branchData(network.domain) = std::move(result.branch_data);
        this->invalidNodes(network.domain) = std::move(result.invalid_nodes);
    }
    this->syncLegacyProductionState_();

    // here, the network imbalance is the difference between the previous nodal pressure and the new nodal pressure
    Scalar network_imbalance = 0.;
    if (!this->active()) {
        return network_imbalance;
    }

    for (const auto& network : details::activeNetworks(well_model_.schedule(), reportStepIdx)) {
        auto& domain_pressures = this->nodePressures(network.domain);
        const auto& invalid = this->invalidNodes(network.domain);
        const auto& previous_domain_pressures = previous_node_pressures[details::domainIndex(network.domain)];

        if (!invalid.empty()) {
            // The VFP tables gave no pressure for these nodes (rate/pressure outside what
            // the tables can deliver). Keep the previous value and report the network as
            // unbalanced so that the wells get another chance to move into range.
            network_imbalance = std::max(network_imbalance, upper_update_bound);
            if (this->invalid_nodes_report_step_ != reportStepIdx) {
                this->invalid_nodes_report_step_ = reportStepIdx;
                OpmLog::warning(fmt::format("Network: no VFP solution for node(s) {} at report step {}; "
                                            "keeping the previous node pressure(s).",
                                            fmt::join(invalid, ", "), reportStepIdx + 1));
            }
        }

        if (!previous_domain_pressures.empty()) {
            for (auto& [name, new_pressure]: domain_pressures) {
                if (previous_domain_pressures.count(name) <= 0) {
                    if (std::abs(new_pressure) > network_imbalance) {
                        network_imbalance = std::abs(new_pressure);
                    }
                    continue;
                }

                const auto pressure = previous_domain_pressures.at(name);
                if (invalid.count(name) > 0) {
                    new_pressure = pressure;
                    continue;
                }
                const Scalar change = (new_pressure - pressure);
                if (std::abs(change) > network_imbalance) {
                    network_imbalance = std::abs(change);
                }
                // We dampen the nodal pressure change during one iteration since our nodal pressure calculation
                // is somewhat explicit. There is a relative dampening factor applied to the update value, and also
                // the maximum update is limited (to 5 bar by default, can be changed with --network-max-pressure-update-in-bars).
                const Scalar damped_change = std::min(damping_factor * std::abs(change), upper_update_bound);
                const Scalar sign = change > 0 ? 1. : -1.;
                new_pressure = pressure + sign * damped_change;
            }
            continue;
        }

        for (const auto& [name, pressure]: domain_pressures) {
            if (std::abs(pressure) > network_imbalance) {
                network_imbalance = std::abs(pressure);
            }
        }
    }
    this->syncLegacyProductionState_();

    for (auto& well : well_model_.genericWells()) {

        if (!well->wellEcl().predictionMode()) {
            continue;
        }

        const auto domain = details::domainForWell(*well);

        if (!domain.has_value()) {
            continue;
        }

        const auto it = this->nodePressures(*domain).find(well->wellEcl().groupName());
        if (it != this->nodePressures(*domain).end()) {
            if (this->invalidNodes(*domain).count(well->wellEcl().groupName()) > 0) {
                // No valid leaf pressure this iteration; keep the well's current THP limit.
                continue;
            }
            if (well->isProducer()) {
                // For producers the leaf-node pressure is the group wellhead THP;
                // apply it directly as a dynamic THP constraint.
                const Scalar new_limit = it->second;
                well->setDynamicThpLimit(new_limit);
                SingleWellState<Scalar, IndexTraits>& ws = well_model_.wellState()[well->indexOfWell()];
                const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                // TODO: not sure why the thp is NOT updated properly elsewhere
                if (thp_is_limit) {
                    ws.thp = well->getTHPConstraint(well_model_.summaryState());
                }
            } else if (well->isInjector() && well->wellEcl().vfp_table_number() > 0) {
                // For injectors, apply the network leaf-node pressure as a dynamic THP only
                // if it falls within the individual well's VFPINJ table THP range.
                // If P_leaf is outside the range, computeBhpAtThpLimitInj would extrapolate
                // to invalid values and mark the well inoperable, causing a rate collapse.
                const auto& inj_vfp = *well_model_.getVFPProperties().getInj();
                const int table_id = well->wellEcl().injectionControls(
                    well_model_.summaryState()).vfp_table_number;
                if (inj_vfp.hasTable(table_id)) {
                    const auto& thp_axis = inj_vfp.getTable(table_id).getTHPAxis();
                    const Scalar min_thp = static_cast<Scalar>(thp_axis.front());
                    const Scalar max_thp = static_cast<Scalar>(thp_axis.back());
                    const Scalar new_limit = it->second;
                    if (new_limit >= min_thp && new_limit <= max_thp) {
                        well->setDynamicThpLimit(new_limit);
                        SingleWellState<Scalar, IndexTraits>& ws =
                            well_model_.wellState()[well->indexOfWell()];
                        const bool thp_is_limit =
                            ws.injection_cmode == Well::InjectorCMode::THP;
                        if (thp_is_limit) {
                            ws.thp = well->getTHPConstraint(well_model_.summaryState());
                        }
                    }
                }
            }
        }
    }
    return network_imbalance;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
assignNodeAndBranchValues(std::map<std::string, data::NodeData>& nodevalues,
                          std::map<std::string, data::BranchData>& branchvalues,
                          std::map<std::string, data::BranchData>& converged_branchvalues,
                          const int reportStepIdx) const
{
    nodevalues.clear();
    branchvalues.clear();
    converged_branchvalues.clear();
    if (reportStepIdx < 0) return;
    for (const auto& [node, pressure] : node_pressures_) {
        nodevalues.emplace(node, data::NodeData{pressure});
        // Assign node values of well groups to GPR:WELLNAME
        const auto& sched = well_model_.schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            nodevalues.emplace(wellname, data::NodeData{pressure});
        }
    }
    for (const auto& [branch, branch_data] : branch_data_) {
        branchvalues.emplace(branch, branch_data);
        // Skip wells (do not consider well->group a branch, at least not for now)
    }

    const auto& network = well_model_.schedule()[reportStepIdx].network();
    if (!network.active()) {
        return;
    }

    auto converged = this->computePressures(network,
                                            *well_model_.getVFPProperties().getProd(),
                                            well_model_.schedule().getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
    const auto& converged_pressures = converged.node_pressures;
    converged_branchvalues = std::move(converged.branch_data);
    for (const auto& [node, converged_pressure] : converged_pressures) {
        auto it = nodevalues.find(node);
        assert(it != nodevalues.end() );
        it->second.converged_pressure = converged_pressure;
        // Assign node values of group to GPR:WELLNAME
        const auto& sched = well_model_.schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            auto it2 = nodevalues.find(wellname);
            assert(it2 != nodevalues.end());
            it2->second.converged_pressure = converged_pressure;
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
initialize(const int report_step)
{
    if (details::anyNetworkActive(well_model_.schedule(), report_step)) {
        this->syncProductionDomainState_();
        for (auto& well : well_model_.genericWells()) {
            initializeWell(*well);
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
initializeWell(WellInterfaceGeneric<Scalar,IndexTraits>& well)
{
    const auto domain = details::domainForWell(well);

    if (domain.has_value() && !this->nodePressures(*domain).empty()) {
        const auto it = this->nodePressures(*domain).find(well.wellEcl().groupName());
        if (it != this->nodePressures(*domain).end()) {
            // For producers, carry forward the previous step's converged network pressure
            // so that prepareTimeStep() starts with the correct THP constraint.
            // For injectors, do NOT set dynamic_thp_limit_ here.  Setting it before
            // prepareTimeStep() causes solveWellEquation() to switch the injector to THP
            // mode; the resulting rate change propagates through the Newton loop and
            // produces large network imbalances that fail to converge in the allowed
            // iterations.  The injection network THP is applied for the first time during
            // the Newton loop via updatePressures(), where the stale-potential bypass in
            // WellConstraints::activeInjectionConstraint ensures correct switching.
            if (well.isProducer()) {
                well.setDynamicThpLimit(it->second);
            }
        }
    }
}

template <typename Scalar, typename IndexTraits>
typename BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::NetworkPressures
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
computePressures(const Network::ExtNetwork& network,
                 const VFPProdProperties<Scalar>& vfp_prod_props,
                 const UnitSystem& unit_system,
                 const int reportStepIdx,
                 const Parallel::Communication& comm) const
{
    OPM_TIMEFUNCTION();
    if (!network.active()) {
        return {};
    }

    NetworkPressureComputation<BlackoilWellModelGeneric<Scalar, IndexTraits>,
                               VFPProdProperties<Scalar>>
        network_pressure_computation(
            well_model_, network, vfp_prod_props, unit_system, reportStepIdx, comm);

    auto [node_pressures, branch_data] = network_pressure_computation.run();
    return {std::move(node_pressures), std::move(branch_data), network_pressure_computation.invalidNodes()};
}

template <typename Scalar, typename IndexTraits>
typename BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::NetworkPressures
BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
computePressures(const Network::ExtNetwork& network,
                 const VFPInjProperties<Scalar>& vfp_inj_props,
                 const UnitSystem& unit_system,
                 const int reportStepIdx,
                 const Parallel::Communication& comm,
                 const Phase injectionPhase) const
{
    OPM_TIMEFUNCTION();
    if (!network.active()) {
        return {};
    }

    NetworkPressureComputation<BlackoilWellModelGeneric<Scalar, IndexTraits>,
                               VFPInjProperties<Scalar>>
        network_pressure_computation(
            well_model_, network, vfp_inj_props, unit_system, reportStepIdx, comm, injectionPhase);

    auto [node_pressures, branch_data] = network_pressure_computation.run();
    return {std::move(node_pressures), std::move(branch_data), network_pressure_computation.invalidNodes()};
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const
{
    return
           this->active_ == rhs.active_
        && this->node_pressures_ == rhs.node_pressures_
        && this->last_valid_node_pressures_ == rhs.last_valid_node_pressures_
        && this->branch_data_ == rhs.branch_data_
        && this->last_valid_branch_data_ == rhs.last_valid_branch_data_
        && this->domain_node_pressures_ == rhs.domain_node_pressures_
        && this->last_valid_domain_node_pressures_ == rhs.last_valid_domain_node_pressures_
        && this->domain_branch_data_ == rhs.domain_branch_data_
        && this->last_valid_domain_branch_data_ == rhs.last_valid_domain_branch_data_;
}

template class BlackoilWellModelNetworkGeneric<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelNetworkGeneric<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
