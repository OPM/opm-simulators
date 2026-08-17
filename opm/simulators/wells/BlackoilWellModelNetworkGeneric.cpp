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
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
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

    template<class Well>
    std::optional<NetworkDomain> domainForWell(const Well& well)
    {
        if (well.isProducer()) {
            return NetworkDomain::Production;
        }
        if (well.isInjector()) {
            if (well.wellEcl().injectorType() == InjectorType::GAS) {
                return NetworkDomain::InjectionGas;
            }
            if (well.wellEcl().injectorType() == InjectorType::WATER) {
                return NetworkDomain::InjectionWater;
            }
        }
        return std::nullopt;
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
                const Scalar upper_update_bound,
                const bool use_secant)
{
    OPM_TIMEFUNCTION();
    if (!details::anyNetworkActive(well_model_.schedule(), reportStepIdx)) {
        return 0.0;
    }

    this->syncProductionDomainState_();
    const auto previous_node_pressures = this->domain_node_pressures_;

    // Per domain and leaf: the lowest wellhead pressure of the wells that are open,
    // flowing and not on THP control. The leaf rate does not depend on the node
    // pressure above it (see NodePressureUpdater::next).
    std::array<std::map<std::string, Scalar>, details::domainIndex(details::NetworkDomain::Count)> plateau_floor;
    for (const auto& well : well_model_.genericWells()) {
        const auto domain = details::domainForWell(*well);
        if (!domain.has_value() || !well->wellEcl().predictionMode()) {
            continue;
        }
        const auto& ws = well_model_.wellState()[well->indexOfWell()];
        const bool on_thp = well->isProducer() ? ws.production_cmode == Well::ProducerCMode::THP
                                               : ws.injection_cmode == Well::InjectorCMode::THP;
        const bool flowing = ws.status == WellStatus::OPEN
            && std::any_of(ws.surface_rates.begin(), ws.surface_rates.end(),
                           [](const Scalar q) { return q != Scalar{0}; });
        if (on_thp || !flowing || ws.thp <= 0.0) {
            continue;
        }
        auto& floors = plateau_floor[details::domainIndex(*domain)];
        auto [it, inserted] = floors.try_emplace(well->wellEcl().groupName(), ws.thp);
        if (!inserted) {
            it->second = std::min(it->second, ws.thp);
        }
    }

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
            // the tables can deliver).
            if (this->invalid_nodes_report_step_ != reportStepIdx) {
                this->invalid_nodes_report_step_ = reportStepIdx;
                OpmLog::warning(fmt::format("Network: no VFP solution for node(s) {} at report step {}; "
                                            "treating them as too high.",
                                            fmt::join(invalid, ", "), reportStepIdx + 1));
            }
        }

        if (!previous_domain_pressures.empty()) {
            auto& updaters = this->pressure_updaters_[details::domainIndex(network.domain)];
            for (auto& [name, computed_pressure]: domain_pressures) {
                if (previous_domain_pressures.count(name) <= 0) {
                    if (std::abs(computed_pressure) > network_imbalance) {
                        network_imbalance = std::abs(computed_pressure);
                    }
                    continue;
                }

                // pressure is what the wells were last solved with, computed_pressure what
                // the network gives for the resulting rates.
                const auto pressure = previous_domain_pressures.at(name);
                const bool valid = invalid.count(name) == 0;
                if (use_secant) {
                    auto& updater = updaters[name];
                    const auto& floors = plateau_floor[details::domainIndex(network.domain)];
                    std::optional<Scalar> floor;
                    if (const auto f = floors.find(name); f != floors.end()) {
                        floor = f->second;
                    }
                    computed_pressure = updater.next(pressure, computed_pressure, valid,
                                                     damping_factor, upper_update_bound, floor);
                    // The residual is amplified by the well response; judge convergence on
                    // the remaining pressure uncertainty instead.
                    network_imbalance = std::max(network_imbalance, updater.error());
                } else if (!valid) {
                    // Keep the previous value; report as unbalanced so the wells get another
                    // chance to move into range.
                    network_imbalance = std::max(network_imbalance, upper_update_bound);
                    computed_pressure = pressure;
                } else {
                    network_imbalance = std::max(network_imbalance, std::abs(computed_pressure - pressure));
                    // We dampen the nodal pressure change during one iteration since our nodal pressure calculation
                    // is somewhat explicit. There is a relative dampening factor applied to the update value, and also
                    // the maximum update is limited (to 5 bar by default, can be changed with --network-max-pressure-update-in-bars).
                    computed_pressure = NodePressureUpdater<Scalar>::damped(pressure, computed_pressure - pressure,
                                                                            damping_factor, upper_update_bound);
                }
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

        std::optional<details::NetworkDomain> domain;
        if (well->isProducer()) {
            domain = details::NetworkDomain::Production;
        } else if (well->isInjector()) {
            if (well->wellEcl().injectorType() == InjectorType::GAS) {
                domain = details::NetworkDomain::InjectionGas;
            } else if (well->wellEcl().injectorType() == InjectorType::WATER) {
                domain = details::NetworkDomain::InjectionWater;
            }
        }

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
                // For injectors the leaf-node pressure is the available wellhead pressure.
                // Clamp it to the well's VFPINJ THP axis: outside the axis the well's own
                // THP->BHP lookup would extrapolate, and a limit at the axis end is the
                // closest statement the table can make (the BHP/rate limits then bind).
                const auto& inj_vfp = *well_model_.getVFPProperties().getInj();
                const int table_id = well->wellEcl().injectionControls(
                    well_model_.summaryState()).vfp_table_number;
                if (inj_vfp.hasTable(table_id)) {
                    const auto& thp_axis = inj_vfp.getTable(table_id).getTHPAxis();
                    const Scalar min_thp = static_cast<Scalar>(thp_axis.front());
                    const Scalar max_thp = static_cast<Scalar>(thp_axis.back());
                    const Scalar new_limit = std::clamp(it->second, min_thp, max_thp);
                    well->setDynamicThpLimit(new_limit);
                    SingleWellState<Scalar, IndexTraits>& ws =
                        well_model_.wellState()[well->indexOfWell()];
                    if (ws.injection_cmode == Well::InjectorCMode::THP) {
                        ws.thp = well->getTHPConstraint(well_model_.summaryState());
                    }
                }
            }
        }
    }
    return network_imbalance;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
assignNodeAndBranchValues(data::GroupAndNetworkValues& values,
                          const int reportStepIdx) const
{
    auto& nodevalues = values.nodeData;
    auto& branchvalues = values.branchData;
    auto& converged_branchvalues = values.convergedBranchData;
    nodevalues.clear();
    branchvalues.clear();
    converged_branchvalues.clear();
    values.gasInjNodeData.clear();
    values.gasInjBranchData.clear();
    values.waterInjNodeData.clear();
    values.waterInjBranchData.clear();
    if (reportStepIdx < 0) return;

    const auto& sched = well_model_.schedule();
    // Node values are also assigned to the wells of the node's group (GPR:WELLNAME).
    auto assign_nodes = [&sched, reportStepIdx](const std::map<std::string, Scalar>& pressures,
                                                std::map<std::string, data::NodeData>& out)
    {
        for (const auto& [node, pressure] : pressures) {
            out.emplace(node, data::NodeData{pressure});
            if (!sched.hasGroup(node, reportStepIdx)) {
                continue;
            }
            for (const std::string& wellname : sched.getGroup(node, reportStepIdx).wells()) {
                out.emplace(wellname, data::NodeData{pressure});
            }
        }
    };

    assign_nodes(node_pressures_, nodevalues);
    for (const auto& [branch, branch_data] : branch_data_) {
        branchvalues.emplace(branch, branch_data);
        // Skip wells (do not consider well->group a branch, at least not for now)
    }

    // Injection networks: current pressures only (no converged variant reported).
    assign_nodes(this->nodePressures(details::NetworkDomain::InjectionGas), values.gasInjNodeData);
    values.gasInjBranchData = this->branchData(details::NetworkDomain::InjectionGas);
    assign_nodes(this->nodePressures(details::NetworkDomain::InjectionWater), values.waterInjNodeData);
    values.waterInjBranchData = this->branchData(details::NetworkDomain::InjectionWater);

    const auto& network = sched[reportStepIdx].network();
    if (!network.active()) {
        return;
    }

    auto converged = this->computePressures(network,
                                            *well_model_.getVFPProperties().getProd(),
                                            sched.getUnits(),
                                            reportStepIdx,
                                            well_model_.comm());
    const auto& converged_pressures = converged.node_pressures;
    converged_branchvalues = std::move(converged.branch_data);
    for (const auto& [node, converged_pressure] : converged_pressures) {
        auto it = nodevalues.find(node);
        assert(it != nodevalues.end() );
        it->second.converged_pressure = converged_pressure;
        // Assign node values of group to GPR:WELLNAME
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
    std::optional<details::NetworkDomain> domain;
    if (well.isProducer()) {
        domain = details::NetworkDomain::Production;
    } else if (well.isInjector()) {
        if (well.wellEcl().injectorType() == InjectorType::GAS) {
            domain = details::NetworkDomain::InjectionGas;
        } else if (well.wellEcl().injectorType() == InjectorType::WATER) {
            domain = details::NetworkDomain::InjectionWater;
        }
    }

    if (domain.has_value() && !this->nodePressures(*domain).empty()) {
        const auto it = this->nodePressures(*domain).find(well.wellEcl().groupName());
        if (it != this->nodePressures(*domain).end()) {
            // Carry forward the previous step's converged network pressure so that
            // prepareTimeStep() starts with the right THP constraint. Without it an
            // injector starts a new report step at its WCONINJE THP, injects at its
            // rate limit, and starves the other leaves of the network before the
            // first network balance.
            Scalar limit = it->second;
            if (well.isInjector()) {
                const auto& inj_vfp = *well_model_.getVFPProperties().getInj();
                const int table_id = well.wellEcl().injectionControls(
                    well_model_.summaryState()).vfp_table_number;
                if (!inj_vfp.hasTable(table_id)) {
                    return;
                }
                const auto& thp_axis = inj_vfp.getTable(table_id).getTHPAxis();
                limit = std::clamp(limit, static_cast<Scalar>(thp_axis.front()),
                                   static_cast<Scalar>(thp_axis.back()));
            }
            well.setDynamicThpLimit(limit);
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
