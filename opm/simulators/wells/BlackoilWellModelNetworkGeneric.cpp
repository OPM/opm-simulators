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
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellGroupHelper.hpp>

#include <cassert>

namespace Opm {

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
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
updateActiveState(const int report_step)
{
    const auto& network = well_model_.schedule()[report_step].network();
    if (!network.active()) {
        this->active_ = false;
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
    this->active_ = well_model_.comm().max(network_active);
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
needPreStepRebalance(const int report_step) const
{
    const auto& network = well_model_.schedule()[report_step].network();
    bool network_rebalance_necessary = false;
    for (const auto& well : well_model_.genericWells()) {
        const bool is_partof_network = network.has_node(well->wellEcl().groupName());
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
shouldBalance(const int reportStepIdx, const int iterationIdx) const
{
    // if network is not active, we do not need to balance the network
    const auto& network = well_model_.schedule()[reportStepIdx].network();
    if (!network.active()) {
        return false;
    }

    const auto& balance = well_model_.schedule()[reportStepIdx].network_balance();
    if (balance.mode() == Network::Balance::CalcMode::TimeStepStart) {
        return iterationIdx == 0;
    } else if (balance.mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = well_model_.schedule()[reportStepIdx].nupcol();
        return iterationIdx < nupcol;
    } else {
        // We do not support any other rebalancing modes,
        // i.e. TimeInterval based rebalancing is not available.
        // This should be warned about elsewhere, so we choose to
        // avoid spamming with a warning here.
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
    // Get the network and return if inactive (no wells in network at this time)
    const auto& network = well_model_.schedule()[reportStepIdx].network();
    if (!network.active()) {
        return 0.0;
    }

    const auto previous_node_pressures = node_pressures_;

    node_pressures_ = well_model_.wgHelper().
                        computeNetworkPressures(network,
                                                *well_model_.getVFPProperties().getProd(),
                                                well_model_.comm());

    // here, the network imbalance is the difference between the previous nodal pressure and the new nodal pressure
    Scalar network_imbalance = 0.;
    if (!this->active()) {
        return network_imbalance;
    }

    if (!previous_node_pressures.empty()) {
        for (const auto& [name, new_pressure]: node_pressures_) {
            if (previous_node_pressures.count(name) <= 0) {
                if (std::abs(new_pressure) > network_imbalance) {
                    network_imbalance = std::abs(new_pressure);
                }
                continue;
            }
            const auto pressure = previous_node_pressures.at(name);
            const Scalar change = (new_pressure - pressure);
            if (std::abs(change) > network_imbalance) {
                network_imbalance = std::abs(change);
            }
            // We dampen the nodal pressure change during one iteration since our nodal pressure calculation
            // is somewhat explicit. There is a relative dampening factor applied to the update value, and also
            // the maximum update is limited (to 5 bar by default, can be changed with --network-max-pressure-update-in-bars).
            const Scalar damped_change = std::min(damping_factor * std::abs(change), upper_update_bound);
            const Scalar sign = change > 0 ? 1. : -1.;
            node_pressures_[name] = pressure + sign * damped_change;
        }
    } else {
        for (const auto& [name, pressure]: node_pressures_) {
            if (std::abs(pressure) > network_imbalance) {
                network_imbalance = std::abs(pressure);
            }
        }
    }

    for (auto& well : well_model_.genericWells()) {

        // Producers only, since we so far only support the
        // "extended" network model (properties defined by
        // BRANPROP and NODEPROP) which only applies to producers.
        if (well->isProducer() && well->wellEcl().predictionMode()) {
            const auto it = node_pressures_.find(well->wellEcl().groupName());
            if (it != node_pressures_.end()) {
                // The well belongs to a group with has a network pressure constraint,
                // set the dynamic THP constraint of the well accordingly.
                const Scalar new_limit = it->second;
                well->setDynamicThpLimit(new_limit);
                SingleWellState<Scalar, IndexTraits>& ws = well_model_.wellState()[well->indexOfWell()];
                const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                // TODO: not sure why the thp is NOT updated properly elsewhere
                if (thp_is_limit) {
                    ws.thp = well->getTHPConstraint(well_model_.summaryState());
                }
            }
        }
    }
    return network_imbalance;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
assignNodeValues(std::map<std::string, data::NodeData>& nodevalues,
                 const int reportStepIdx) const
{
    nodevalues.clear();
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

    const auto& network = well_model_.schedule()[reportStepIdx].network();
    if (!network.active()) {
        return;
    }

    auto converged_pressures =
        well_model_.wgHelper().computeNetworkPressures(network,
                                                       *well_model_.getVFPProperties().getProd(),
                                                       well_model_.comm());
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
    const auto& network = well_model_.schedule()[report_step].network();
    if (network.active() && !node_pressures_.empty()) {
        for (auto& well : well_model_.genericWells()) {
            // Producers only, since we so far only support the
            // "extended" network model (properties defined by
            // BRANPROP and NODEPROP) which only applies to producers.
            if (well->isProducer()) {
                const auto it = this->node_pressures_.find(well->wellEcl().groupName());
                if (it != this->node_pressures_.end()) {
                    // The well belongs to a group which has a network nodal pressure,
                    // set the dynamic THP constraint based on the network nodal pressure
                    const Scalar nodal_pressure = it->second;
                    well->setDynamicThpLimit(nodal_pressure);
                }
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelNetworkGeneric<Scalar, IndexTraits>::
operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const
{
    return
           this->node_pressures_ == rhs.node_pressures_
        && this->last_valid_node_pressures_ == rhs.last_valid_node_pressures_;
}

template class BlackoilWellModelNetworkGeneric<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelNetworkGeneric<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
