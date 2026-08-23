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

#ifndef OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>

#include <opm/output/data/Groups.hpp>

#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/wells/NetworkAndersonAcceleration.hpp>
#include <opm/simulators/wells/NetworkNodePressureUpdater.hpp>
#include <opm/simulators/wells/NetworkSystem.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>

namespace Opm {
    class Schedule;
    class UnitSystem;
    template<class Scalar, class IndexTraits> class BlackoilWellModelGeneric;
    template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
    template<typename Scalar> class VFPInjProperties;
    template<typename Scalar> class VFPProdProperties;
}

namespace Opm {

namespace details {

    enum class NetworkDomain : std::size_t {
        Production = 0,
        InjectionGas,
        InjectionWater,
        Count
    };

    constexpr std::size_t domainIndex(const NetworkDomain domain)
    {
        return static_cast<std::size_t>(domain);
    }

    struct ActiveNetworkDescriptor {
        NetworkDomain domain;
        std::reference_wrapper<const Network::ExtNetwork> network;
    };

    /// The network a well belongs to, or nullopt for a well that is in none of them.
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

    /// The injected phase of an injection network domain, nullopt for the production one.
    std::optional<Phase> injectionPhaseForDomain(const NetworkDomain domain);

    /// Helper to check if any network (production, gas injection, water injection) is active at a given time step.
    bool anyNetworkActive(const Schedule& schedule, const int timeStepIdx);

    /// Helper to get all active networks (production, gas injection, water injection) at a given time step.
    std::vector<ActiveNetworkDescriptor>
    activeNetworks(const Schedule& schedule, const int timeStepIdx);

} // namespace details


/// Class for handling the blackoil well network model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelNetworkGeneric
{
public:
    BlackoilWellModelNetworkGeneric(BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model);

    virtual ~BlackoilWellModelNetworkGeneric() = default;

    /// return true if network is active (at least one network well in prediction mode)
    bool active() const
    { return active_; }

    const std::map<std::string, Scalar>&
    nodePressures() const { return node_pressures_; }

    // do not use, only needed for serialization testing
    void setNodePressures(const std::map<std::string, Scalar>& values)
    { node_pressures_ = values; }

    void setFromRestart(const std::optional<std::map<std::string, double>>& restart_pressures);

    //! \brief Initialize wells according to network configuration.
    void initialize(const int report_step);

    //! \brief Initialize a single well according to network configuration.
    void initializeWell(WellInterfaceGeneric<Scalar,IndexTraits>& well);

    /// Checks if network is active (at least one network well on prediction).
    void updateActiveState(const int report_step);

    /// Checks if there are reasons to perform a pre-step network re-balance.
    /// (Currently, the only reasons are network well status changes.)
    /// (TODO: Consider if adding network change events would be helpful.)
    bool needPreStepRebalance(const int report_step) const;

    /// Checks if we shall perform a network re-balance.
    /// This is typically controlled by the NETBALAN keyword.
    bool shouldBalance(const int reportStepIndex) const;
    /// Checks if we will perform a network re-balance on the next Newton iteration.
    bool willBalanceOnNextIteration(const int reportStepIndex) const;

    /// Recompute the node pressures from the current leaf rates and move the applied
    /// pressures towards them. With use_secant the per-node update is a secant step on
    /// r(P) = P_computed(P) - P (clipped to the bracket [P, P_computed], which contains
    /// the fixed point when the wells respond monotonically); otherwise, and as fallback,
    /// the change is damped by damping_factor and capped at update_upper_bound.
    /// Returns the largest |r| over the nodes.
    Scalar updatePressures(const int reportStepIdx,
                           const Scalar damping_factor,
                           const Scalar update_upper_bound,
                           const bool use_secant = false,
                           const bool secant_for_production = false,
                           const int anderson_depth = 0);

    /// Forget the secant history; call at the start of every time step.
    void beginTimeStep()
    {
        for (auto& u : pressure_updaters_) {
            u.clear();
        }
        for (auto& a : pressure_accelerators_) {
            a.clear();
        }
    }

    /// Fill the production node/branch values (GPR, GPRB, ..., with the converged
    /// pressures recomputed from the current rates) and the gas/water injection
    /// network node and branch values (GPRG, GPRW) of `values`.
    void assignNodeAndBranchValues(data::GroupAndNetworkValues& values,
                                   const int reportStepIdx) const;

    void commitState()
    {
        this->last_valid_node_pressures_ = this->node_pressures_;
        this->last_valid_branch_data_ = this->branch_data_;
        this->last_valid_domain_node_pressures_ = this->domain_node_pressures_;
        this->last_valid_domain_branch_data_ = this->domain_branch_data_;
    }

    void resetState()
    {
        this->node_pressures_ = this->last_valid_node_pressures_;
        this->branch_data_ = this->last_valid_branch_data_;
        this->domain_node_pressures_ = this->last_valid_domain_node_pressures_;
        this->domain_branch_data_ = this->last_valid_domain_branch_data_;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(node_pressures_);
        serializer(last_valid_node_pressures_);
        serializer(branch_data_);
        serializer(last_valid_branch_data_);
        serializer(domain_node_pressures_);
        serializer(last_valid_domain_node_pressures_);
        serializer(domain_branch_data_);
        serializer(last_valid_domain_branch_data_);
    }

    bool operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const;


    /// Solve an injection network simultaneously in its pressures and rates
    /// instead of relaxing the node pressures against the wells. Off by default;
    /// --network-solver=newton turns it on. Requires the injectors' implicit IPR
    /// to have been refreshed, which the templated caller does.
    void useNewtonSolver(const bool on) { newton_solver_ = on; }
    bool usesNewtonSolver() const { return newton_solver_; }

    /// Assemble the network Jacobian from the VFP table derivatives instead of
    /// differencing the residual.
    void useAnalyticJacobian(const bool on) { analytic_jacobian_ = on; }

    /// Let the network hold the group's total and place the split itself, rather
    /// than taking each group-controlled well's rate as fixed.
    void useNetworkGroupControl(const bool on) { network_group_control_ = on; }
    void useNetworkAutochoke(const bool on) { network_autochoke_ = on; }
    /// Per local well, the hydrostatic correction its tubing table needs;
    /// computed on the typed side, where the well's density lives.
    void setWellVfpDp(const std::string& well, const Scalar dp) { well_vfp_dp_[well] = dp; }
    bool networkAutochoke() const { return network_autochoke_; }
    /// Answer the gas lift optimiser's trials from the network instead of
    /// from a well solve at a fixed thp.
    void useGasLiftNetworkResponse(const bool on) { gaslift_network_response_ = on; }
    bool gasLiftNetworkResponse() const { return gaslift_network_response_; }

    /// What a well would produce, and at what bhp, with its lift gas set to
    /// alq -- with every node pressure responding. Water, oil, gas, bhp, all
    /// SI and production positive; nullopt if the well is in no solved tree
    /// or the trial does not converge.
    std::optional<std::array<Scalar, 4>>
    gasLiftTrial(const std::string& well, const Scalar alq) const;

    /// Write each network system that fails to converge, for replay in
    /// tests/test_networksolve.cpp. Empty disables it.
    void dumpNetworkFailuresTo(const std::string& prefix) { network_dump_prefix_ = prefix; }

protected:
    /// Result of one network pressure evaluation for one network (domain).
    struct NetworkPressures
    {
        std::map<std::string, Scalar> node_pressures;
        std::map<std::string, data::BranchData> branch_data;
        // Nodes (and their descendants) whose VFP lookup has no solution; their
        // node_pressures entries are placeholders and must not be used.
        std::set<std::string> invalid_nodes;
    };

    NetworkPressures
    computePressures(const Network::ExtNetwork& network,
                     const VFPProdProperties<Scalar>& vfp_prod_props,
                     const UnitSystem& unit_system,
                     const int reportStepIdx,
                     const Parallel::Communication& comm) const;

    NetworkPressures
    computePressures(const Network::ExtNetwork& network,
                     const VFPInjProperties<Scalar>& vfp_inj_props,
                     const UnitSystem& unit_system,
                     const int reportStepIdx,
                     const Parallel::Communication& comm,
                     const Phase injectionPhase) const;

    void updateActiveStateImpl(const Network::ExtNetwork& network);

    static constexpr details::NetworkDomain productionNetworkDomain()
    {
        return details::NetworkDomain::Production;
    }

    const std::map<std::string, Scalar>& nodePressures(const details::NetworkDomain domain) const
    {
        return domain_node_pressures_[details::domainIndex(domain)];
    }

    std::map<std::string, Scalar>& nodePressures(const details::NetworkDomain domain)
    {
        return domain_node_pressures_[details::domainIndex(domain)];
    }

    const std::map<std::string, data::BranchData>& branchData(const details::NetworkDomain domain) const
    {
        return domain_branch_data_[details::domainIndex(domain)];
    }

    std::map<std::string, data::BranchData>& branchData(const details::NetworkDomain domain)
    {
        return domain_branch_data_[details::domainIndex(domain)];
    }

    const std::set<std::string>& invalidNodes(const details::NetworkDomain domain) const
    {
        return domain_invalid_nodes_[details::domainIndex(domain)];
    }

    std::set<std::string>& invalidNodes(const details::NetworkDomain domain)
    {
        return domain_invalid_nodes_[details::domainIndex(domain)];
    }

    void syncLegacyProductionState_()
    {
        this->node_pressures_ = this->nodePressures(productionNetworkDomain());
        this->branch_data_ = this->branchData(productionNetworkDomain());
    }

    void syncProductionDomainState_()
    {
        this->nodePressures(productionNetworkDomain()) = this->node_pressures_;
        this->branchData(productionNetworkDomain()) = this->branch_data_;
    }

    bool active_{false};
    BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model_;

    // Network pressures for output and initialization
    std::map<std::string, Scalar> node_pressures_;
    // Network branch pressure drops and flow rates for output (outlet branch for production network, inlet branch for injection network)
    std::map<std::string, data::BranchData> branch_data_;
    // Domain-scoped pressure state to avoid collisions between production and injection networks.
    std::array<std::map<std::string, Scalar>, details::domainIndex(details::NetworkDomain::Count)> domain_node_pressures_;
    std::array<std::map<std::string, data::BranchData>, details::domainIndex(details::NetworkDomain::Count)> domain_branch_data_;
    // Nodes without a valid VFP solution in the last evaluation (per domain); not serialized,
    /// Node pressures from the simultaneous solve, or nullopt if it did not
    /// converge -- in which case the caller keeps the fixed-point result.
    std::optional<std::map<std::string, Scalar>>
    newtonNodePressures(const Network::ExtNetwork& network,
                        const Phase injection_phase,
                        const int reportStepIdx,
                        const Network::Node& root) const;

    /// The same for a production network. A rate is three numbers instead of
    /// one and the wells are producers, which is the whole difference.
    std::optional<std::map<std::string, Scalar>>
    newtonProductionNodePressures(const Network::ExtNetwork& network,
                                  const int reportStepIdx,
                                  const Network::Node& root) const;

    bool newton_solver_ = false;
    bool analytic_jacobian_ = false;
    bool network_group_control_ = false;
    bool network_autochoke_ = false;
    std::map<std::string, Scalar> well_vfp_dp_;
    /// Last production solve per tree root: the inputs it was built from and
    /// what it gave. Inside a network sub-loop the wells are frozen, so the
    /// same inputs come back sub-iteration after sub-iteration.
    struct SolvedTree {
        std::vector<Scalar> inputs;
        std::map<std::string, Scalar> pressures;
        std::vector<std::string> order;                     // node names by index
        std::shared_ptr<const NetworkSolve::ProductionSystem<Scalar>> system;
    };
    mutable std::map<std::string, SolvedTree> last_production_solve_;
    bool gaslift_network_response_ = false;
    std::string network_dump_prefix_;
    mutable int network_dumps_written_ = 0;

    // recomputed on every updatePressures().
    std::array<std::set<std::string>, details::domainIndex(details::NetworkDomain::Count)> domain_invalid_nodes_;
    int invalid_nodes_report_step_{-1};
    // Per node: state of the bracketing/secant pressure update. Not serialized.
    std::array<std::map<std::string, NodePressureUpdater<Scalar>>,
               details::domainIndex(details::NetworkDomain::Count)> pressure_updaters_;
    // Optional whole-vector acceleration, one per domain (off by default).
    std::array<NetworkAndersonAccelerator<Scalar>,
               details::domainIndex(details::NetworkDomain::Count)> pressure_accelerators_;
    // Valid network pressures for output and initialization for safe restart after failed iterations
    std::map<std::string, Scalar> last_valid_node_pressures_;
    // Valid network branch pressure drops and flow rates for output (outlet branch for production network, inlet branch for injection network) for safe restart after failed iterations
    std::map<std::string, data::BranchData> last_valid_branch_data_;
    std::array<std::map<std::string, Scalar>, details::domainIndex(details::NetworkDomain::Count)> last_valid_domain_node_pressures_;
    std::array<std::map<std::string, data::BranchData>, details::domainIndex(details::NetworkDomain::Count)> last_valid_domain_branch_data_;
};

} // namespace Opm

#endif
