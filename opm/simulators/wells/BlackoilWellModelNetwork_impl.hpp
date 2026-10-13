/*
  Copyright 2016 - 2019 SINTEF Digital, Mathematics & Cybernetics.
  Copyright 2016 - 2018 Equinor ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 Norce AS

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

#ifndef OPM_BLACKOILWELLMODEL_NETWORK_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NETWORK_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_BLACKOILWELLMODEL_NETWORK_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelNetwork.hpp>
#endif

#include <opm/common/TimingMacros.hpp>
#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <optional>
#include <tuple>

namespace Opm {

template<typename TypeTag>
BlackoilWellModelNetwork<TypeTag>::
BlackoilWellModelNetwork(BlackoilWellModel<TypeTag>& well_model)
    : BaseType(well_model)
    , well_model_(well_model)
{}

template<typename TypeTag>
void
BlackoilWellModelNetwork<TypeTag>::
doPreStepRebalance(DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const double dt = well_model_.simulator().timeStepSize();
    // TODO: should we also have the group and network backed-up
    //       here in case the solution did not get converged?
    auto& well_state = well_model_.wellState();

    const bool changed_well_group =
        well_model_.updateWellControlsAndNetwork(/*mandatory_network_balance=*/true,
                                                 dt,
                                                 deferred_logger);
    well_model_.assembleWellEqWithoutIteration(dt);
    const bool converged =
        well_model_.getWellConvergence(well_model_.B_avg(), true).converged() &&
        !changed_well_group;

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    for (auto& well : this->well_model_) {
        well->solveEqAndUpdateWellState(well_model_.simulator(),
                                        well_model_.groupStateHelper(),
                                        well_state);
    }
    OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModelNetwork::doPreStepRebalance() failed: ",
                                well_model_.simulator().vanguard().grid().comm());

    if (!converged) {
        deferred_logger.warning("Initial (pre-step) network balance did not converge.");
    }
}

template<typename TypeTag>
std::tuple<bool, typename BlackoilWellModelNetwork<TypeTag>::Scalar>
BlackoilWellModelNetwork<TypeTag>::
update(const bool mandatory_network_balance,
       DeferredLogger& deferred_logger,
       const bool relax_network_tolerance)
{
    OPM_TIMEFUNCTION();
    const int episodeIdx = well_model_.simulator().episodeIndex();
    if (!well_model_.wellsActive() && !details::anyNetworkActive(well_model_.schedule(), episodeIdx)) {
        return {/*more_network_update=*/false, /*network_imbalance=*/0.0};
    }

    const auto& comm = well_model_.simulator().vanguard().grid().comm();

    // network related
    Scalar network_imbalance = 0.0;
    bool more_network_update = false;
    if (this->shouldBalance(episodeIdx) || mandatory_network_balance) {
        OPM_TIMEBLOCK(BalanceNetwork);
        const double dt = well_model_.simulator().timeStepSize();
        // Calculate common THP for subsea manifold well group (item 3 of NODEPROP set to YES)
        const bool well_group_thp_updated = computeWellGroupThp(dt, deferred_logger);
        if (well_group_thp_updated) {
            // The network pressures must follow from the rates of the wells
            // just solved against the new common THP.
            well_model_.updateAndCommunicateGroupData(episodeIdx, /*update_wellgrouptarget*/ true);
        }
        const int max_number_of_sub_iterations =
            well_model_.param().network_max_sub_iterations_;
        const Scalar network_pressure_update_damping_factor =
            well_model_.param().network_pressure_update_damping_factor_;
        const Scalar network_max_pressure_update =
            well_model_.param().network_max_pressure_update_in_bars_ * unit::barsa;
        bool more_network_sub_update = false;
        for (int i = 0; i < max_number_of_sub_iterations; i++) {
            const auto local_network_imbalance =
                this->updatePressures(episodeIdx,
                                      network_pressure_update_damping_factor,
                                      network_max_pressure_update);
            network_imbalance = comm.max(local_network_imbalance);
            const auto& balance = well_model_.schedule()[episodeIdx].network_balance();
            constexpr Scalar relaxation_factor = 10.0;
            const Scalar tolerance =
                relax_network_tolerance ? relaxation_factor * balance.pressure_tolerance()
                                        : balance.pressure_tolerance();
            more_network_sub_update = this->active() && network_imbalance > tolerance;
#ifdef RESERVOIR_COUPLING_ENABLED
            if (well_model_.isReservoirCouplingMaster()) {
                // Tight reservoir-coupling: exchange node pressures and slave rates with the
                // slaves per inner sub-iteration. Must run before the break below so the converged
                // sub-iteration still feeds back the slaves' rates.
                well_model_.rescoupHelper().maybeExchangeNetworkSubIterationWithSlaves();
            }
#endif
            if (!more_network_sub_update) {
                break;
            }

            for (const auto& well : well_model_) {
                if (!well->wellEcl().predictionMode()) {
                     continue;
                }

                const auto domain = details::domainForWell(*well);

                if (!domain.has_value()) {
                    continue;
                }

                const auto it = this->nodePressures(*domain).find(well->wellEcl().groupName());
                if (it != this->nodePressures(*domain).end()) {
                    well->prepareWellBeforeAssembling(well_model_.simulator(),
                                                      dt,
                                                      well_model_.groupStateHelper(),
                                                      well_model_.wellState());
                    // Option B: after re-solving at the current network THP, update
                    // ws.well_potentials for injection wells.  The rate_less_than_potential
                    // check in WellConstraints::activeInjectionConstraint compares current
                    // injection rates against ws.well_potentials to decide whether switching
                    // to THP mode would increase or decrease injection.  The potentials are
                    // normally computed once per timestep at the static WCONINJE THP and are
                    // stale during network iterations.  Refreshing them here (from the rate
                    // the well just solved to under the current network THP) makes the check
                    // accurate for subsequent outer iterations.
                    if (well->isInjector()) {
                        auto& ws = well_model_.wellState().well(well->indexOfWell());
                        if (ws.injection_cmode == Well::InjectorCMode::THP) {
                            const int np = well_model_.numPhases();
                            for (int p = 0; p < np; ++p) {
                                ws.well_potentials[p] =
                                    std::max(Scalar{0.0}, ws.surface_rates[p]);
                            }
                        }
                    }
                }
            }
            well_model_.updateAndCommunicateGroupData(episodeIdx, /*update_wellgrouptarget*/ true);
        }
        more_network_update = more_network_sub_update || well_group_thp_updated;
    }
    return { more_network_update, network_imbalance };
}

template <typename TypeTag>
bool
BlackoilWellModelNetwork<TypeTag>::
computeWellGroupThp(const double dt, DeferredLogger& local_deferredLogger)
{
    OPM_TIMEFUNCTION();
    const int reportStepIdx = well_model_.simulator().episodeIndex();
    // This function is only relevant for auto-choke groups, and
    // therefore as of now only relevant for the production network.
    // \TODO: If we later also want to support auto-choke groups in the
    // injection network, we should change this function also.
    const auto& network = well_model_.schedule()[reportStepIdx].network();
    const auto& balance = well_model_.schedule()[reportStepIdx].network_balance();
    const Scalar thp_tolerance = balance.thp_tolerance();

    if (!network.active()) {
        return false;
    }

    const auto& simulator = well_model_.simulator();
    const auto& summary_state = simulator.vanguard().summaryState();
    const auto& comm = simulator.vanguard().grid().comm();
    const auto& group_state_helper = well_model_.groupStateHelper();
    auto& well_state = well_model_.wellState();
    auto& group_state = well_model_.groupState();

    bool well_group_thp_updated = false;
    for (const std::string& nodeName : network.node_names()) {
        const bool has_choke = network.node(nodeName).as_choke();
        if (has_choke) {
            const Group& group = well_model_.schedule().getGroup(nodeName, reportStepIdx);

            //TODO: Auto choke combined with RESV control is not supported
            std::vector<Scalar> resv_coeff(Indices::numPhases, 1.0);

            // The group's wells are on THP control, so its current control mode
            // may have been reset to NONE; the mode set in GCONPROD then applies.
            auto cmode = group_state.production_control(nodeName);
            if (cmode == Group::ProductionCMode::FLD || cmode == Group::ProductionCMode::NONE) {
                cmode = group.productionControls(summary_state).cmode;
            }
            Scalar target{0.0};
            if (cmode == Group::ProductionCMode::FLD || cmode == Group::ProductionCMode::NONE) {
                // Target is set for an ancestor group. Target for autochoke group to be
                // derived via group guide rates
                const Scalar efficiencyFactor = 1.0;
                const Group& parentGroup = well_model_.schedule().getGroup(group.parent(), reportStepIdx);
                std::tie(target, cmode) = group_state_helper.
                    getAutoChokeGroupProductionTargetRate(group,
                                                          parentGroup,
                                                          resv_coeff,
                                                          efficiencyFactor);
            }
            else {
                // Target is set for the autochoke group itself
                target = group_state_helper.getProductionGroupTargetForMode(group, cmode);
            }
            constexpr auto rate_modes = std::array{Group::ProductionCMode::ORAT,
                                                   Group::ProductionCMode::WRAT,
                                                   Group::ProductionCMode::GRAT,
                                                   Group::ProductionCMode::LRAT,
                                                   Group::ProductionCMode::RESV};
            const bool has_target = target > 0.0
                && std::ranges::find(rate_modes, cmode) != rate_modes.end();

            // The wells cannot produce against less than the node's own network
            // pressure, which does not include the manifold's choke.
            Scalar nodal_pressure{0.0};
            if (const auto it = this->node_pressures_.find(nodeName);
                it != this->node_pressures_.end())
            {
                nodal_pressure = it->second;
            }
            else {
                // No network solution yet; start from the fixed pressure node above.
                std::string node_name = nodeName;
                while (!network.node(node_name).terminal_pressure().has_value()) {
                    node_name = network.uptree_branch(node_name).value().uptree_node();
                }
                nodal_pressure = network.node(node_name).terminal_pressure().value();
            }

            using TargetCalculatorType = GroupStateHelpers::TargetCalculator<Scalar, IndexTraits>;
            const TargetCalculatorType tcalc{group_state_helper, resv_coeff, cmode};
            // Group rate with the producers at a common THP, from their inflow at
            // the current reservoir state. The well state is left untouched.
            const auto group_rate = [&](const Scalar thp) {
                Scalar rate{0.0};
                for (const auto& well : well_model_) {
                    if (!well->isProducer() || !group.hasWell(well->name())) {
                        continue;
                    }
                    const auto& ws = well_state.well(well->indexOfWell());
                    if (ws.status != Well::Status::OPEN) {
                        continue;
                    }
                    // A well without a VFP table cannot be put under THP
                    // control, so the common THP is not imposed on it.
                    const auto rates = well->wellEcl().vfp_table_number() > 0
                                       && well->wellEcl().predictionMode()
                        ? well->computeWellRatesAtThpLimitProd(simulator, group_state_helper,
                                                               well_state, thp)
                        : ws.surface_rates;
                    if (well_state.wellIsOwned(well->indexOfWell(), well->name())) {
                        rate -= well->wellEcl().getEfficiencyFactor() * ws.efficiency_scaling_factor
                            * tcalc.calcModeRateFromRates(rates);
                    }
                }
                return comm.sum(rate);
            };

            std::optional<Scalar> manifold_thp;
            if (has_target) {
                const auto mismatch = [&](const Scalar thp) {
                    return (group_rate(thp) - target) / target;
                };
                // The group rate falls as the common THP rises. With the choke
                // open the wells produce against the nodal pressure.
                Scalar low = nodal_pressure;
                if (mismatch(low) <= thp_tolerance) {
                    manifold_thp = low;
                }
                else {
                    // Bracket the THP meeting the target, from the last solution if any.
                    // Steps doubling from 1 bar reach a thousand bar within ten attempts.
                    constexpr Scalar initial_step = 1.0 * unit::barsa;
                    constexpr int max_bracket_attempts = 10;
                    const auto cached = this->well_group_thp_calc_.find(nodeName);
                    Scalar high = (cached != this->well_group_thp_calc_.end() && cached->second > low)
                        ? cached->second
                        : low + initial_step;
                    Scalar step = std::max(high - low, initial_step);
                    Scalar mismatch_high = mismatch(high);
                    for (int attempt = 0; mismatch_high > 0.0 && attempt < max_bracket_attempts; ++attempt) {
                        low = high;
                        step *= 2;
                        high = low + step;
                        mismatch_high = mismatch(high);
                    }
                    if (mismatch_high > 0.0) {
                        local_deferredLogger.debug(fmt::format("Sub-sea manifold {}: no common THP "
                                                               "brings the group down to its target",
                                                               nodeName));
                    }
                    else {
                        int iterations = 0;
                        manifold_thp = RegulaFalsiBisection<ContinueOnError>::
                            solve(mismatch, low, high, balance.thp_max_iter(),
                                  thp_tolerance, iterations);
                    }
                }
            }

            Scalar well_group_thp = nodal_pressure;
            if (manifold_thp.has_value()) {
                this->well_group_thp_calc_[nodeName] = *manifold_thp;
                well_group_thp = std::max(*manifold_thp, nodal_pressure);
            }
            local_deferredLogger.debug(fmt::format("Sub-sea manifold {}: common THP {:.3f} bar, "
                                                   "nodal pressure {:.3f} bar, {} target {:.3f} per day",
                                                   nodeName,
                                                   well_group_thp / unit::barsa,
                                                   nodal_pressure / unit::barsa,
                                                   Group::ProductionCMode2String(cmode),
                                                   target * unit::day));

            for (auto& well : well_model_) {
                if (!well->isProducer() || !well->wellEcl().predictionMode()
                    || !group.hasWell(well->name())) {
                    continue;
                }
                // As above: a well without a VFP table cannot be put under THP
                // control, so the auto-choke group THP is not imposed on it.
                if (well->wellEcl().vfp_table_number() > 0) {
                    well->setDynamicThpLimit(well_group_thp);
                }
                // Wells on rate control may have to switch to the new THP limit
                // before the network sees their rates.
                well->prepareWellBeforeAssembling(simulator, dt, group_state_helper, well_state);
            }

            // The wells of the group produce against this THP, see producerThpLimit().
            const auto& current_well_group_thp = group_state.is_autochoke_group(nodeName)
                ? group_state.well_group_thp(nodeName)
                : 1e30;
            if (std::abs(current_well_group_thp - well_group_thp) > balance.pressure_tolerance()) {
                well_group_thp_updated = true;
                group_state.update_well_group_thp(nodeName, well_group_thp);
            }
        }
    }
    return well_group_thp_updated;
}

} // namespace Opm

#endif // OPM_BLACKOILWELLMODEL_NETWORK_IMPL_HEADER_INCLUDED
