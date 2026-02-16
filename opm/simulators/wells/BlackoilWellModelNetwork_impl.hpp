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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellGroupControls.hpp>

#include <fmt/format.h>

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
        well_model_.updateWellControlsAndNetwork(true, dt, deferred_logger);
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
    const auto& network = well_model_.schedule()[episodeIdx].network();
    if (!well_model_.wellsActive() && !network.active()) {
        return {false, 0.0};
    }

    const int iterationIdx = well_model_.simulator().model().newtonMethod().numIterations();
    const auto& comm = well_model_.simulator().vanguard().grid().comm();

    // network related
    Scalar network_imbalance = 0.0;
    bool more_network_update = false;
    if (this->shouldBalance(episodeIdx, iterationIdx) || mandatory_network_balance) {
        OPM_TIMEBLOCK(BalanceNetwork);
        const double dt = well_model_.simulator().timeStepSize();
        // Calculate common THP for subsea manifold well group (item 3 of NODEPROP set to YES)
        const bool well_group_thp_updated = computeWellGroupThp(dt, deferred_logger);
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
            if (!more_network_sub_update) {
                break;
            }

            for (const auto& well : well_model_) {
                if (well->isInjector() || !well->wellEcl().predictionMode()) {
                     continue;
                }

                const auto it = this->node_pressures_.find(well->wellEcl().groupName());
                if (it != this->node_pressures_.end()) {
                    well->prepareWellBeforeAssembling(well_model_.simulator(),
                                                      dt,
                                                      well_model_.groupStateHelper(),
                                                      well_model_.wellState());
                }
            }
            well_model_.updateAndCommunicateGroupData(episodeIdx,
                                                      iterationIdx,
                                                      well_model_.param().nupcol_group_rate_tolerance_,
                                                      /*update_wellgrouptarget*/ true);
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
    const auto& network = well_model_.schedule()[reportStepIdx].network();
    const auto& balance = well_model_.schedule()[reportStepIdx].network_balance();
    const Scalar thp_tolerance = balance.thp_tolerance();

    if (!network.active()) {
        return false;
    }

    auto& well_state = well_model_.wellState();
    auto& group_state = well_model_.groupState();

    bool well_group_thp_updated = false;
    for (const std::string& nodeName : network.node_names()) {
        const bool has_choke = network.node(nodeName).as_choke();
        if (has_choke) {
            const auto& summary_state = well_model_.simulator().vanguard().summaryState();
            const Group& group = well_model_.schedule().getGroup(nodeName, reportStepIdx);

            //TODO: Auto choke combined with RESV control is not supported
            std::vector<Scalar> resv_coeff(Indices::numPhases, 1.0);

            const auto ctrl = group.productionControls(summary_state);
            auto cmode_tmp = ctrl.cmode;
            Scalar target_tmp{0.0};
            bool fld_none = false;
            if (cmode_tmp == Group::ProductionCMode::FLD || cmode_tmp == Group::ProductionCMode::NONE) {
                fld_none = true;
                // Target is set for an ancestor group. Target for autochoke group to be
                // derived via group guide rates
                const Scalar efficiencyFactor = 1.0;
                const Group& parentGroup = well_model_.schedule().getGroup(group.parent(), reportStepIdx);
                auto target = WellGroupControls<Scalar, IndexTraits>::
                    getAutoChokeGroupProductionTargetRate(group,
                                                          parentGroup,
                                                          well_model_.groupStateHelper(),
                                                          resv_coeff,
                                                          efficiencyFactor);
                target_tmp = target.first;
                cmode_tmp = target.second;
            }
            using TargetCalculatorType =  GroupStateHelpers::TargetCalculator<Scalar, IndexTraits>;
            TargetCalculatorType tcalc{well_model_.groupStateHelper(), resv_coeff, group};
            if (!fld_none)
            {
                // Target is set for the autochoke group itself
                target_tmp = tcalc.groupTarget();
            }

            const Scalar orig_target = target_tmp;

            auto mismatch = [&] (auto group_thp) {
                Scalar group_rate(0.0);
                Scalar rate(0.0);
                for (auto& well : well_model_) {
                    std::string well_name = well->name();
                    auto& ws = well_state.well(well_name);
                    if (group.hasWell(well_name)) {
                        well->setDynamicThpLimit(group_thp);
                        const Well& well_ecl = well_model_.eclWells()[well->indexOfWell()];
                        const auto inj_controls = Well::InjectionControls(0);
                        const auto prod_controls = well_ecl.productionControls(summary_state);
                        well->iterateWellEqWithSwitching(well_model_.simulator(),
                                                         dt,
                                                         inj_controls,
                                                         prod_controls,
                                                         well_model_.groupStateHelper(),
                                                         well_state,
                                                         /*fixed_control=*/false,
                                                         /*fixed_status=*/false,
                                                         /*solving_with_zero_rate=*/false);
                        rate = -tcalc.calcModeRateFromRates(ws.surface_rates);
                        group_rate += rate;
                    }
                }
                return (group_rate - orig_target)/orig_target;
            };

            const auto upbranch = network.uptree_branch(nodeName);
            const auto it = this->node_pressures_.find((*upbranch).uptree_node());
            const Scalar nodal_pressure = it->second;
            Scalar well_group_thp = nodal_pressure;

            std::optional<Scalar> autochoke_thp;
            if (auto iter = this->well_group_thp_calc_.find(nodeName);
                iter != this->well_group_thp_calc_.end())
            {
                autochoke_thp = this->well_group_thp_calc_.at(nodeName);
            }

            using WellBhpThpCalculatorType = WellBhpThpCalculator<Scalar, IndexTraits>;
            //Find an initial bracket
            std::array<Scalar, 2> range_initial;
            if (!autochoke_thp.has_value()){
                Scalar min_thp, max_thp;
                // Retrieve the terminal pressure of the associated root of the manifold group
                std::string node_name =  nodeName;
                while (!network.node(node_name).terminal_pressure().has_value()) {
                    auto branch = network.uptree_branch(node_name).value();
                    node_name = branch.uptree_node();
                }
                min_thp = network.node(node_name).terminal_pressure().value();
                WellBhpThpCalculatorType::bruteForceBracketCommonTHP(mismatch, min_thp, max_thp);
                // Narrow down the bracket
                Scalar low1, high1;
                std::array<Scalar, 2> range = {Scalar{0.9}*min_thp, Scalar{1.1}*max_thp};
                std::optional<Scalar> appr_sol;
                WellBhpThpCalculatorType::bruteForceBracketCommonTHP(mismatch,
                                                                     range,
                                                                     low1,
                                                                     high1,
                                                                     appr_sol,
                                                                     0.0,
                                                                     local_deferredLogger);
                min_thp = low1;
                max_thp = high1;
                range_initial = {min_thp, max_thp};
            }

            if (!autochoke_thp.has_value() || autochoke_thp.value() > nodal_pressure) {
                // The bracket is based on the initial bracket or
                // on a range based on a previous calculated group thp
                std::array<Scalar, 2> range = autochoke_thp.has_value() ?
                    std::array<Scalar, 2>{Scalar{0.9} * autochoke_thp.value(),
                                          Scalar{1.1} * autochoke_thp.value()} : range_initial;
                Scalar low, high;
                std::optional<Scalar> approximate_solution;
                const Scalar tolerance1 = thp_tolerance;
                local_deferredLogger.debug("Using brute force search to bracket the group THP");
                const bool finding_bracket = WellBhpThpCalculatorType::
                    bruteForceBracketCommonTHP(mismatch,
                                               range,
                                               low,
                                               high,
                                               approximate_solution,
                                               tolerance1,
                                               local_deferredLogger);

                if (approximate_solution.has_value()) {
                    autochoke_thp = *approximate_solution;
                    local_deferredLogger.debug("Approximate group THP value found: "  +
                                               std::to_string(autochoke_thp.value()));
                } else if (finding_bracket) {
                    const Scalar tolerance2 = thp_tolerance;
                    const int max_iteration_solve = 100;
                    int iteration = 0;
                    autochoke_thp = RegulaFalsiBisection<ThrowOnError>::
                                     solve(mismatch,
                                           low,
                                           high,
                                           max_iteration_solve,
                                           tolerance2,
                                           iteration);
                    local_deferredLogger.debug(" bracket = [" + std::to_string(low) + ", " +
                                               std::to_string(high) + "], " +
                                               "iteration = " + std::to_string(iteration));
                    local_deferredLogger.debug("Group THP value = " + std::to_string(autochoke_thp.value()));
                } else {
                    autochoke_thp.reset();
                    local_deferredLogger.debug("Group THP solve failed due to bracketing failure");
                }
            }
             if (autochoke_thp.has_value()) {
                well_group_thp_calc_[nodeName] = autochoke_thp.value();
                // Note: The node pressure of the auto-choke node is set
                // to well_group_thp in computeNetworkPressures()
                // and must be larger or equal to the pressure of the uptree node of its branch.
                well_group_thp = std::max(autochoke_thp.value(), nodal_pressure);
            }

            for (auto& well : well_model_) {
                std::string well_name = well->name();

                if (well->isInjector() || !well->wellEcl().predictionMode())
                    continue;

                if (group.hasWell(well_name)) {
                    well->setDynamicThpLimit(well_group_thp);
                }
                const auto& ws = well_model_.wellState().well(well->indexOfWell());
                const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                if (thp_is_limit) {
                    well->prepareWellBeforeAssembling(well_model_.simulator(),
                                                      dt,
                                                      well_model_.groupStateHelper(),
                                                      well_model_.wellState());
                }
            }

            // Use the group THP in computeNetworkPressures().
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
