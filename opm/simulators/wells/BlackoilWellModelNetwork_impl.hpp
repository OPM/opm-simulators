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
#include <opm/simulators/wells/WellHelpers.hpp>
#endif

#include <opm/common/TimingMacros.hpp>
#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <optional>

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
        this->noteNetworkTimeStep(well_model_.simulator().time(), dt);
        // Calculate common THP for subsea manifold well group (item 3 of NODEPROP set to YES)
        const bool well_group_thp_updated = computeWellGroupThp(dt, deferred_logger);
        const int max_number_of_sub_iterations =
            well_model_.param().network_max_sub_iterations_;
        const Scalar network_pressure_update_damping_factor =
            well_model_.param().network_pressure_update_damping_factor_;
        const Scalar network_max_pressure_update =
            well_model_.param().network_max_pressure_update_in_bars_ * unit::barsa;
        const auto& secant_mode = well_model_.param().network_pressure_update_secant_;
        if (secant_mode != "injection" && secant_mode != "all" && secant_mode != "none") {
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Invalid value '" + secant_mode + "' for --network-pressure-update-secant; "
                             "expected injection, all or none", deferred_logger);
        }
        const bool use_secant = secant_mode != "none";
        const bool secant_production = secant_mode == "all";
        const auto& accel_mode = well_model_.param().network_pressure_update_acceleration_;
        if (accel_mode != "none" && accel_mode != "anderson") {
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Invalid value '" + accel_mode + "' for --network-pressure-update-acceleration; "
                             "expected none or anderson", deferred_logger);
        }
        const int anderson_depth = (accel_mode == "anderson")
            ? well_model_.param().network_anderson_depth_ : 0;
        // Only a deck with both a production and an injection network can have the
        // producers re-solved here feed an injection group target in the same
        // sub-iteration; nothing else pays for the extra group update.
        const auto active_networks = details::activeNetworks(well_model_.schedule(), episodeIdx);
        const auto has_domain = [&active_networks](const bool production)
        {
            return std::any_of(active_networks.begin(), active_networks.end(),
                               [production](const auto& n)
                               { return (n.domain == details::NetworkDomain::Production) == production; });
        };
        const bool refresh_group_data_between =
            has_domain(/*production=*/true) && has_domain(/*production=*/false);
        const auto& proxy_mode = well_model_.param().network_well_proxy_;
        if (proxy_mode != "none" && proxy_mode != "ipr") {
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Invalid value '" + proxy_mode + "' for --network-well-proxy; "
                             "expected none or ipr", deferred_logger);
        }
        if (proxy_mode == "ipr") {
            // Get the network close to balance against the frozen well linearisation first;
            // the loop below then does the real well solves from a much better starting point.
            this->proxyBalance(episodeIdx, dt,
                               well_model_.param().network_well_proxy_max_iterations_,
                               network_pressure_update_damping_factor,
                               network_max_pressure_update,
                               use_secant, secant_production, deferred_logger);
        }
        const auto& solver_mode = well_model_.param().network_solver_;
        if (solver_mode != "fixedpoint" && solver_mode != "newton") {
            OPM_DEFLOG_THROW(std::runtime_error,
                             "Invalid value '" + solver_mode + "' for --network-solver; "
                             "expected fixedpoint or newton", deferred_logger);
        }
        this->useNewtonSolver(solver_mode == "newton");
        this->useAnalyticJacobian(well_model_.param().network_analytic_jacobian_);
        this->useNetworkGroupControl(well_model_.param().network_group_control_);
        this->useNetworkAutochoke(well_model_.param().network_autochoke_);
        this->useNetworkComplementarity(well_model_.param().network_complementarity_);
        this->useGasLiftNetworkResponse(well_model_.param().gaslift_network_response_);
        this->dumpNetworkFailuresTo(well_model_.param().network_dump_failures_);
        if (solver_mode == "newton") {
            // The simultaneous solve needs every well's rate response to its own
            // bhp. That is the implicit IPR, which the well solve maintains only
            // where its own control logic happens to need it -- never for
            // injectors, and for producers only on some paths. Refresh it here
            // for all of them, or a network solve arrives with a well it cannot
            // linearise and hands the whole network back to the relaxed update.
            for (const auto& well : well_model_) {
                if (well->wellEcl().predictionMode()) {
                    well->updateIPRImplicit(well_model_.simulator(),
                                            well_model_.groupStateHelper(),
                                            well_model_.wellState());
                    // The tubing table's datum is not the well's reference
                    // depth; the well's thp evaluation corrects for it and
                    // the network system has to apply the same.
                    const int table = well->wellEcl().vfp_table_number();
                    Scalar dp = 0.0;
                    if (table > 0) {
                        const auto& vfp = well_model_.getVFPProperties();
                        const Scalar datum = well->isInjector()
                            ? vfp.getInj()->getTable(table).getDatumDepth()
                            : vfp.getProd()->getTable(table).getDatumDepth();
                        // wellhelpers::computeHydrostaticCorrection, inline.
                        dp = well->refDensity() * well->gravity() * (datum - well->refDepth());
                    }
                    this->setWellVfpDp(well->name(), dp);
                }
            }
        }

        bool more_network_sub_update = false;
        for (int i = 0; i < max_number_of_sub_iterations; i++) {
            const auto local_network_imbalance =
                this->updatePressures(episodeIdx,
                                      network_pressure_update_damping_factor,
                                      network_max_pressure_update,
                                      use_secant,
                                      secant_production,
                                      anderson_depth);
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

            // Re-solve the producers before the injectors: an injection group target
            // (GCONINJE VREP/REIN) is a function of the produced voidage, so the
            // injectors must see this sub-iteration's production, not the previous
            // one's. The intermediate group update is skipped when no injection
            // network is active, which keeps production-only runs unchanged.
            const auto resolve = [&](const bool injectors)
            {
                for (const auto& well : well_model_) {
                    if (well->isInjector() != injectors || !well->wellEcl().predictionMode()) {
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
                    }
                }
            };
            resolve(/*injectors=*/false);
            if (refresh_group_data_between) {
                well_model_.updateAndCommunicateGroupData(episodeIdx, /*update_wellgrouptarget*/ true);
            }
            resolve(/*injectors=*/true);
            well_model_.updateAndCommunicateGroupData(episodeIdx, /*update_wellgrouptarget*/ true);
        }
        more_network_update = more_network_sub_update || well_group_thp_updated;
    }
    return { more_network_update, network_imbalance };
}

template<typename TypeTag>
std::optional<typename BlackoilWellModelNetwork<TypeTag>::Scalar>
BlackoilWellModelNetwork<TypeTag>::
proxyInjectionRate(WellInterface<TypeTag>& well,
                   const int phase_pos,
                   DeferredLogger& deferred_logger) const
{
    const auto& summary_state = well_model_.simulator().vanguard().summaryState();
    const auto& ws = well_model_.wellState().well(well.indexOfWell());
    const auto& ipr_a = ws.implicit_ipr_a;
    const auto& ipr_b = ws.implicit_ipr_b;
    if (ipr_a.empty() || ipr_b.empty()) {
        return std::nullopt;
    }

    // The well index linearisation of the converged well equation: rates linear in bhp,
    // exact at the bhp the well was last solved at. Both arrays are phase-indexed.
    auto frates = [&ipr_a, &ipr_b](const Scalar bhp)
    {
        std::vector<Scalar> rates(ipr_a.size(), 0.0);
        for (std::size_t p = 0; p < rates.size(); ++p) {
            rates[p] = ipr_b[p] * bhp - ipr_a[p];
        }
        return rates;
    };

    // getTHPConstraint() returns the dynamic limit updatePressures() has just applied,
    // so this is the rate at the current node pressure.
    const auto bhp = WellBhpThpCalculator(well)
        .computeBhpAtThpLimitInj(frates, summary_state, well.refDensity(),
                                 1e-6, 50, /*throwOnError=*/false, deferred_logger);
    if (!bhp.has_value()) {
        return std::nullopt;
    }
    const auto controls = well.wellEcl().injectionControls(summary_state);
    const Scalar rate = frates(std::min(*bhp, static_cast<Scalar>(controls.bhp_limit)))[phase_pos];
    return std::max(rate, Scalar{0});
}

template<typename TypeTag>
typename BlackoilWellModelNetwork<TypeTag>::Scalar
BlackoilWellModelNetwork<TypeTag>::
proxyBalance(const int episodeIdx,
             const double,
             const int max_iterations,
             const Scalar damping_factor,
             const Scalar max_pressure_update,
             const bool use_secant,
             const bool secant_production,
             DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const auto& comm = well_model_.simulator().vanguard().grid().comm();
    const auto& balance = well_model_.schedule()[episodeIdx].network_balance();
    auto& group_state = well_model_.groupStateHelper().groupState();

    // Refresh the well index linearisation once, at the state the wells were last solved
    // in. Only injectors need it: producers keep the rates the well solve gave them.
    for (const auto& well : well_model_) {
        if (well->isInjector() && well->wellEcl().predictionMode()) {
            well->updateIPRImplicit(well_model_.simulator(),
                                    well_model_.groupStateHelper(),
                                    well_model_.wellState());
        }
    }

    Scalar imbalance = 0.0;
    for (int it = 0; it < max_iterations; ++it) {
        imbalance = comm.max(this->updatePressures(episodeIdx, damping_factor,
                                                  max_pressure_update, use_secant,
                                                  secant_production));
        if (!this->active() || imbalance <= balance.pressure_tolerance()) {
            break;
        }
        // Predict the leaf rates at the pressures just applied. Producers and the
        // production network are left alone; only the injection leaves are refreshed.
        for (const auto& network : details::activeNetworks(well_model_.schedule(), episodeIdx)) {
            const auto phase = details::injectionPhaseForDomain(network.domain);
            if (!phase.has_value()) {
                continue;
            }
            const int phase_pos = (*phase == Phase::GAS)
                ? well_model_.phaseUsage().canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)
                : well_model_.phaseUsage().canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            std::map<std::string, Scalar> leaf_rate;
            for (const auto& well : well_model_) {
                if (!well->isInjector() || !well->wellEcl().predictionMode()) {
                    continue;
                }
                if (details::domainForWell(*well) != network.domain) {
                    continue;
                }
                const auto& node = well->wellEcl().groupName();
                if (!network.network.get().has_node(node)) {
                    continue;
                }
                const auto& ws = well_model_.wellState().well(well->indexOfWell());
                const Scalar current = ws.surface_rates[phase_pos];
                Scalar rate = current;
                if (const auto q = this->proxyInjectionRate(*well, phase_pos, deferred_logger)) {
                    // A well not on THP control is held by its group or rate target: it
                    // follows the node pressure only once the THP limit bites.
                    rate = (ws.injection_cmode == Well::InjectorCMode::THP)
                        ? *q : std::min(*q, current);
                }
                leaf_rate[node] += rate * well->wellEcl().getEfficiencyFactor(/*network=*/true);
            }
            for (const auto& [node, rate] : leaf_rate) {
                auto rates = group_state.has_network_leaf_node_injection_rates(node, *phase)
                    ? group_state.network_leaf_node_injection_rates(node, *phase)
                    : std::vector<Scalar>(well_model_.numPhases(), 0.0);
                rates[phase_pos] = comm.sum(rate);
                group_state.update_network_leaf_node_injection_rates(node, *phase, rates);
            }
        }
    }
    return imbalance;
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

    // With the simultaneous solve owning the choke nodes, the group thp is
    // already in the group state and the search below would fight it. Read
    // the parameters, not the flags: this runs before the flags are set on
    // the first pass of a step, and one pass of the search is enough to
    // register the group with a pressure nothing else will overwrite.
    if (well_model_.param().network_solver_ == "newton"
        && well_model_.param().network_autochoke_) {
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
                auto target = well_model_.groupStateHelper().
                    getAutoChokeGroupProductionTargetRate(group,
                                                          parentGroup,
                                                          resv_coeff,
                                                          efficiencyFactor);
                target_tmp = target.first;
                cmode_tmp = target.second;
            }
            using TargetCalculatorType =  GroupStateHelpers::TargetCalculator<Scalar, IndexTraits>;
            // Built on the control decided on above, not the group state's:
            // the state says NONE for a group under its limit, and the
            // calculator asserts on NONE.
            TargetCalculatorType tcalc{well_model_.groupStateHelper(), resv_coeff, cmode_tmp};
            if (!fld_none)
            {
                // Target is set for the autochoke group itself. Read it off the
                // deck control decided on above -- the group *state* may say NONE
                // when the group is under its limit, and asking the state for a
                // target then throws (NETWORK_MODEL5_STDW_AUTOCHK, day 3.2). A
                // group under its limit is simply a choke that ends up open.
                switch (cmode_tmp) {
                case Group::ProductionCMode::ORAT: target_tmp = ctrl.oil_target;    break;
                case Group::ProductionCMode::WRAT: target_tmp = ctrl.water_target;  break;
                case Group::ProductionCMode::GRAT: target_tmp = ctrl.gas_target;    break;
                case Group::ProductionCMode::LRAT: target_tmp = ctrl.liquid_target; break;
                case Group::ProductionCMode::RESV: target_tmp = ctrl.resv_target;   break;
                default:
                    target_tmp = well_model_.groupStateHelper().getProductionGroupTarget(group);
                }
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
            // Empty on the first pass of a run; dereferencing end() here
            // handed the group a garbage pressure.
            const Scalar nodal_pressure = (it != this->node_pressures_.end())
                ? it->second
                : network.node(nodeName).terminal_pressure().value_or(Scalar{0});
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
                                                                     local_deferredLogger,
                                                                     well_model_.param().network_autochoke_bracket_samples_);
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
                                               local_deferredLogger,
                                               well_model_.param().network_autochoke_bracket_samples_);

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
