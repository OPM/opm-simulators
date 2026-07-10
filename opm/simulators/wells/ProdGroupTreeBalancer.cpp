/*
  Copyright 2026 SINTEF Digital

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
#include <opm/simulators/wells/ProdGroupTreeBalancer.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>

namespace Opm::ProdGroupTreeBalancer {

namespace {

// ---------------------------------------------------------------------------
// Internal constants
// ---------------------------------------------------------------------------

/// Tolerance used for capping rate-vs-limit comparisons to avoid floating-point noise.
template<class Scalar>
constexpr Scalar kFeasibilityTolerance = Scalar(1e-10);

/// OIL / WATER / GAS canonical indices (always 0/1/2 in the tree).
constexpr int kOil   = 0;
constexpr int kWater = 1;
constexpr int kGas   = 2;

// ---------------------------------------------------------------------------
// Mode conversion helpers
// ---------------------------------------------------------------------------

/// Convert Group::ProductionCMode to Well::ProducerCMode for balancing.
/// FLD and VREP are not valid for wells, so they map to CMODE_UNDEFINED.
inline Well::ProducerCMode groupModeToWellMode(Group::ProductionCMode groupMode)
{
    switch (groupMode) {
        case Group::ProductionCMode::ORAT: return Well::ProducerCMode::ORAT;
        case Group::ProductionCMode::WRAT: return Well::ProducerCMode::WRAT;
        case Group::ProductionCMode::GRAT: return Well::ProducerCMode::GRAT;
        case Group::ProductionCMode::LRAT: return Well::ProducerCMode::LRAT;
        case Group::ProductionCMode::RESV: return Well::ProducerCMode::RESV;
        case Group::ProductionCMode::NONE: return Well::ProducerCMode::NONE;
        default: return Well::ProducerCMode::CMODE_UNDEFINED;  // FLD, VREP, etc.
    }
}

/// Convert Well::ProducerCMode to Group::ProductionCMode (for well's preferred mode).
/// BHP, THP, GRUP, CRAT are well-specific and map to NONE for group representation.
inline Group::ProductionCMode wellModeToGroupMode(Well::ProducerCMode wellMode)
{
    switch (wellMode) {
        case Well::ProducerCMode::ORAT: return Group::ProductionCMode::ORAT;
        case Well::ProducerCMode::WRAT: return Group::ProductionCMode::WRAT;
        case Well::ProducerCMode::GRAT: return Group::ProductionCMode::GRAT;
        case Well::ProducerCMode::LRAT: return Group::ProductionCMode::LRAT;
        case Well::ProducerCMode::RESV: return Group::ProductionCMode::RESV;
        case Well::ProducerCMode::NONE: return Group::ProductionCMode::NONE;
        default: return Group::ProductionCMode::NONE;  // BHP, THP, GRUP, etc.
    }
}

// ---------------------------------------------------------------------------
// Phase mapping helpers
// ---------------------------------------------------------------------------

/// Return the active phase index for canonical index \p canonical using pu,
/// or -1 if the phase is not active.
template<typename IndexTraits>
int activeIdx(const PhaseUsageInfo<IndexTraits>& pu, int canonical)
{
    if (canonical == kOil)   return pu.phaseIsActive(IndexTraits::oilPhaseIdx)   ? pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)   : -1;
    if (canonical == kWater) return pu.phaseIsActive(IndexTraits::waterPhaseIdx) ? pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx) : -1;
    if (canonical == kGas)   return pu.phaseIsActive(IndexTraits::gasPhaseIdx)   ? pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)   : -1;
    return -1;
}

/// Copy surface_rates (active-phase vector) into a canonical [oil,water,gas] array,
/// placing zero for inactive phases.  The sign convention is preserved
/// (negative = production).
template<class Scalar, typename IndexTraits>
std::array<Scalar, 3> toCanonical3(const std::vector<Scalar>& activeRates,
                                   const PhaseUsageInfo<IndexTraits>& pu)
{
    std::array<Scalar, 3> r{};
    for (int c = 0; c < 3; ++c) {
        const int a = activeIdx(pu, c);
        if (a >= 0 && a < static_cast<int>(activeRates.size())) {
            r[c] = activeRates[a];
        }
    }
    return r;
}

/// Convert canonical [oil,water,gas] array back into an active-phase vector.
template<class Scalar, typename IndexTraits>
std::vector<Scalar> toActive(const std::array<Scalar, 3>& canonical3,
                              const PhaseUsageInfo<IndexTraits>& pu)
{
    std::vector<Scalar> v(pu.numPhases, Scalar(0));
    for (int c = 0; c < 3; ++c) {
        const int a = activeIdx(pu, c);
        if (a >= 0) {
            v[a] = canonical3[c];
        }
    }
    return v;
}

// ---------------------------------------------------------------------------
// Rate projections onto individual limit types
// ---------------------------------------------------------------------------

/// Project vector \p v onto the scalar axis defined by \p cmode.
/// No sign change is performed; callers are responsible for applying the
/// correct sign to match their input convention.
///
/// Canonical rates are stored as negative = production, so for those arrays
/// the production rate is  -projectOnMode(rates, mode, resvCoeff).
/// Accumulated sums (rateSums, guideRateSums) are stored positive = production,
/// so for those arrays the production rate is  projectOnMode(sums, mode, resvCoeff).
template<class Scalar>
Scalar projectOnMode(const std::array<Scalar, 3>& v,
                     Well::ProducerCMode cmode,
                     const std::array<Scalar, 3>& resvCoeff)
{
    switch (cmode) {
    case Well::ProducerCMode::ORAT:  return v[kOil];
    case Well::ProducerCMode::WRAT:  return v[kWater];
    case Well::ProducerCMode::GRAT:  return v[kGas];
    case Well::ProducerCMode::LRAT:  return v[kOil] + v[kWater];
    case Well::ProducerCMode::BHP:
    case Well::ProducerCMode::THP:   return v[kOil] + v[kWater] + v[kGas];
    case Well::ProducerCMode::RESV:
    {
        Scalar r = Scalar(0);
        for (int c = 0; c < 3; ++c) r += v[c] * resvCoeff[c];
        return r;
    }
    default:
        return Scalar(0);
    }
}

inline int canonicalRateIndexForMode(Well::ProducerCMode mode)
{
    switch (mode) {
        case Well::ProducerCMode::ORAT: return kOil;
        case Well::ProducerCMode::WRAT: return kWater;
        case Well::ProducerCMode::GRAT: return kGas;
        default: return -1;
    }
}

// ---------------------------------------------------------------------------
// Guide-rate helpers
// ---------------------------------------------------------------------------

/// Convert Group::ProductionCMode to the corresponding GuideRateModel::Target.
/// Mirrors GroupStateHelper::getProductionGuideTargetModeFromControlMode.
GuideRateModel::Target productionCModeToGuideTarget(Group::ProductionCMode cmode)
{
    switch (cmode) {
    case Group::ProductionCMode::ORAT: return GuideRateModel::Target::OIL;
    case Group::ProductionCMode::WRAT: return GuideRateModel::Target::WAT;
    case Group::ProductionCMode::GRAT: return GuideRateModel::Target::GAS;
    case Group::ProductionCMode::LRAT: return GuideRateModel::Target::LIQ;
    case Group::ProductionCMode::RESV: return GuideRateModel::Target::RES;
    default:                           return GuideRateModel::Target::NONE;
    }
}

/// Look up the guide rate for node \p name in mode \p ctrlMode.
template<class Scalar>
Scalar getGuideRateForMode(const std::string& name,
                           const std::array<Scalar, 3>& rates,
                           Group::ProductionCMode ctrlMode,
                           const GuideRate& guideRate)
{
    const auto target = productionCModeToGuideTarget(ctrlMode);
    if (target == GuideRateModel::Target::NONE) return -1.0;
    // GuideRate::RateVector is {oil, gas, water} (not canonical [oil,water,gas])
    const Scalar oilRate   = -rates[kOil];
    const Scalar gasRate   = -rates[kGas];
    const Scalar waterRate = -rates[kWater];
    const GuideRate::RateVector rv{oilRate, gasRate, waterRate};
    const Scalar gr = guideRate.get(name, target, rv);
    return gr >= Scalar(0) ? gr : 0.0;
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void populateWellNode(ProdGroupTreeNode<Scalar>& node,
                      const std::string& wellName,
                      const Schedule& schedule,
                      const WellState<Scalar, IndexTraits>& wellState,
                      const GroupState<Scalar>& /*groupState*/,
                      const GuideRate& /*guideRate*/,
                      const SummaryState& /*summaryState*/,
                      int reportStep,
                      int fipnum,
                      int pvtreg,
                      const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                      const std::unordered_map<std::string, std::pair<int, Scalar>>& limits)
{
    // Use globally available state so this function produces the same result on
    // every MPI rank regardless of which rank owns the well.
    //
    // - well_rates / currentWellRates: populated for ALL wells on ALL ranks.
    // - isProductionGrup(): from GlobalWellInfo, communicated via updateGlobalIsGrup().
    // - getGlobalEfficiencyScalingFactor(): from GlobalWellInfo (comm.min reduction).
    // - limits: passed in from prepareWellsForBalancing_*() — globally consistent.

    const auto& eclWell = schedule.getWell(wellName, reportStep);
    const auto& pu      = wellModel.phaseUsage();

    node.name   = wellName;
    node.type   = ProdNodeType::Well;
    node.parent = eclWell.groupName();
    node.availableForGroupControl = eclWell.isAvailableForGroupControl();
    node.hasGuideRate = true;
    node.efficiencyFactor = eclWell.getEfficiencyFactor()
                          * wellState.getGlobalEfficiencyScalingFactor(wellName);

    // well_rates/currentWellRates stores positive = production; the balancer uses negative = production.
    // Convert via the same active→canonical mapping as toCanonical3.
    {
        const auto& wr = wellState.currentWellRates(wellName);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.rates[c] = (a >= 0 && a < static_cast<int>(wr.size())) ? -wr[a] : Scalar(0);
        }
    }
    node.initialRates = node.rates; // snapshot of current rates

    // RESV conversion coefficients — computed from globally consistent rates.
    {
        std::vector<Scalar> coeffVec(pu.numPhases, Scalar(0));
        std::vector<Scalar> posRates(pu.numPhases, Scalar(0));
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            if (a >= 0) posRates[a] = -node.rates[c]; // node.rates is negative = production
        }
        wellModel.calcResvCoeff(fipnum, pvtreg, posRates, coeffVec);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.resvCoeff[c] = (a >= 0) ? coeffVec[a] : Scalar(0);
        }
    }

    // Populate limits and control mode from the globally consistent limits map.
    // Single lookup used for both node.mode and node.Limits below.
    const auto limIt = limits.find(wellName);
    if (wellState.isProductionGrup(wellName)) {
        // We don't know actual ORAT/WRAT/GRAT/RESV now, so set to GRUP
        // Actual mode is set by parent when we traverse at end of buildTree().  
        node.mode         = Well::ProducerCMode::GRUP; 
        node.modeCategory = ProdNodeModeCategory::Group;
    } else {
        // Individually-controlled well: mode equals the strictest limit mode.
        node.modeCategory = ProdNodeModeCategory::Individual;
        node.mode = (limIt != limits.end())
            ? static_cast<Well::ProducerCMode>(limIt->second.first)
            : Well::ProducerCMode::CMODE_UNDEFINED;
    }

    if (limIt == limits.end()) return;

    const auto cachedMode  = static_cast<Well::ProducerCMode>(limIt->second.first);
    const Scalar limitValue = limIt->second.second;

    if (cachedMode == Well::ProducerCMode::CMODE_UNDEFINED) return;

    if (cachedMode == Well::ProducerCMode::THP && limitValue <= Scalar(0)) {
        // No stable operating point at THP limit: well cannot produce.
        assert(false);
        node.rates = {};
        node.initialRates = {};
        return;
    }

    node.Limits[cachedMode] = limitValue;
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
void populateGroupNode(ProdGroupTreeNode<Scalar>& node,
                       const std::string& groupName,
                       const Schedule& schedule,
                       const WellState<Scalar, IndexTraits>& /*wellState*/,
                       const GroupState<Scalar>& groupState,
                       const PhaseUsageInfo<IndexTraits>& pu,
                       const SummaryState& summaryState,
                       int reportStep,
                       int fipnum,
                       int pvtreg,
                       const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel)
{
    const auto& group = schedule.getGroup(groupName, reportStep);

    node.name   = groupName;
    node.type   = ProdNodeType::Group;
    node.parent = (groupName == "FIELD") ? "" : group.parent();

    // FIELD group is never available for group control (has no parent to control it)
    if (groupName == "FIELD") {
        node.availableForGroupControl = false;
        node.efficiencyFactor = Scalar(1);
    } else {
        node.availableForGroupControl = group.productionGroupControlAvailable();
        node.efficiencyFactor = group.getGroupEfficiencyFactor();
    }

    // Children
    node.children = group.groups();
    for (const auto& w : group.wells()) {
        node.children.push_back(w);
    }
    // If no children, group is either a satellite group with fixed rates from GSATPROD, 
    // or inactive group. In the latter case, treat it as a satellite group with zero rates (hack).
    if (node.children.empty()) {
        node.isSatellite = true; 
        node.availableForGroupControl = false;
        node.modeCategory = ProdNodeModeCategory::Individual;
        if (group.hasSatelliteProduction()) {
            const auto& gsat_prod = schedule[reportStep].gsatprod();
            if (gsat_prod.has(groupName)) {
                const auto& sat_rates = gsat_prod.get(groupName, summaryState);
                // GSATPROD rates map: rate[phase_enum_value] where phase is Phase::OIL, WATER, GAS
                node.rates[0] = -sat_rates.rate[static_cast<int>(Phase::OIL)];
                node.rates[1] = -sat_rates.rate[static_cast<int>(Phase::WATER)];
                node.rates[2] = -sat_rates.rate[static_cast<int>(Phase::GAS)];
            }
        }
    } else if (groupState.has_production_rates(groupName)) {
        const auto& gr = groupState.production_rates(groupName);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.rates[c] = (a >= 0 && a < static_cast<int>(gr.size())) ? -gr[a] : Scalar(0);
        }
    }
    node.initialRates = node.rates; // snapshot of beginning-of-timestep rates

    // Group control status and own control mode
    const auto ctrl = groupState.has_production_control(groupName)
        ? groupState.production_control(groupName)
        : Group::ProductionCMode::NONE;

    node.mode = groupModeToWellMode(ctrl); // Convert to well mode for balancing (FLD/VREP→CMODE_UNDEFINED)

    // RESV conversion coefficients for this group (should check if needed)
    {
        std::vector<Scalar> coeffVec(pu.numPhases, Scalar(0));
        std::vector<Scalar> posRates(pu.numPhases, Scalar(0));
        // Convert group rates (negative=production) to positive rates for RESV calculation
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            if (a >= 0) {
                posRates[a] = -node.rates[c];
            }
        }
        wellModel.calcResvCoeff(fipnum, pvtreg, posRates, coeffVec);
        for (int c = 0; c < 3; ++c) {
            const int a = activeIdx(pu, c);
            node.resvCoeff[c] = (a >= 0) ? coeffVec[a] : Scalar(0);
        }
    }

    // Populate individual limits from group production controls
    if (group.isProductionGroup()) {
        const auto controls = group.productionControls(summaryState);

        // Set preferredMode from schedule control mode (store Group::ProductionCMode directly)
        // FLD and NONE will be resolved via inheritance later
        node.preferredMode = controls.cmode;

        // A group participates in guide-rate balancing if item 10 is present or if
        // it has active individual limits that keep it non-transparent.
        node.hasGuideRate = (controls.guide_rate_def != Group::GuideRateProdTarget::NO_GUIDE_RATE);
        
        const auto& action = controls.group_limit_action; 
        const bool actionAllIsRate = (action.allRates == Opm::Group::ExceedAction::RATE);
        // Populate individual limits from schedule controls
        if (group.has_control(Group::ProductionCMode::ORAT)) {
            if (actionAllIsRate || action.oil == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::ORAT] = controls.oil_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::WRAT)) {
            if (actionAllIsRate || action.water == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::WRAT] = controls.water_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::GRAT)) {
            if (actionAllIsRate || action.gas == Opm::Group::ExceedAction::RATE) {
                node.Limits[Well::ProducerCMode::GRAT] = controls.gas_target;
            }
        }
        if (group.has_control(Group::ProductionCMode::LRAT)) {
            // Skip degenerate LRAT == ORAT case (no water production): same guard as
            // GroupStateHelper::checkGroupProductionConstraints.
            if ((actionAllIsRate || action.liquid == Opm::Group::ExceedAction::RATE)
                && controls.liquid_target != controls.oil_target)
            {
                node.Limits[Well::ProducerCMode::LRAT] = controls.liquid_target;
            }
        }
        // GRAT: prefer GCONSALE-adjusted target when present.
        // Mirrors GroupStateHelper::getProductionGroupTargetForMode_(GRAT).
        if (node.Limits.count(Well::ProducerCMode::GRAT) > 0) {
            if (groupState.has_grat_sales_target(groupName)) {
                const Scalar salesTarget = groupState.grat_sales_target(groupName);
                if (salesTarget > Scalar(0)) {
                    node.Limits[Well::ProducerCMode::GRAT] = salesTarget;
                }
            }
        }
        // RESV: mirror GroupStateHelper's GPMAINT-override logic.
        // Mirrors GroupStateHelper::getProductionGroupTargetForMode_(RESV).
        if (group.has_control(Group::ProductionCMode::RESV)) {
            Scalar resv_target = Scalar(0);
            if (group.has_gpmaint_control(Group::ProductionCMode::RESV)
                && groupState.has_gpmaint_target(groupName))
            {
                resv_target = groupState.gpmaint_target(groupName);
            } else {
                resv_target = controls.resv_target;
            }
            if (resv_target > Scalar(0)) {
                node.Limits[Well::ProducerCMode::RESV] = resv_target;
            }
        }
    }
}

/// Top-down pass over the tree: propagate the effective group control mode,
/// assign guide rates, resolve preferred control modes, and mark nodes that
/// have a limited ancestor.  Recurses depth-first from \p nodeName.
///
/// \param[in,out] tree                   The production group tree.
/// \param[in]     guideRate              Guide-rate data.
/// \param[in]     nodeName               Node to process.
/// \param[in]     inheritedMode          Control mode inherited from the parent.
/// \param[in]     inheritedPreferredMode Preferred mode inherited from the parent.
/// \param[in]     parentSeesLimits       True if any ancestor carries a rate limit.
template<class Scalar>
void propagateGuideRatesAndMode(Tree<Scalar>& tree,
                                const GuideRate& guideRate,
                                const std::string& nodeName,
                                Well::ProducerCMode inheritedMode,
                                Group::ProductionCMode inheritedPreferredMode,
                                bool parentSeesLimits)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);
    if (node.isSatellite) return;

    node.hasLimitedAncestor = parentSeesLimits;
    if (!node.availableForGroupControl || !node.hasLimitedAncestor) {
        // Not available for group control: no inheritance from above.
        // If preferred mode is FLD/NONE and limits exist, revert to individual.
        if (node.preferredMode == Group::ProductionCMode::FLD ||
            node.preferredMode == Group::ProductionCMode::NONE) {
            if (!node.Limits.empty()) {
                node.mode = node.Limits.begin()->first;
            }
        }
    } else {
        // Available for group control: inherit preferred mode if unset.
        if (node.preferredMode == Group::ProductionCMode::FLD ||
            node.preferredMode == Group::ProductionCMode::NONE) {
            node.preferredMode = inheritedPreferredMode;
        }
        if (node.modeCategory != ProdNodeModeCategory::Individual) {
            node.mode = inheritedMode;
        }
    }

    // Assign guide rates using the resolved control mode.
    if (node.hasGuideRate) {
        const bool validMode = (node.mode != Well::ProducerCMode::CMODE_UNDEFINED) &&
                               (node.mode != Well::ProducerCMode::NONE);
        const auto groupMode = validMode ? wellModeToGroupMode(node.mode) : node.preferredMode;
        node.groupTarget.ctrlMode  = groupMode;
        node.groupTarget.guideRate = getGuideRateForMode(nodeName, node.initialRates, groupMode, guideRate);

        // Also fill groupTargetFallback using the preferred mode.
        node.groupTargetFallback.ctrlMode  = node.preferredMode;
        node.groupTargetFallback.guideRate = getGuideRateForMode(
            nodeName, node.initialRates, node.preferredMode, guideRate);

        // If the preferred-phase guide rate is zero the node cannot participate
        // in guide-rate balancing (would make the problem infeasible); exclude it.
        if (node.groupTargetFallback.guideRate <= Scalar(0)) {
            node.availableForGroupControl = false;
            node.modeCategory = ProdNodeModeCategory::Individual;
            if (!node.Limits.empty()) {
                node.mode = node.Limits.begin()->first;
            }
        }
    }

    // Recurse into children.
    for (const auto& childName : node.children) {
        propagateGuideRatesAndMode(tree, guideRate, childName,
                                   node.mode, node.preferredMode,
                                   node.hasLimitedAncestor || !node.Limits.empty());
    }
}

} // anonymous namespace

// ===========================================================================
// Public API
// ===========================================================================

template<class Scalar, typename IndexTraits>
Tree<Scalar> buildTree(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                       const SummaryState& summaryState,
                       int reportStep,
                       const std::unordered_map<std::string, std::pair<int, Scalar>>& limits)
{
    const auto& schedule   = wellModel.schedule();
    const auto& wellState  = wellModel.wellState();
    const auto& groupState = wellModel.groupState();
    const auto& guideRate  = wellModel.guideRate();
    const auto& pu         = wellModel.phaseUsage();

    // Get RESV conversion parameters (fipnum, pvtreg)
    const auto [fipnum, pvtreg] = wellModel.getGroupFipnumAndPvtreg();

    Tree<Scalar> tree;

    // A well is valid for the tree iff it appears in the limits map — which is
    // globally consistent on all ranks after gatherWellLimits_().  This covers:
    //   • non-local wells: limits is from comm.sum() so every rank has the same data
    //   • stopped wells: prepareWellsForBalancing_*() never adds them to limits
    //   • injectors: filtered out before building localLimits
    auto validWellNode = [&](const std::string& wname) -> bool {
        return limits.count(wname) > 0;
    };

    // Use a stack-based traversal starting from FIELD
    std::vector<std::string> toVisit = {"FIELD"};
    while (!toVisit.empty()) {
        const std::string name = toVisit.back();
        toVisit.pop_back();

        if (schedule.hasWell(name, reportStep)) {
            // Well node: visible to all ranks via currentWellRates().
            if (validWellNode(name)) {
                auto& node = tree[name];
                populateWellNode(node, name, schedule, wellState, groupState,
                                 guideRate, summaryState, reportStep,
                                 fipnum, pvtreg, wellModel, limits);
            }
        } else {
            // Group node: populate and push children onto the stack
            const auto& group = schedule.getGroup(name, reportStep);
            auto& node = tree[name];
            populateGroupNode(node, name, schedule, wellState, groupState, pu,
                            summaryState, reportStep, fipnum, pvtreg, wellModel);

            // Push children onto the stack
            for (const auto& child : group.groups()) {
                toVisit.push_back(child);
            }
            for (const auto& well : group.wells()) {
                toVisit.push_back(well);
            }
        }
    }
    
    // Top-down pass: propagate effective group control modes, compute guide rates,
    // propagate preferred control, and detect transparent groups.
    propagateGuideRatesAndMode(tree, guideRate, "FIELD",
                               Well::ProducerCMode::CMODE_UNDEFINED,
                               Group::ProductionCMode::NONE,
                               /* parentSeesLimits */ false);
    return tree;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string> getSubTreeOrdering(const Tree<Scalar>& tree,
                                            const std::string& rootName)
{
    // Collect all rootnodes of subtrees that can be balanced individually. This 
    // includes:
    // - top-most node with a given limit (often the FIELD group, but can be a subgroup if FIELD has no limits)
    // - groups that are not available for group control (e.g., have GCONPROD item 8 = NO)
    // - wells that are not available for group control (e.g., have WGRUPCON item 2 = NO)
    // - satellite groups are skipped since they have fixed rates (requires no balancing)
    // Post-order traversal: leaves first, root last.
    std::vector<std::string> ordering;
    std::vector<std::string> stack;
    stack.push_back(rootName);

    // Safety limit against malformed trees (e.g. cycles); large enough for any
    // realistic FIELD hierarchy.
    constexpr int kMaxIterations = 10000;
    int iterCount = 0;
    while (!stack.empty()) {
        if (++iterCount > kMaxIterations) {
            break;
        }
        const std::string name = stack.back();
        stack.pop_back();
        if (tree.count(name) == 0) continue;
        const auto& node = tree.at(name);
        if (node.isSatellite) continue;
        // has limit and cannot be controlled by parent
        if (!node.Limits.empty() &&
            (!node.hasLimitedAncestor || !node.availableForGroupControl)) {
            ordering.push_back(name);
        }
        // Push children before the current node so they appear first in output
        for (const auto& child : node.children) {
            if (tree.count(child) > 0) {
                stack.push_back(child);
            }
        }
    }
    // Reverse so children appear before parents (post-order)
    std::ranges::reverse(ordering);
    return ordering;
}

// ---------------------------------------------------------------------------
// Sorting-based algorithm helpers
// ---------------------------------------------------------------------------

// Forward declaration: balanceGroupTree is defined later but called by
// updateTransparentGroups and distributeFallbackRates below.
template<class Scalar>
void balanceGroupTree(Tree<Scalar>& tree,
                      const std::string& nodeName,
                      const GuideRate& guideRate,
                      Well::ProducerCMode targetMode,
                      Scalar targetRate,
                      Scalar tol,
                      DeferredLogger& logger);

/// Compute the product of efficiency factors from \p fromName (inclusive) up to
/// \p toAncestorName (exclusive).  Mirrors GroupStateHelper::accumulateGroupEfficiencyFactor:
/// each group's own efficiency is multiplied in as the traversal moves up the tree,
/// so the result is the accumulated efficiency seen by the ancestor.
template<class Scalar>
Scalar accumulatedEfficiency(const Tree<Scalar>& tree,
                             const std::string& fromName,
                             const std::string& toAncestorName)
{
    Scalar eff = Scalar(1);
    std::string cur = fromName;
    while (cur != toAncestorName && tree.count(cur) > 0) {
        eff *= tree.at(cur).efficiencyFactor;
        cur = tree.at(cur).parent;
    }
    return eff;
}

template<class Scalar>
bool hasFreePath(const Tree<Scalar>& tree,
                 const std::string& nodeName,
                 const std::string& originName)
{
    // Check that path from node up to originName is "free" (all transparent)
    if (nodeName == originName) return true;
    if (tree.count(nodeName) == 0) return false;

    const auto& node = tree.at(nodeName);
    if (node.parent == originName) return true;

    if (node.parent.empty() || tree.count(node.parent) == 0) return false;

    const auto& parent = tree.at(node.parent);
    // !parent.visited: once a transparent group has been balanced by updateTransparentGroups
    // (which sets visited=true via categorizeBalancedNode), paths through it are no longer
    // free — preventing subsequent c-children from routing through an already-committed node.
    const bool isFreeStep = (parent.modeCategory == ProdNodeModeCategory::Transparent) && !parent.visited;

    return isFreeStep && hasFreePath(tree, node.parent, originName);
}

template<class Scalar>
void incrementParentRateSums(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::string& originName,
                             const std::optional<std::array<Scalar, 3>>& rates = std::nullopt)
{
    // Increment rate_sums from node up to origin
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    // Determine rates to add (convert to positive production for rateSums,
    // scaled by the node's own efficiency factor so the parent sees the
    // efficiency-adjusted contribution).
    std::array<Scalar, 3> ratesToAdd;
    if (rates.has_value()) {
        // Caller already provides the delta in parent-frame units.
        ratesToAdd = rates.value();
    } else {
        // node.rates are negative (production); convert to positive and scale.
        for (int c = 0; c < 3; ++c) {
            ratesToAdd[c] = node.efficiencyFactor * (-node.rates[c]);
        }
    }

    // Walk from nodeName up to originName.  ratesToAdd starts in nodeName's
    // direct-parent frame (it already includes nodeName.efficiencyFactor).
    // Each intermediate transparent group must apply its own efficiency factor
    // to convert from its native frame to its parent's frame before the value
    // is added to the grandparent's rateSums — mirroring the accumulated
    // efficiency convention used in GroupStateHelper.
    std::string cur = nodeName;
    while (cur != originName && tree.count(cur) > 0) {
        const auto& curNode = tree.at(cur);
        // For every group beyond the starting node, scale by that group's
        // efficiency factor to step into its parent's frame.
        if (cur != nodeName) {
            for (int c = 0; c < 3; ++c) {
                ratesToAdd[c] *= curNode.efficiencyFactor;
            }
        }
        const std::string next = curNode.parent;
        if (!next.empty() && tree.count(next) > 0) {
            auto& nextNode = tree.at(next);
            for (int c = 0; c < 3; ++c) {
                nextNode.rateSums[c] += ratesToAdd[c];
            }
        }
        cur = next;
    }
}

template<class Scalar>
void decrementParentGuideRateSums(Tree<Scalar>& tree,
                                   const std::string& nodeName,
                                   const std::string& originName,
                                   const std::array<Scalar, 3>& guideRateSums)
{
    // Subtract guide_rate_sums from node up to and including origin
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);

    std::string parent = node.parent;
    while (!parent.empty() && tree.count(parent) > 0) {
        auto& parentNode = tree.at(parent);
        for (int c = 0; c < 3; ++c) {
            parentNode.guideRateSums[c] -= guideRateSums[c];
        }
        if (parent == originName) break;  // Stop after updating origin
        parent = parentNode.parent;
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
getLocalTreeDescendants(const Tree<Scalar>& tree, const std::string& nodeName)
{
    // Collect local subtree descendants of nodeName, categorised as:
    // d: nodes with guide-rate available for group control
    // df: nodes not available for group control (fixed)
    // dt: nodes without guide-rate in between d and node (transparent)
    // TODO: can split dt into those with and without limits for improved
    // efficiency

    std::vector<std::string> d, df, dt;

    if (tree.count(nodeName) == 0) return {d, df, dt};
    const auto& node = tree.at(nodeName);

    for (const auto& childName : node.children) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        if (!child.availableForGroupControl) {
            df.push_back(childName);
        } else if (child.hasGuideRate) {
            d.push_back(childName);
        } else {
            // Recursive descent for transparent groups
            auto [d_child, df_child, dt_child] = getLocalTreeDescendants(tree, childName);
            d.insert(d.end(), d_child.begin(), d_child.end());
            df.insert(df.end(), df_child.begin(), df_child.end());
            dt.insert(dt.end(), dt_child.begin(), dt_child.end());
            dt.push_back(childName);
        }
    }

    return {d, df, dt};
}

// ---------------------------------------------------------------------------
template<class Scalar>
void updateGuideRatesForMode(Tree<Scalar>& tree,
                             const std::vector<std::string>& c,
                             Well::ProducerCMode mode,
                             const GuideRate& guideRate)
{
    const auto modeAsGroupMode = wellModeToGroupMode(mode);
    // Update guide rates for children if needed (e.g., if mode changed)
    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        if (child.groupTarget.ctrlMode != modeAsGroupMode) {
            child.groupTarget.ctrlMode = modeAsGroupMode;
            child.groupTarget.guideRate = getGuideRateForMode(childName, child.initialRates, modeAsGroupMode, guideRate);
        }
    }
}

template<class Scalar>
void resetRatesAndGuideRateSums(Tree<Scalar>& tree,
                                const std::vector<std::string>& c,
                                const std::vector<std::string>& c_fixed,
                                const std::vector<std::string>& c_trans,
                                const std::string& origin,
                                Well::ProducerCMode mode)
{
    // Reset guide-rate-sums and rate-sums for all relevant nodes
    // rateSums are set to zero for all nodes in c, c_fixed, c_trans, and origin
    // guideRateSums:
    //  c: 1x3 array with guideRateSums(mode) = guiderate, and 
    //     proportional to c.rates
    //  c_fixed: just zeros
    //  c_trans/origin: sum of guideRateSums from children 

    std::vector<std::string> allNodes = {origin};
    allNodes.insert(allNodes.end(), c.begin(), c.end());
    allNodes.insert(allNodes.end(), c_fixed.begin(), c_fixed.end());
    allNodes.insert(allNodes.end(), c_trans.begin(), c_trans.end());

    for (const auto& nodeName : allNodes) {
        if (tree.count(nodeName) == 0) continue;
        auto& node = tree.at(nodeName);
        node.rateSums = {Scalar(0), Scalar(0), Scalar(0)};
        node.guideRateSums = {Scalar(0), Scalar(0), Scalar(0)};
        node.visited = false;
    }

    // Set transparent nodes category
    for (const auto& nodeName : c_trans) {
        if (tree.count(nodeName) > 0) {
            tree.at(nodeName).modeCategory = ProdNodeModeCategory::Transparent;
        }
    }

    if (tree.count(origin) == 0) return;
    const auto& originNode = tree.at(origin);
    const Group::ProductionCMode modePreferred = originNode.preferredMode;
    // Convert mode to Group::ProductionCMode for comparison
    const Group::ProductionCMode modeAsGroupMode = wellModeToGroupMode(mode);
    const bool modeIsPref = (modeAsGroupMode == modePreferred);

    // Compute total guide rate for mode if not preferred
    Scalar guideSumOrigin = Scalar(0);
    if (!modeIsPref) {
        for (const auto& childName : c) {
            if (tree.count(childName) > 0) {
                const auto& child = tree.at(childName);
                guideSumOrigin += child.groupTarget.guideRate;
            }
        }
    }

    // Set guide_rate_sums for each child in c
    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        auto& child = tree.at(childName);
        child.useFallback = false;

        const auto& rates = child.rates;

        // Check for fallback condition (non-preferred mode with tiny fraction)
        if (!modeIsPref) {
            const Scalar rateSum = std::accumulate(rates.begin(), rates.end(), Scalar(0),
                                                   [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar rateForCurrentMode = -projectOnMode(rates, mode, child.resvCoeff);
            const Scalar rateFracMode = (rateSum > Scalar(0)) ? (rateForCurrentMode / rateSum) : Scalar(0);
            const Scalar guideRateForMode = child.groupTarget.guideRate;
            const Scalar guideRatioMode = (guideSumOrigin > Scalar(0))
                ? (guideRateForMode / guideSumOrigin) : Scalar(0);

            if (rateFracMode < std::sqrt(std::numeric_limits<Scalar>::epsilon()) ||
                guideRatioMode < std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
                child.useFallback = true;
                child.guideRateSums = {Scalar(0), Scalar(0), Scalar(0)};
                continue;
            }
        }

        // Set guide-rate sums proportional to rates, scaled by the child's
        // efficiency factor (hence guide_rate_sums carries both the guiderate
        // and the linear term). Parent-frame guide-rate sums and rate-sums are
        // consistently in the same (parent) frame.
        const Scalar guideRateValue = child.groupTarget.guideRate;
        const Scalar rateMagnitude = -projectOnMode(rates, mode, child.resvCoeff);

        // Scale by accumulated efficiency from this child up to the origin,
        // mirroring FractionCalculator::guideRate() which multiplies by
        // group.getGroupEfficiencyFactor() at each level.  This ensures gk/gsum
        // fractions are correct when promoted wells pass through transparent
        // groups with efficiency factors different from 1.
        const Scalar accEff = accumulatedEfficiency(tree, childName, origin);
        if (rateMagnitude > Scalar(0)) {
            for (int cmp = 0; cmp < 3; ++cmp) {
                child.guideRateSums[cmp] =
                    accEff * (guideRateValue / rateMagnitude) * (-rates[cmp]);
            }
        } else {
            child.guideRateSums = {Scalar(0), Scalar(0), Scalar(0)};
        }

        // Propagate up to origin
        std::string parent = childName;
        while (parent != origin && tree.count(parent) > 0) {
            const auto& pnode = tree.at(parent);
            parent = pnode.parent;
            if (tree.count(parent) > 0) {
                auto& parentNode = tree.at(parent);
                for (int cmp = 0; cmp < 3; ++cmp) {
                    parentNode.guideRateSums[cmp] += child.guideRateSums[cmp];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<Scalar>
computeRatiosForSorting(const Tree<Scalar>& tree, const std::vector<std::string>& c,
                        const std::string& origin)
{
    // Get limit-to-guide-rate ratios for distribution ordering.
    // For each child, production rates evolve linearly as
    //   rate(a) = rateSums + alpha * guideRateSums,
    // and we find the smallest non-negative alpha that activates any limit.
    // All quantities (rateSums, guideRateSums, effectiveLimit) are expressed in
    // the origin's frame using the accumulated efficiency factor.
    const size_t nc = c.size();
    std::vector<Scalar> ratios(nc, std::numeric_limits<Scalar>::infinity());

    for (size_t k = 0; k < nc; ++k) {
        const auto& childName = c[k];
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);

        if (child.Limits.empty()) {
            continue;  // ratios[k] already infinity
        }
        // For 'None' nodes, the node is already limited by its children.
        // Use current node rates (converted to parent frame via efficiencyFactor) as
        // per-phase limits.  rateSums and guideRateSums are already in parent frame, so
        // the limit must be scaled consistently.
        // Accumulated efficiency from this child to the origin: converts the
        // child's native-frame limit to the origin frame consistently with the
        // accEff-scaled guideRateSums set by resetRatesAndGuideRateSums.
        const Scalar accEff = accumulatedEfficiency(tree, childName, origin);

        if (child.modeCategory == ProdNodeModeCategory::None) {
            Scalar minAlpha = std::numeric_limits<Scalar>::infinity();
            for (int cmp = 0; cmp < 3; ++cmp) {
                const Scalar limit   = accEff * (-child.rates[cmp]);  // origin frame
                const Scalar current = child.rateSums[cmp];   // positive production at alpha = 0
                const Scalar slope   = child.guideRateSums[cmp];
                if (slope <= Scalar(0)) continue;
                minAlpha = std::min(minAlpha, (limit - current) / slope);
            }
            ratios[k] = std::max(minAlpha, Scalar(0));
            continue;
        }

        // rateSums and guideRateSums are in the origin frame (scaled by accEff).
        // Limits are native (child frame). Convert the limit to the origin frame.
        Scalar minAlpha = std::numeric_limits<Scalar>::infinity();
        for (const auto& [mode, limit] : child.Limits) {
            // Effective limit in the origin frame.
            const Scalar effectiveLimit = accEff * limit;
            const Scalar current = projectOnMode(child.rateSums, mode, child.resvCoeff);
            const Scalar slope   = projectOnMode(child.guideRateSums, mode, child.resvCoeff);

            // This mode cannot become more restrictive along the current distribution direction.
            if (slope <= Scalar(0)) continue;

            minAlpha = std::min(minAlpha, (effectiveLimit - current) / slope);
        }

        ratios[k] = std::max(minAlpha, Scalar(0));
    }
    return ratios;
}

// ---------------------------------------------------------------------------

template<class Scalar>
std::vector<std::string>
updateTransparentGroups(Tree<Scalar>& tree,
                        const std::string& nodeName,
                        std::vector<std::string>& c_trans,
                        const GuideRate& guideRate,
                        Scalar nextRatio,
                        Scalar qm,
                        Well::ProducerCMode mode,
                        Scalar tol,
                        DeferredLogger& logger)
{
    // Check if any transparent groups have limit-to-guide ratio less than nextRatio
    std::vector<std::string> c_trans_update;

    bool foundLessThanNext = true;
    while (foundLessThanNext && !c_trans.empty()) {
        const auto ratios = computeRatiosForSorting(tree, c_trans, nodeName);

        // Find minimum ratio
        Scalar ratioMin = std::numeric_limits<Scalar>::max();
        size_t minIdx = 0;
        bool foundMin = false;
        for (size_t i = 0; i < ratios.size(); ++i) {
            if (ratios[i] < ratioMin) {
                ratioMin = ratios[i];
                minIdx = i;
                foundMin = true;
            }
        }

        if (!foundMin || ratioMin >= nextRatio) {
            foundLessThanNext = false;
            break;
        }

        if (tree.count(nodeName) == 0) break;
        const auto& originNode = tree.at(nodeName);
        const Scalar gsum = projectOnMode(originNode.guideRateSums, mode, originNode.resvCoeff);
        const Scalar qmRemain = qm - projectOnMode(originNode.rateSums, mode, originNode.resvCoeff);

        if (ratioMin*gsum >= qmRemain) {
            // ratioMin < nextRatio but not enough remaining production to reach it
            foundLessThanNext = false;
            break;
        }

        const auto& ckName = c_trans[minIdx];
        if (!hasFreePath(tree, ckName, nodeName)) {
            c_trans.erase(c_trans.begin() + static_cast<int>(minIdx));
            continue;
        }

        if (tree.count(ckName) == 0) {
            c_trans.erase(c_trans.begin() + static_cast<int>(minIdx));
            continue;
        }

        c_trans_update.push_back(ckName);
        auto& ck = tree.at(ckName);

        // Save original values
        const auto rateSumsOrig = ck.rateSums;
        const auto guideRateSumsOrig = ck.guideRateSums;
        const Scalar gk = projectOnMode(ck.guideRateSums, mode, ck.resvCoeff);
        // Accumulated efficiency from ck to nodeName (origin).
        // rateSumsOrig is in ck's native frame (first-hop unscaled from
        // incrementParentRateSums); multiply by accEff to convert to origin frame
        // before adding to the guide-rate-proportional budget share.
        const Scalar accEff = accumulatedEfficiency(tree, ckName, nodeName);
        // qkParent = guide-rate-proportional share of remaining budget
        //           + existing allocation already accumulated in rateSumsOrig.
        const Scalar qkParent = (gsum > Scalar(0) ? (gk / gsum) * qmRemain : Scalar(0))
                              + accEff * projectOnMode(rateSumsOrig, mode, ck.resvCoeff);
        const Scalar qk = qkParent / accEff;

        // Recursive balance
        balanceGroupTree(tree, ckName, guideRate, mode, qk, tol, logger);

        // ratesDelta is in ck's direct-parent frame: ck.efficiencyFactor converts
        // ck's post-balance native rates to that frame, and rateSumsOrig is also in
        // ck's native frame, so the difference is consistent.  incrementParentRateSums
        // then propagates further up with the correct intermediate-eff scaling.
        std::array<Scalar, 3> ratesDelta;
        for (int c = 0; c < 3; ++c) {
            ratesDelta[c] = ck.efficiencyFactor * ((-ck.rates[c]) - rateSumsOrig[c]);
        }
        incrementParentRateSums(tree, ckName, nodeName, std::make_optional(ratesDelta));
        decrementParentGuideRateSums(tree, ckName, nodeName, guideRateSumsOrig);

        // Remove from c_trans
        c_trans.erase(c_trans.begin() + static_cast<std::ptrdiff_t>(minIdx));
    }

    return c_trans_update;
}

// ---------------------------------------------------------------------------

template<class Scalar>
void distributeFallbackRates(Tree<Scalar>& tree,
                             const std::string& nodeName,
                             const std::vector<std::string>& c,
                             const GuideRate& guideRate,
                             Well::ProducerCMode mode,
                             Scalar tol,
                             DeferredLogger& logger)
{
    // Distribute rates to children that need fallback (tiny fractions in non-preferred mode)
    if (tree.count(nodeName) == 0) return;
    const auto& node = tree.at(nodeName);
    const Group::ProductionCMode modePreferred = node.preferredMode;
    const int modeIdx = canonicalRateIndexForMode(mode);
    const Well::ProducerCMode modePreferredAsWell = groupModeToWellMode(modePreferred);

    // Collect fallback children and compute the preferred-mode rate/guide ratio
    // from already-balanced GRUP children.
    std::vector<std::string> cFallback;
    std::vector<Scalar> cFallbackGuideRates;
    Scalar grupRateSum = Scalar(0);
    Scalar grupGuideRateSum = Scalar(0);

    for (const auto& childName : c) {
        if (tree.count(childName) == 0) continue;
        const auto& child = tree.at(childName);
        const Scalar preferredGuideRate = getGuideRateForMode(childName, child.initialRates, modePreferred, guideRate);

        if (child.useFallback) {
            cFallback.push_back(childName);
            cFallbackGuideRates.push_back(preferredGuideRate);
        } else if (child.modeCategory == ProdNodeModeCategory::Group) {
            grupRateSum -= projectOnMode(child.rates, modePreferredAsWell, child.resvCoeff);
            grupGuideRateSum += preferredGuideRate;
        }
    }

    if (cFallback.empty()) return;

    const bool anyGroupWells = (grupRateSum > Scalar(0));
    const Scalar ratio = (grupGuideRateSum > Scalar(0)) ? (grupRateSum / grupGuideRateSum) : Scalar(1);

    for (size_t i = 0; i < cFallback.size(); ++i) {
        const auto& ckName = cFallback[i];
        if (tree.count(ckName) == 0) continue;
        auto& ck = tree.at(ckName);

        Scalar qk;
        if (anyGroupWells) {
            qk = cFallbackGuideRates[i] * ratio;
        } else {
            // No group-controlled wells: use limit
            const auto it = ck.Limits.find(modePreferredAsWell);
            qk = (it != ck.Limits.end()) ? it->second : Scalar(0);
        }
        // Update fallback-target
        tree.at(ckName).groupTargetFallback.value = qk;
        // qk is in the origin frame (derived from guide-rate ratios against
        // other origin-frame rates); convert to child native frame using
        // accumulated efficiency.
        const Scalar accEff = accumulatedEfficiency(tree, ckName, nodeName);
        balanceGroupTree(tree, ckName, guideRate, modePreferredAsWell, qk / accEff, tol, logger);

        // Add efficiency-scaled rates, zeroing the control-mode component to
        // avoid triggering re-balancing for this fallback child.
        std::array<Scalar, 3> scaledRates;
        for (int ph = 0; ph < 3; ++ph) {
            scaledRates[ph] = ck.efficiencyFactor * (-ck.rates[ph]);
        }
        if (modeIdx >= 0) {
            scaledRates[modeIdx] = Scalar(0);
        }
        incrementParentRateSums(tree, ckName, nodeName, std::make_optional(scaledRates));
    }
}

// ---------------------------------------------------------------------------
// Private algorithm sub-helpers (not in the header; only called from balanceGroupTree)
// ---------------------------------------------------------------------------
namespace {

/// Tighten the balancing mode and target to the most violated active limit when
/// scaling current rates to \p targetRate.  If no limit is violated the inputs
/// are returned unchanged.
///
/// \param[in] node        Tree node being balanced.
/// \param[in] mode        Initial control mode.
/// \param[in] targetRate  Initial target rate.
/// \param[in] tol         Convergence tolerance.
/// \return    Effective {mode, target} after limit tightening.
template<class Scalar>
std::pair<Well::ProducerCMode, Scalar>
tightenModeAndTarget(const ProdGroupTreeNode<Scalar>& node,
                     Well::ProducerCMode mode,
                     Scalar targetRate,
                     Scalar tol)
{
    const auto& rates = node.rates;
    const Scalar rateForModeVal = -projectOnMode(rates, mode, node.resvCoeff);
    if (rateForModeVal <= Scalar(0)) {
        return {mode, targetRate};
    }

    // Scale rates to target and find the most violated limit.
    std::array<Scalar, 3> scaledRates;
    for (int cmp = 0; cmp < 3; ++cmp) {
        scaledRates[cmp] = rates[cmp] * (targetRate / rateForModeVal);
    }

    Scalar maxViolation = Scalar(0);
    Well::ProducerCMode violatingMode = mode;
    for (const auto& [limitMode, limit] : node.Limits) {
        const Scalar violation = -projectOnMode(scaledRates, limitMode, node.resvCoeff) / limit;
        if (violation > maxViolation) {
            maxViolation = violation;
            violatingMode = limitMode;
        }
    }

    if (maxViolation > Scalar(1) - tol) {
        const Scalar newTarget = -projectOnMode(scaledRates, violatingMode, node.resvCoeff) / maxViolation;
        return {violatingMode, newTarget};
    }

    return {mode, targetRate};
}

// ---------------------------------------------------------------------------

/// Scale a well node's rates proportionally to meet \p targetRate in \p mode.
template<class Scalar>
void balanceWellNode(ProdGroupTreeNode<Scalar>& node,
                     Well::ProducerCMode mode,
                     Scalar targetRate)
{
    const Scalar currentRate = -projectOnMode(node.rates, mode, node.resvCoeff);
    if (currentRate > Scalar(0)) {
        const Scalar scale = targetRate / currentRate;
        for (int cmp = 0; cmp < 3; ++cmp) {
            node.rates[cmp] *= scale;
        }
    }
    node.isBalanced = true;
}

// ---------------------------------------------------------------------------

/// Status flags returned by runSingleDistributionPass.
template<class Scalar>
struct DistributionPassResult {
    bool needsResorting          = false;
    bool anyGroupChildren        = false;
    bool anyChildrenNeedFallback = false;
};

/// Execute one iteration of the guide-rate distribution loop for a group node.
/// Distributes rates to all sorted children, handles transparent groups, and
/// propagates residual rateSums back into transparent-group rates.
///
/// \param[in,out] tree           The production group tree.
/// \param[in]     nodeName       Origin (balancing root) node.
/// \param[in]     guideRate      Guide-rate data.
/// \param[in]     c              Children available for group control (with guide rates).
/// \param[in,out] c_trans        Transparent children; modified as nodes are balanced.
/// \param[in]     sortedIndices  Indices into \p c sorted by limit-to-guide ratio.
/// \param[in]     ratios         Limit-to-guide ratios for each element of \p c.
/// \param[in]     mode           Active control mode.
/// \param[in]     qm             Total target rate for the origin node.
/// \param[in]     tol            Convergence tolerance.
/// \param[in,out] logger         Deferred logger.
/// \return  Status flags: needsResorting, anyGroupChildren, anyChildrenNeedFallback.
template<class Scalar>
DistributionPassResult<Scalar>
runSingleDistributionPass(Tree<Scalar>& tree,
                          const std::string& nodeName,
                          const GuideRate& guideRate,
                          const std::vector<std::string>& c,
                          std::vector<std::string>& c_trans,
                          const std::vector<size_t>& sortedIndices,
                          const std::vector<Scalar>& ratios,
                          Well::ProducerCMode mode,
                          Scalar qm,
                          Scalar tol,
                          DeferredLogger& logger)
{
    DistributionPassResult<Scalar> result{};
    auto& node = tree.at(nodeName);

    for (size_t k = 0; k < sortedIndices.size(); ++k) {
        const size_t idx = sortedIndices[k];
        const auto& ckName = c[idx];

        if (tree.count(ckName) == 0) continue;
        auto& ck = tree.at(ckName);

        if (ck.useFallback) {
            // skip and flag for later fallback distribution
            result.anyChildrenNeedFallback = true;
            continue;
        }

        // Skip if path to origin was blocked by transparent node -> active limit
        if (!hasFreePath(tree, ckName, nodeName)) continue;

        // Check and update transparent groups before processing this child
        if (!c_trans.empty() && !result.anyGroupChildren) {
            auto c_trans_update = updateTransparentGroups(tree, nodeName, c_trans,
                                                          guideRate,
                                                          ratios[idx], qm, mode, tol, logger);

            if (!c_trans_update.empty()) {
                /* 
                // THIS BLOCK IS NOT NEEDED - KEEP UNTIL WE ARE SURE
                bool anyTransparent = false;
                for (const auto& ctkName : c_trans_update) {
                    if (tree.count(ctkName) > 0 && hasFreePath(tree, ctkName, nodeName)) {
                        const auto& ctk = tree.at(ctkName);
                        anyTransparent = anyTransparent || (ctk.modeCategory == ProdNodeModeCategory::Transparent);
                    }
                }
                result.anyGroupChildren = result.anyGroupChildren || anyTransparent;
                */
                if (!hasFreePath(tree, ckName, nodeName)) continue;
            }
        }

        const Scalar gsum = projectOnMode(node.guideRateSums, mode, node.resvCoeff);
        const Scalar gk = projectOnMode(ck.guideRateSums, mode, ck.resvCoeff);
        const Scalar qmRemain = qm - projectOnMode(node.rateSums, mode, node.resvCoeff);
        // qkParent is in origin frame; divide by accumulated efficiency
        // (product of all eff factors from ck up to nodeName, inclusive)
        // to convert to ck's native frame before recursing.
        const Scalar accEff = accumulatedEfficiency(tree, ckName, nodeName);
        const Scalar qkParent = gsum > Scalar(0) ? (gk / gsum) * qmRemain : Scalar(0);
        const Scalar qk = qkParent / accEff;

        const auto guideRateSumsOrig = ck.guideRateSums;

        balanceGroupTree(tree, ckName, guideRate, mode, qk, tol, logger);

        if (tree.count(ckName) > 0 && qk > Scalar(0)) {
            const auto& ckBalanced = tree.at(ckName);
            result.anyGroupChildren = result.anyGroupChildren ||
                (ckBalanced.modeCategory == ProdNodeModeCategory::Group);
            if (ckBalanced.modeCategory != ProdNodeModeCategory::Group && result.anyGroupChildren) {
                // todo: in special cases with zero rates, the categorization may be ambiguous,
                // so here we might trigger a re-sort even if it's not neccesasary. In such cases
                // we will hit the maxResortingCount limit before exiting the loop.
                result.needsResorting = true;
            }

            // Update parent sums (all nodes in c have hasGuideRate == true by construction)
            decrementParentGuideRateSums(tree, ckName, nodeName, guideRateSumsOrig);
            incrementParentRateSums(tree, ckName, nodeName);
        }
    }

    // Propagate rateSums back into rates for remaining transparent groups.
    for (const auto& ctName : c_trans) {
        if (tree.count(ctName) > 0) {
            auto& ct = tree.at(ctName);
            for (int phase = 0; phase < 3; ++phase) {
                ct.rates[phase] = -ct.rateSums[phase];
            }
            for (const auto& [limitMode, limit] : ct.Limits) {
                if (-projectOnMode(ct.rates, limitMode, ct.resvCoeff) > limit * (Scalar(1) + tol)) {
                    result.needsResorting = true;
                    break;
                }
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------

/// Result type for checkAndSwitchMode.
template<class Scalar>
struct ModeSwitchResult {
    bool limitViolated    = false;
    Well::ProducerCMode newMode;
    Scalar newTarget      = Scalar(0);
    int newTopSwitchCount = 0;
};

/// Check whether a stricter limit than the current mode is violated after one
/// distribution pass; switch mode and target if so.  Also updates the node's
/// groupTarget control mode and guide rate on a switch.
///
/// \param[in,out] node               Origin node (groupTarget fields may be updated).
/// \param[in]     nodeName           Node name (for guide-rate lookup).
/// \param[in]     currentMode        Active control mode.
/// \param[in]     targetMode         Original target mode passed to balanceGroupTree.
/// \param[in]     targetRate         Original target rate passed to balanceGroupTree.
/// \param[in]     guideRate          Guide-rate data.
/// \param[in]     topSwitchCount     Current mode-switch iteration count.
/// \param[in]     maxTopSwitchCount  Maximum allowed mode switches.
/// \param[in]     tol                Convergence tolerance.
/// \return  {limitViolated, newMode, newTarget, updatedTopSwitchCount}.
template<class Scalar>
ModeSwitchResult<Scalar>
checkAndSwitchMode(ProdGroupTreeNode<Scalar>& node,
                   const std::string& nodeName,
                   Well::ProducerCMode currentMode,
                   Well::ProducerCMode targetMode,
                   Scalar targetRate,
                   const GuideRate& guideRate,
                   int topSwitchCount,
                   int maxTopSwitchCount,
                   Scalar tol)
{
    ModeSwitchResult<Scalar> result{false, currentMode, Scalar(0), topSwitchCount};

    if (topSwitchCount >= maxTopSwitchCount) {
        return result;
    }

    const auto& rateSums = node.rateSums;
    Scalar maxViolation = Scalar(0);
    Well::ProducerCMode newMode = currentMode;

    for (const auto& [limitMode, limit] : node.Limits) {
        const Scalar effectiveLimit = (limitMode == targetMode)
            ? std::min(limit, targetRate) : limit;
        const Scalar violation = projectOnMode(rateSums, limitMode, node.resvCoeff) / effectiveLimit;
        if (violation > maxViolation) {
            maxViolation = violation;
            newMode = limitMode;
        }
    }

    if (maxViolation > Scalar(1) - tol && currentMode != newMode) {
        // if we are switching back to the original target mode, we do one final switch and then stop
        result.newTopSwitchCount = (newMode == targetMode) ? maxTopSwitchCount : topSwitchCount + 1;
        result.limitViolated = true;
        result.newMode = newMode;
        result.newTarget = projectOnMode(rateSums, newMode, node.resvCoeff) / maxViolation;
        const auto groupMode = wellModeToGroupMode(newMode);
        node.groupTarget.ctrlMode = groupMode;
        node.groupTarget.guideRate = getGuideRateForMode(nodeName, node.initialRates, groupMode, guideRate);
    }

    return result;
}

// ---------------------------------------------------------------------------

/// Assign the final modeCategory (and mode) to a node after its rates have been
/// set by balanceGroupTree.  A node is:
///   Individual  — at one of its own rate limits (within tol);
///   Group       — consuming its guide-rate share of the parent target;
///   None        — a group with no group-controlled children;
///   Transparent — no guide rate/no active individual limit.
/// If none of the above applies the node has group-controlled children yet is
/// not at a limit and not on target; a warning is logged and the node (and its
/// Group children) are reset to Individual / CMODE_UNDEFINED.
///
/// \param[in,out] tree              The production group tree.
/// \param[in]     nodeName          Node to categorize.
/// \param[in]     mode              Final balancing mode.
/// \param[in]     qm               Final balancing target rate.
/// \param[in]     anyGroupChildren  True if any child ended up Group-controlled.
/// \param[in]     tol              Convergence tolerance.
/// \param[in,out] logger            Deferred logger.
template<class Scalar>
void categorizeBalancedNode(Tree<Scalar>& tree,
                            const std::string& nodeName,
                            Well::ProducerCMode mode,
                            Scalar qm,
                            bool anyGroupChildren,
                            Scalar tol,
                            DeferredLogger& logger)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);

    node.mode = mode;
    node.visited = true;

    bool atLimit = false;
    for (const auto& [limitMode, limit] : node.Limits) {
        if (std::abs(-projectOnMode(node.rates, limitMode, node.resvCoeff) - limit) <= tol * limit) {
            atLimit = true;
            break;
        }
    }
    const bool atTarget = std::abs(-projectOnMode(node.rates, mode, node.resvCoeff) - qm) <= tol * qm;

    if (atLimit) {
        node.modeCategory = ProdNodeModeCategory::Individual;
    } else if (node.hasGuideRate && atTarget) {
        node.modeCategory = ProdNodeModeCategory::Group;
    } else if (node.type == ProdNodeType::Group && !anyGroupChildren) {
        node.modeCategory = ProdNodeModeCategory::None;
    } else if (!node.hasGuideRate) {
        // Not intended usage (calling balanceGroupTree for node without a guide rate should
        // imply active individual limit), but include for completeness
        node.modeCategory = ProdNodeModeCategory::Transparent;
    } else {
        // Problematic case: node has group controlled children, but is not at a limit and is not 
        // consuming its guide-rate share. If this occurs set modeCategory to None and reset any 
        // group-controlled children to individual to avoid subsequent issues with undefined group-control.
        logger.warning("ProdGroupTreeBalancer",
            fmt::format("categorizeBalancedNode: {} has group-controlled children but is not at a limit and is not consuming its guide-rate share.", nodeName));
        node.modeCategory = ProdNodeModeCategory::None;
        node.mode = Well::ProducerCMode::CMODE_UNDEFINED;
        for (const auto& childName : node.children) {
            if (tree.count(childName) == 0) continue;
            auto& child = tree.at(childName);
            if (child.modeCategory == ProdNodeModeCategory::Group) {
                child.modeCategory = ProdNodeModeCategory::Individual;
                if (!child.Limits.empty()) {
                    // Reset to the first limit mode if available
                    child.mode = child.Limits.begin()->first;
                } else {
                    child.mode = Well::ProducerCMode::CMODE_UNDEFINED;
                }
            }
        }
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------

template<class Scalar>
void balanceGroupTree(Tree<Scalar>& tree,
                      const std::string& nodeName,
                      const GuideRate& guideRate,
                      Well::ProducerCMode targetMode,
                      Scalar targetRate,
                      Scalar tol,
                      DeferredLogger& logger)
{
    if (tree.count(nodeName) == 0) return;
    auto& node = tree.at(nodeName);
    if (node.isSatellite) return;

    // These limits are not expected to be hit. In theory, the switch limit could
    // be reached, e.g., orat -> wrat -> grat -> orat, but very unlikely (requires
    // phase guiderate ratios sufficiently far from actual phase rate ratios.
    const int maxResortingCount = 5;
    const int maxTopSwitchCount = 3;

    auto [mode, qm] = tightenModeAndTarget(node, targetMode, targetRate, tol);

    node.balancingCount++;
    bool anyGroupChildren = false;

    if (node.type == ProdNodeType::Well) {
        balanceWellNode(node, mode, qm);
    } else {
        // Group node: recursive balancing
        bool balanced = false;
        int resortingCount = 0;
        int topSwitchCount = 0;
        int iterationCount = 0;

        // Main distribution: typically a single pass is sufficient
        while (!balanced && resortingCount <= maxResortingCount && topSwitchCount <= maxTopSwitchCount) {
            ++iterationCount;
            auto [c, c_fixed, c_trans] = getLocalTreeDescendants(tree, nodeName);

            updateGuideRatesForMode(tree, c, mode, guideRate);
            resetRatesAndGuideRateSums(tree, c, c_fixed, c_trans, nodeName, mode);

            // Add in fixed children to parent sums (e.g., satellites)
            for (const auto& cfName : c_fixed) {
                incrementParentRateSums(tree, cfName, nodeName);
            }

            const auto ratios = computeRatiosForSorting(tree, c, nodeName);

            std::vector<size_t> sortedIndices(c.size());
            std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
            std::sort(sortedIndices.begin(), sortedIndices.end(),
                     [&ratios](size_t i1, size_t i2) { return ratios[i1] < ratios[i2]; });

            auto passResult = runSingleDistributionPass(
                tree, nodeName, guideRate, c, c_trans,
                sortedIndices, ratios, mode, qm, tol, logger);

            anyGroupChildren = passResult.anyGroupChildren;

            if (passResult.needsResorting && resortingCount < maxResortingCount) {
                resortingCount++;
                continue;
            }

            resortingCount = 0;

            if (passResult.anyChildrenNeedFallback) {
                distributeFallbackRates(tree, nodeName, c, guideRate, mode, tol, logger);

                // Propagate fallback rates into transparent groups.
                for (const auto& ctName : c_trans) {
                    if (tree.count(ctName) > 0) {
                        auto& ct = tree.at(ctName);
                        for (int phase = 0; phase < 3; ++phase) {
                            ct.rates[phase] = -ct.rateSums[phase];
                        }
                    }
                }
            }

            auto switchResult = checkAndSwitchMode(
                node, nodeName, mode, targetMode, targetRate,
                guideRate, topSwitchCount, maxTopSwitchCount, tol);

            if (switchResult.limitViolated) {
                mode           = switchResult.newMode;
                qm             = switchResult.newTarget;
                topSwitchCount = switchResult.newTopSwitchCount;
            }

            balanced = !switchResult.limitViolated;
            node.isBalanced = balanced;
            for (int phase = 0; phase < 3; ++phase) {
                node.rates[phase] = -node.rateSums[phase];
            }
        }
        node.lastIterationCount = iterationCount;
    }

    categorizeBalancedNode(tree, nodeName, mode, qm, anyGroupChildren, tol, logger);
}

// ---------------------------------------------------------------------------

// Set group targets for all descendants of topName, recursively.
template<class Scalar>
void setTargets(Tree<Scalar>& tree, const std::string& topName)
{
    if (tree.count(topName) == 0) return;
    const auto& top = tree.at(topName);
    if (top.children.empty()) return;

    // Collect non-transparent descendants by expanding through transparent groups
    // (iterative DFS, same traversal order as a recursive descent).
    std::vector<std::string> cList;
    {
        std::vector<std::string> toExpand(top.children.rbegin(), top.children.rend());
        while (!toExpand.empty()) {
            const auto name = std::move(toExpand.back());
            toExpand.pop_back();
            if (tree.count(name) == 0) continue;
            const auto& node = tree.at(name);
            if (node.modeCategory == ProdNodeModeCategory::Transparent) {
                toExpand.insert(toExpand.end(), node.children.rbegin(), node.children.rend());
            } else {
                cList.push_back(name);
            }
        }
    }

    if (top.modeCategory == ProdNodeModeCategory::None) {
        for (const auto& cName : cList) {
            if (tree.count(cName) == 0) continue;
            auto& c = tree.at(cName);
            c.groupTarget.groupName = topName;
            c.groupTarget.ctrlMode  = Group::ProductionCMode::NONE;
            c.groupTarget.value     = std::numeric_limits<Scalar>::max();
            setTargets(tree, cName);
        }
        return;
    }

    const Well::ProducerCMode    mode          = top.mode;
    const Group::ProductionCMode modePreferred = top.preferredMode;
    const bool modeIsPref = (wellModeToGroupMode(mode) == modePreferred);

    const bool anyFallback = !modeIsPref &&
        std::any_of(cList.begin(), cList.end(), [&](const std::string& n) {
            return tree.count(n) > 0 && tree.at(n).useFallback;
        });

    // Group name to report in each child's target
    const std::string& groupName = (top.modeCategory == ProdNodeModeCategory::Group)
        ? top.groupTarget.groupName : topName;

    // -----------------------------------------------------------------------
    // First pass: compute aggregates needed for guideRateRatio and the
    // Individual-node hypothetical GRUP target.
    //
    // Active mode:
    //   guideSum     - total active-mode guide rates (0 for useFallback children)
    //   guideSumGrup - active-mode guide rates of GRUP non-fallback children
    //   targetSum    - top rate minus each non-GRUP child's parent-frame
    //                  contribution = budget available to GRUP children.
    //                  Used to compute the hypothetical GRUP target for
    //                  Individual children (limit used to check mode switch).
    //
    // Preferred mode (only when !modeIsPref):
    //   guideSumFallback - total preferred-mode guide rates (all children);
    //                      needed only for guideRateRatio.
    // -----------------------------------------------------------------------
    Scalar guideSum     = Scalar(0);
    Scalar guideSumGrup = Scalar(0);
    Scalar targetSum    = -projectOnMode(top.rates, mode, top.resvCoeff);

    Scalar guideSumFallback = Scalar(0);
    const Well::ProducerCMode modePreferredAsWell = groupModeToWellMode(modePreferred);

    for (const auto& cName : cList) {
        if (tree.count(cName) == 0) continue;
        const auto& c = tree.at(cName);
        const bool isGrup     = (c.modeCategory == ProdNodeModeCategory::Group);
        const bool isFallback = anyFallback && c.useFallback;
        const Scalar gr = isFallback ? Scalar(0) : c.groupTarget.guideRate;

        guideSum += gr;
        if (isGrup) {
            guideSumGrup += gr;
        } else if (!isFallback) {
            // Individual child: remove its parent-frame contribution so GRUP
            // children share only the remaining budget.
            targetSum += c.efficiencyFactor * projectOnMode(c.rates, mode, c.resvCoeff);
        }

        if (!modeIsPref) {
            guideSumFallback += c.groupTargetFallback.guideRate;
        }
    }

    // -----------------------------------------------------------------------
    // Second pass: assign group targets and recurse.
    // -----------------------------------------------------------------------
    for (const auto& cName : cList) {
        if (tree.count(cName) == 0) continue;
        auto& c = tree.at(cName);
        if (c.isSatellite) continue; // skip satellites
        const bool isGrup     = (c.modeCategory == ProdNodeModeCategory::Group);
        const bool isFallback = anyFallback && c.useFallback;
        const Scalar gr = isFallback ? Scalar(0) : c.groupTarget.guideRate;

        // Active-mode target
        c.groupTarget.groupName = groupName;
        c.groupTarget.ctrlMode  = wellModeToGroupMode(mode);
        c.groupTarget.guideRate = gr;

        if (!modeIsPref) {
            c.groupTarget.guideRateRatio = (guideSum > Scalar(0)) ? gr / guideSum : Scalar(0);
        }

        if (isFallback) {
            // Controlled in preferred mode via distributeFallbackRates; there is no
            // meaningful active-mode target for this node.
            c.groupTarget.value = std::numeric_limits<Scalar>::quiet_NaN();
        } else if (isGrup) {
            // GRUP child: actual balanced rate IS the allocated target.
            c.groupTarget.value = -projectOnMode(c.rates, mode, c.resvCoeff);
        } else {
            // Individual child: compute the hypothetical GRUP target — i.e. what
            // this child would receive if it were group-controlled.  This value is
            // used downstream to check whether the well should switch to GRUP mode.
            const Scalar totalGuide  = guideSumGrup + gr;
            const Scalar totalTarget = targetSum
                - c.efficiencyFactor * projectOnMode(c.rates, mode, c.resvCoeff);
            const Scalar alloc = (totalGuide > Scalar(0))
                ? (gr / totalGuide) * totalTarget : Scalar(0);
            c.groupTarget.value = (c.efficiencyFactor > Scalar(0))
                ? alloc / c.efficiencyFactor : Scalar(0);
        }

        // Preferred-mode fallback target (only when active mode differs from preferred).
        // Use the actual balanced rate in preferred mode — the fractions have already
        // been applied by distributeFallbackRates / balanceGroupTree.
        if (!modeIsPref) {
            const Scalar grFb = c.groupTargetFallback.guideRate;
            c.groupTargetFallback.groupName      = groupName;
            c.groupTargetFallback.ctrlMode       = modePreferred;
            c.groupTargetFallback.value          = -projectOnMode(c.rates, modePreferredAsWell, c.resvCoeff);
            c.groupTargetFallback.guideRateRatio = (guideSumFallback > Scalar(0))
                ? grFb / guideSumFallback : Scalar(0);
        }

        setTargets(tree, cName);
    }
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
bool runBalancingAlgorithm(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                           Tree<Scalar>& tree, Scalar tol, DeferredLogger& logger)
{
    // Main algorithm entry point
    // Get sub-tree ordering - each "top" node will be balanced independently, 
    // starting from the bottom of the tree
    const auto topnodes = getSubTreeOrdering(tree, "FIELD");

    if (topnodes.empty()) {
        // No nodes to balance
        return true;
    }

    for (const auto& nodeName : topnodes) {
        if (tree.count(nodeName) == 0) continue;
        auto& node = tree.at(nodeName);
        if (node.isSatellite) continue; // Skip satellite groups

        // Determine mode and target
        // A top node is either individual or NONE (in which case it has a preferredMode)
        Well::ProducerCMode mode;
        if (node.modeCategory == ProdNodeModeCategory::Individual || node.type == ProdNodeType::Well) {
            mode = node.mode;
        } else {
            // Convert preferredMode (Group::ProductionCMode) to Well::ProducerCMode
            mode = groupModeToWellMode(node.preferredMode);
        }
        if (mode == Well::ProducerCMode::CMODE_UNDEFINED || mode == Well::ProducerCMode::NONE) {
            logger.warning("ProdGroupTreeBalancer",
                fmt::format("runBalancingAlgorithm: top node '{}' has undefined mode, returning", nodeName));
            return false;
        }

        // Get limit for mode
        Scalar qm = 0.0;
        if (node.Limits.count(mode) > 0) {
            qm = node.Limits.at(mode);
        } else if (wellModel.comm().rank() == 0) {
            logger.warning("ProdGroupTreeBalancer",
                fmt::format("runBalancingAlgorithm: top node '{}' has no limit for given mode, returning", nodeName));
            return false;
        }

        // Balance this subtree
        balanceGroupTree(tree, nodeName, wellModel.guideRate(), mode, qm, tol, logger);

        // setTargets is intentionally not called here; group target values
        // (groupName, value, guideRateRatio) are populated by getWellGroupTargetProducer.
        //setTargets(tree, nodeName);

        // Report completion of this top node's balancing
        if (tree.at(nodeName).type == ProdNodeType::Group && wellModel.comm().rank() == 0) {
            const auto& n = tree.at(nodeName);
            logger.debug("ProdGroupTreeBalancer",
                fmt::format("Balancer: Completed for top node '{}': balanced={}, iterations={}",
                            nodeName,
                            n.isBalanced ? "yes" : "no",
                            n.lastIterationCount));
        }
    }

    // Update reporting-only groups that sit above the balanced subtrees.
    // These are group nodes that have no limits of their own and therefore
    // never entered the main balancing loop.
    //
    // Traversal: post-order DFS from FIELD, stopping at (already balanced) topnodes.
    {
        auto isTopNode = [&topnodes](const std::string& name) {
            return std::ranges::find(topnodes, name) != topnodes.end();
        };
        std::function<void(const std::string&)> updateReportingGroups =
            [&](const std::string& name)
        {
            if (tree.count(name) == 0) return;
            if (isTopNode(name)) return; // already balanced — do not touch
            auto& node = tree.at(name);
            if (node.isSatellite) return;
            if (node.type == ProdNodeType::Well) return;
            // Recurse into children first so their rates are up-to-date.
            for (const auto& childName : node.children) {
                updateReportingGroups(childName);
            }
            // Sum up efficiency-weighted child contributions (same convention
            // as incrementParentRateSums): child produces (-child.rates) in its
            // own frame; the parent sees child.efficiencyFactor of that.
            node.rates = {Scalar(0), Scalar(0), Scalar(0)};
            for (const auto& childName : node.children) {
                if (tree.count(childName) == 0) continue;
                const auto& child = tree.at(childName);
                for (int c = 0; c < 3; ++c) {
                    node.rates[c] -= child.efficiencyFactor * (-child.rates[c]);
                }
            }
            node.modeCategory = ProdNodeModeCategory::None;
            node.mode         = Well::ProducerCMode::CMODE_UNDEFINED;
        };
        updateReportingGroups("FIELD");
    }

    return true;
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

template<class Scalar>
bool checkTreeValidity(const Tree<Scalar>& tree,
                       const std::string& topName,
                       Scalar tol,
                       DeferredLogger& logger)
{
    bool valid = true;

    std::function<void(const std::string&)> check = [&](const std::string& name) {
        if (tree.count(name) == 0) return;
        const auto& node = tree.at(name);

        const auto category = node.modeCategory;

        // Check that Individual nodes are at their limit
        if (category == ProdNodeModeCategory::Individual) {
            const auto it = node.Limits.find(node.mode);
            if (it != node.Limits.end()) {
                const Scalar limit = it->second;
                const Scalar current = (node.mode == Well::ProducerCMode::BHP)
                    ? std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                       [](Scalar s, Scalar r){ return s + (-r); })
                    : -projectOnMode(node.rates, node.mode, node.resvCoeff);
                const Scalar relErr = std::abs(current - limit) / (std::abs(limit) + kFeasibilityTolerance<Scalar>);
                if (relErr > tol) {
                    logger.warning("ProdGroupTreeBalancer",
                        fmt::format("Node '{}' is Individual but current rate ({:.4g}) "
                                    "differs from limit ({:.4g}) by {:.2g}%%",
                                    name, current, limit, 100.0 * relErr));
                    valid = false;
                }
            }
        }

        // Check that Group nodes satisfy their group target
        // skipped - only relevant if setTargets has been called
        /*
        if (node.modeCategory == ProdNodeModeCategory::Group &&
            node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
            const Scalar target = node.groupTarget.value;
            const Scalar current = std::accumulate(node.rates.begin(), node.rates.end(), Scalar(0),
                                                    [](Scalar s, Scalar r){ return s + (-r); });
            const Scalar relErr = std::abs(current - target) / (std::abs(target) + kFeasibilityTolerance<Scalar>);
            if (false) {//(relErr > tol) {
                logger.warning("ProdGroupTreeBalancer",
                    fmt::format("Node '{}' is Group but current total rate ({:.4g}) "
                                "differs from group target ({:.4g}) by {:.2g}%%",
                                name, current, target, 100.0 * relErr));
                valid = false;
            }
        }
        */
        // Recurse
        for (const auto& child : node.children) { check(child); }
    };

    check(topName);
    return valid;
}

// ---------------------------------------------------------------------------
// Debug output
// ---------------------------------------------------------------------------

// This function is intended for debugging. Call it after buildTree or after balancing
// to dump the tree state to the DBG log file.
template<class Scalar>
void logTree(const Tree<Scalar>& tree, DeferredLogger& logger)
{
    constexpr Scalar s_to_day = Scalar(86400);

    auto catLetter = [](ProdNodeModeCategory cat) -> char {
        switch (cat) {
        case ProdNodeModeCategory::Individual:  return 'I';
        case ProdNodeModeCategory::Group:       return 'G';
        case ProdNodeModeCategory::Transparent: return 'T';
        case ProdNodeModeCategory::None:        return 'N';
        }
        return '?';
    };

    // Build the full control chain from startName upward, collecting every
    // Group-controlled ancestor and stopping (inclusively) at the first
    // Individual ancestor.  This produces strings like "G2 < G1" so the
    // caller can render "W < G2 < G1".  Non-group/non-individual ancestors
    // (Transparent, None) are skipped without being added to the chain.
    auto buildCtrlChain = [&](const std::string& startName) -> std::string {
        std::vector<std::string> chain;
        std::string cur = startName;
        while (tree.count(cur) > 0) {
            const std::string& par = tree.at(cur).parent;
            if (par.empty() || tree.count(par) == 0) break;
            const auto cat = tree.at(par).modeCategory;
            if (cat == ProdNodeModeCategory::Group ||
                cat == ProdNodeModeCategory::Individual) {
                chain.push_back(par);
                if (cat == ProdNodeModeCategory::Individual) break;
            }
            cur = par;
        }
        if (chain.empty()) return "";
        std::string result;
        for (const auto& n : chain) {
            if (!result.empty()) result += " < ";
            result += n;
        }
        return result;
    };

    // BFS from FIELD so lines follow natural tree order.
    std::vector<std::string> queue;
    if (tree.count("FIELD") > 0) {
        queue.push_back("FIELD");
    }
    while (!queue.empty()) {
        const auto name = std::move(queue.front());
        queue.erase(queue.begin());

        if (tree.count(name) == 0) continue;
        const auto& node = tree.at(name);

        // Rates in m3/day (positive = production).
        const Scalar qo = s_to_day * (-node.rates[kOil]);
        const Scalar qw = s_to_day * (-node.rates[kWater]);
        const Scalar qg = s_to_day * (-node.rates[kGas]);
        if (qo + qw + qg <= Scalar(0) && node.type == ProdNodeType::Group) continue;
        // Line 1: name, type tag, mode category, control info, and phase rates.
        std::string ctrlInfo;
        switch (node.modeCategory) {
        case ProdNodeModeCategory::Individual:
            if (node.mode != Well::ProducerCMode::CMODE_UNDEFINED &&
                node.mode != Well::ProducerCMode::NONE) {
                ctrlInfo = fmt::format("[I] {} {:.4g}", WellProducerCMode2String(node.mode),
                    s_to_day * (-projectOnMode(node.rates, node.mode, node.resvCoeff)));
            } else {
                ctrlInfo = "[I] undefined";
            }
            break;
        case ProdNodeModeCategory::Group: {
            // groupTarget.groupName is only populated when setTargets has been called.
            // If empty, build the full control chain up to the nearest Individual ancestor.
            const std::string ctrlChain = node.groupTarget.groupName.empty()
                ? buildCtrlChain(name)
                : node.groupTarget.groupName;
            ctrlInfo = fmt::format("[G] GRUP < {}", ctrlChain.empty() ? "?" : ctrlChain);
            break;
        }
        case ProdNodeModeCategory::Transparent:
            ctrlInfo = "[T]";
            break;
        case ProdNodeModeCategory::None:
            ctrlInfo = "[N]";
            break;
        }

        const std::string typeTag = (node.type == ProdNodeType::Well) ? "well " : "group";
        logger.debug("ProdGroupTreeBalancer",
            fmt::format("  {:<25} ({})  {:<30}  O={:.1f} W={:.1f} G={:.1f} m3/d",
                        name, typeTag, ctrlInfo, qo, qw, qg));

        // Line 2 (groups only): direct children with mode-category letter.
        if (node.type == ProdNodeType::Group && !node.children.empty()) {
            std::string childrenStr;
            for (const auto& childName : node.children) {
                if (tree.count(childName) == 0) continue;
                if (!childrenStr.empty()) childrenStr += "  ";
                const char letter = (tree.count(childName) > 0)
                    ? catLetter(tree.at(childName).modeCategory)
                    : '?';
                childrenStr += fmt::format("{}({})", childName, letter);
            }
            logger.debug("ProdGroupTreeBalancer",
                fmt::format("    children: {}", childrenStr));
        }

        for (const auto& childName : node.children) {
            queue.push_back(childName);
        }
    }
}

// ---------------------------------------------------------------------------

// Write the results of the tree balancing back into the well and group states.
template<class Scalar, typename IndexTraits>
void applyTreeToState(const Tree<Scalar>& tree,
                      BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                      DeferredLogger& logger)
{
    auto& wellState  = wellModel.wellState();
    auto& groupState = wellModel.groupState();

    for (const auto& [name, node] : tree) {
        if (node.type == ProdNodeType::Well) {
            if (!wellState.has(name)) continue;
            auto& ws = wellState.well(name);
            if (!ws.producer) continue;

            // Convert canonical 3-component rates back to active-phase vector
            const auto activeRates = toActive(node.rates, ws.pu);
            ws.surface_rates = activeRates;

            // Update control mode
            const Well::ProducerCMode oldWellCMode = ws.production_cmode;
            if (node.modeCategory== ProdNodeModeCategory::Group) {
                ws.production_cmode = Well::ProducerCMode::GRUP;
            } else if (node.modeCategory== ProdNodeModeCategory::Individual) {
                ws.production_cmode = node.mode;
            }
            if (ws.production_cmode != oldWellCMode) {
                logger.debug("ProdGroupTreeBalancer",
                    fmt::format("Balancer: Well {}: production_cmode changed from {} to {}",
                                name,
                                WellProducerCMode2String(oldWellCMode),
                                WellProducerCMode2String(ws.production_cmode)));
            }

            // Update the group target stored in the well state
            if (node.groupTarget.ctrlMode != Group::ProductionCMode::NONE) {
                typename SingleWellState<Scalar, IndexTraits>::GroupTarget gt;
                gt.group_name       = node.groupTarget.groupName;
                gt.production_cmode = node.groupTarget.ctrlMode;
                gt.target_value     = node.groupTarget.value;
                gt.guiderate_ratio  = node.groupTarget.guideRateRatio;
                ws.group_target     = gt;
            } else {
                ws.group_target = std::nullopt;
            }

            // Update the fallback group target when the active control mode
            // differs from the preferred mode (set by setTargets when !modeIsPref).
            if (node.groupTargetFallback.ctrlMode != Group::ProductionCMode::NONE &&
                node.groupTargetFallback.ctrlMode != node.groupTarget.ctrlMode) {
                typename SingleWellState<Scalar, IndexTraits>::GroupTarget gtFb;
                gtFb.group_name       = node.groupTargetFallback.groupName;
                gtFb.production_cmode = node.groupTargetFallback.ctrlMode;
                gtFb.target_value     = node.groupTargetFallback.value;
                gtFb.guiderate_ratio  = node.groupTargetFallback.guideRateRatio;
                ws.group_target_fallback = gtFb;
            } else {
                ws.group_target_fallback = std::nullopt;
            }
        } else {
            if (node.isSatellite) continue; // Skip satellite groups
            // Group node: update group state
            // Rates in groupState are stored in active-phase order, positive = production
            const auto& pu = wellModel.phaseUsage();
            std::vector<Scalar> activeRates(pu.numPhases, Scalar(0));
            for (int c = 0; c < 3; ++c) {
                const int a = activeIdx(pu, c);
                if (a >= 0) {
                    // GroupState stores positive = production; our rates are negative = production
                    activeRates[a] = -node.rates[c];
                }
            }
            groupState.update_production_rates(name, activeRates);

            // Update group control mode.
            // applyTreeToState() now runs on all MPI ranks simultaneously (each rank
            // runs the same deterministic algorithm on globally consistent inputs and
            // reaches the same result), so production_control writes are valid.
            const Group::ProductionCMode oldGroupCMode =
                groupState.has_production_control(name)
                    ? groupState.production_control(name)
                    : Group::ProductionCMode::NONE;
            Group::ProductionCMode newGroupCMode = Group::ProductionCMode::NONE;
            if (node.modeCategory == ProdNodeModeCategory::Group && node.hasLimitedAncestor) {
                newGroupCMode = Group::ProductionCMode::FLD;
            } else if (node.modeCategory == ProdNodeModeCategory::Individual) {
                    switch (node.mode) {
                    case Well::ProducerCMode::ORAT: newGroupCMode = Group::ProductionCMode::ORAT; break;
                    case Well::ProducerCMode::WRAT: newGroupCMode = Group::ProductionCMode::WRAT; break;
                    case Well::ProducerCMode::GRAT: newGroupCMode = Group::ProductionCMode::GRAT; break;
                    case Well::ProducerCMode::LRAT: newGroupCMode = Group::ProductionCMode::LRAT; break;
                    case Well::ProducerCMode::RESV: newGroupCMode = Group::ProductionCMode::RESV; break;
                    default: break;
                }
            }
            groupState.production_control(name, newGroupCMode);
            if (newGroupCMode != oldGroupCMode && wellModel.comm().rank() == 0) {
                logger.debug("ProdGroupTreeBalancer",
                    fmt::format("Balancer: Group '{}': production_control changed from {} to {}",
                                name,
                                Group::ProductionCMode2String(oldGroupCMode),
                                Group::ProductionCMode2String(newGroupCMode)));
            }
        }
    }
}

// ---------------------------------------------------------------------------

template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                          DeferredLogger& logger)
{
    OPM_TIMEFUNCTION();

    const auto t0 = std::chrono::steady_clock::now();

    auto tree = buildTree(wellModel, summaryState, reportStep, limits);

    const bool success = runBalancingAlgorithm(wellModel, tree, tol, logger);

    if (wellModel.comm().rank() == 0) {
        logTree(tree, logger);
    }

    if (!success && wellModel.comm().rank() == 0) {
        logger.warning("ProdGroupTreeBalancer",
            "Group tree balancer did not converge for all subtrees. "
            "Results may be inconsistent.");
    }

    const bool valid = checkTreeValidity(tree, "FIELD", tol, logger);

    if (!valid && success && wellModel.comm().rank() == 0) {
        logger.warning("ProdGroupTreeBalancer",
            "Group tree balancer converged but validity check failed. "
            "Some nodes may not satisfy their constraints.");
    }

    applyTreeToState(tree, wellModel, logger);

    if (wellModel.comm().rank() == 0) {
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - t0).count();
        logger.debug("ProdGroupTreeBalancer",
            fmt::format("Group tree balancer completed in {}ms. "
                        "Convergence: {}, Validity: {}",
                        elapsed, success ? "OK" : "FAILED", valid ? "OK" : "FAILED"));
    }

    return valid;
}

// ===========================================================================
// Explicit instantiations
// ===========================================================================

// Only the top-level entry point requires an explicit instantiation; all
// internal helper functions are implicitly instantiated as part of it.

template bool runGroupTreeBalancer<double, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, double,
    const std::unordered_map<std::string, std::pair<int, double>>&,
    DeferredLogger&);

#ifdef FLOW_INSTANTIATE_FLOAT

template bool runGroupTreeBalancer<float, BlackOilDefaultFluidSystemIndices>(
    BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>&,
    const SummaryState&, int, float,
    const std::unordered_map<std::string, std::pair<int, float>>&,
    DeferredLogger&);

#endif

} // namespace Opm::ProdGroupTreeBalancer
