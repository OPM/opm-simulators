/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellInterfaceEval.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Opm
{

template<class FluidSystem>
WellInterfaceEval<FluidSystem>::
WellInterfaceEval(const WellInterfaceFluidSystem<FluidSystem>& baseif)
  : baseif_(baseif)
{}

template<class FluidSystem>
template<class EvalWell>
void
WellInterfaceEval<FluidSystem>::
getGroupInjectionControl(const Group& group,
                         const WellState& well_state,
                         const GroupState& group_state,
                         const Schedule& schedule,
                         const SummaryState& summaryState,
                         const InjectorType& injectorType,
                         const EvalWell& bhp,
                         const EvalWell& injection_rate,
                         EvalWell& control_eq,
                         double efficiencyFactor,
                         DeferredLogger& deferred_logger) const
{
    // Setting some defaults to silence warnings below.
    // Will be overwritten in the switch statement.
    Phase injectionPhase = Phase::WATER;
    switch (injectorType) {
    case InjectorType::WATER:
    {
        injectionPhase = Phase::WATER;
        break;
    }
    case InjectorType::OIL:
    {
        injectionPhase = Phase::OIL;
        break;
    }
    case InjectorType::GAS:
    {
        injectionPhase = Phase::GAS;
        break;
    }
    default:
        // Should not be here.
        assert(false);
    }

    auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
    if (currentGroupControl == Group::InjectionCMode::FLD ||
        currentGroupControl == Group::InjectionCMode::NONE) {
        if (!group.injectionGroupControlAvailable(injectionPhase)) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            const auto& controls = baseif_.wellEcl().injectionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        } else {
            // Inject share of parents control
            const auto& parent = schedule.getGroup( group.parent(), baseif_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            getGroupInjectionControl(parent, well_state, group_state, schedule, summaryState, injectorType, bhp, injection_rate, control_eq, efficiencyFactor, deferred_logger);
            return;
        }
    }

    const auto& well = baseif_.wellEcl();
    const auto pu = baseif_.phaseUsage();

    if (!group.isInjectionGroup()) {
        // use bhp as control eq and let the updateControl code find a valid control
        const auto& controls = well.injectionControls(summaryState);
        control_eq = bhp - controls.bhp_limit;
        return;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(pu.num_phases, 1.0);
    baseif_.rateConverter().calcInjCoeff(0, baseif_.pvtRegionIdx(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    double sales_target = 0;
    if (schedule[baseif_.currentStep()].gconsale().has(group.name())) {
        const auto& gconsale = schedule[baseif_.currentStep()].gconsale().get(group.name(), summaryState);
        sales_target = gconsale.sales_target;
    }
    WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu, resv_coeff, group.name(), sales_target, group_state, injectionPhase, group.has_gpmaint_control(injectionPhase, currentGroupControl), deferred_logger);
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, baseif_.currentStep(), baseif_.guideRate(), tcalc.guideTargetMode(), pu, false, injectionPhase);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child);
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    const double orig_target = tcalc.groupTarget(group.injectionControls(injectionPhase, summaryState), deferred_logger);
    const auto chain = WellGroupHelpers::groupChainTopBot(baseif_.name(), group.name(), schedule, baseif_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || baseif_.guideRate()->has(chain[ii], injectionPhase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto current_rate = injection_rate;

    control_eq = current_rate - target_rate;
}

template<class FluidSystem>
template<class EvalWell>
void
WellInterfaceEval<FluidSystem>::
getGroupProductionControl(const Group& group,
                          const WellState& well_state,
                          const GroupState& group_state,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          const EvalWell& bhp,
                          const std::vector<EvalWell>& rates,
                          EvalWell& control_eq,
                          double efficiencyFactor) const
{
    const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
    if (currentGroupControl == Group::ProductionCMode::FLD ||
        currentGroupControl == Group::ProductionCMode::NONE) {
        if (!group.productionGroupControlAvailable()) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            const auto& controls = baseif_.wellEcl().productionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        } else {
            // Produce share of parents control
            const auto& parent = schedule.getGroup(group.parent(), baseif_.currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            getGroupProductionControl(parent, well_state, group_state, schedule, summaryState, bhp, rates, control_eq, efficiencyFactor);
            return;
        }
    }

    const auto& well = baseif_.wellEcl();
    const auto pu = baseif_.phaseUsage();

    if (!group.isProductionGroup()) {
        // use bhp as control eq and let the updateControl code find a valid control
        const auto& controls = well.productionControls(summaryState);
        control_eq = bhp - controls.bhp_limit;
        return;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(baseif_.phaseUsage().num_phases, 1.0);
    baseif_.rateConverter().calcCoeff(0, baseif_.pvtRegionIdx(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    // gconsale may adjust the grat target.
    // the adjusted rates is send to the targetCalculator
    double gratTargetFromSales = 0.0;
    if (group_state.has_grat_sales_target(group.name()))
        gratTargetFromSales = group_state.grat_sales_target(group.name());

    WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff, gratTargetFromSales, group.name(), group_state, group.has_gpmaint_control(currentGroupControl));
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, baseif_.currentStep(), baseif_.guideRate(), tcalc.guideTargetMode(), pu, true, Phase::OIL);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child);
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    const double orig_target = tcalc.groupTarget(group.productionControls(summaryState));
    const auto chain = WellGroupHelpers::groupChainTopBot(baseif_.name(), group.name(), schedule, baseif_.currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || baseif_.guideRate()->has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    control_eq = current_rate - target_rate;
}

template<class FluidSystem>
template<class EvalWell>
void
WellInterfaceEval<FluidSystem>::
assembleControlEqProd_(const WellState& well_state,
                       const GroupState& group_state,
                       const Schedule& schedule,
                       const SummaryState& summaryState,
                       const Well::ProductionControls& controls,
                       const EvalWell& bhp,
                       const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                       const std::function<EvalWell()>& bhp_from_thp,
                       EvalWell& control_eq,
                       DeferredLogger& deferred_logger) const
{
    const auto current = well_state.well(baseif_.indexOfWell()).production_cmode;
    const auto& pu = baseif_.phaseUsage();
    const double efficiencyFactor = baseif_.wellEcl().getEfficiencyFactor();

    switch (current) {
    case Well::ProducerCMode::ORAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Liquid];
        control_eq = rate - controls.oil_rate;
        break;
    }
    case Well::ProducerCMode::WRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Aqua];
        control_eq = rate - controls.water_rate;
        break;
    }
    case Well::ProducerCMode::GRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Vapour];
        control_eq = rate - controls.gas_rate;
        break;
    }
    case Well::ProducerCMode::LRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        EvalWell rate = -rates[BlackoilPhases::Aqua] - rates[BlackoilPhases::Liquid];
        control_eq = rate - controls.liquid_rate;
        break;
    }
    case Well::ProducerCMode::CRAT: {
        OPM_DEFLOG_THROW(std::runtime_error, "CRAT control not supported " << baseif_.name(), deferred_logger);
    }
    case Well::ProducerCMode::RESV: {
        auto total_rate = rates[0]; // To get the correct type only.
        total_rate = 0.0;
        std::vector<double> convert_coeff(baseif_.numPhases(), 1.0);
        baseif_.rateConverter().calcCoeff(/*fipreg*/ 0, baseif_.pvtRegionIdx(), convert_coeff);
        for (int phase = 0; phase < 3; ++phase) {
            if (pu.phase_used[phase]) {
                const int pos = pu.phase_pos[phase];
                total_rate -= rates[phase] * convert_coeff[pos]; // Note different indices.
            }
        }
        if (controls.prediction_mode) {
            control_eq = total_rate - controls.resv_rate;
        } else {
            std::vector<double> hrates(baseif_.numPhases(), 0.);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                hrates[pu.phase_pos[Water]] = controls.water_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                hrates[pu.phase_pos[Oil]] = controls.oil_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                hrates[pu.phase_pos[Gas]] = controls.gas_rate;
            }
            std::vector<double> hrates_resv(baseif_.numPhases(), 0.);
            baseif_.rateConverter().calcReservoirVoidageRates(/*fipreg*/ 0, baseif_.pvtRegionIdx(), hrates, hrates_resv);
            double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
            control_eq = total_rate - target;
        }
        break;
    }
    case Well::ProducerCMode::BHP: {
        control_eq = bhp - controls.bhp_limit;
        break;
    }
    case Well::ProducerCMode::THP: {
        control_eq = bhp - bhp_from_thp();
        break;
    }
    case Well::ProducerCMode::GRUP: {
        assert(baseif_.wellEcl().isAvailableForGroupControl());
        const auto& group = schedule.getGroup(baseif_.wellEcl().groupName(), baseif_.currentStep());
        // Annoying thing: the rates passed to this function are
        // always of size 3 and in canonical (for PhaseUsage)
        // order. This is what is needed for VFP calculations if
        // they are required (THP controlled well). But for the
        // group production control things we must pass only the
        // active phases' rates.
        std::vector<EvalWell> active_rates(pu.num_phases);
        for (int canonical_phase = 0; canonical_phase < 3; ++canonical_phase) {
            if (pu.phase_used[canonical_phase]) {
                active_rates[pu.phase_pos[canonical_phase]] = rates[canonical_phase];
            }
        }
        getGroupProductionControl(group, well_state, group_state, schedule, summaryState, bhp, active_rates, control_eq, efficiencyFactor);
        break;
    }
    case Well::ProducerCMode::CMODE_UNDEFINED: {
        OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + baseif_.name(), deferred_logger);
    }
    case Well::ProducerCMode::NONE: {
        OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + baseif_.name(), deferred_logger);
    }
    }
}

template<class FluidSystem>
template<class EvalWell>
void
WellInterfaceEval<FluidSystem>::
assembleControlEqInj_(const WellState& well_state,
                      const GroupState& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      const Well::InjectionControls& controls,
                      const EvalWell& bhp,
                      const EvalWell& injection_rate,
                      const std::function<EvalWell()>& bhp_from_thp,
                      EvalWell& control_eq,
                      DeferredLogger& deferred_logger) const
{
    auto current = well_state.well(baseif_.indexOfWell()).injection_cmode;
    const InjectorType injectorType = controls.injector_type;
    const auto& pu = baseif_.phaseUsage();
    const double efficiencyFactor = baseif_.wellEcl().getEfficiencyFactor();

    switch (current) {
    case Well::InjectorCMode::RATE: {
        control_eq = injection_rate - controls.surface_rate;
        break;
    }
    case Well::InjectorCMode::RESV: {
        std::vector<double> convert_coeff(baseif_.numPhases(), 1.0);
        baseif_.rateConverter().calcInjCoeff(/*fipreg*/ 0, baseif_.pvtRegionIdx(), convert_coeff);

        double coeff = 1.0;

        switch (injectorType) {
        case InjectorType::WATER: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Aqua]];
            break;
        }
        case InjectorType::OIL: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Liquid]];
            break;
        }
        case InjectorType::GAS: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Vapour]];
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + baseif_.wellEcl().name());
        }

        control_eq = coeff * injection_rate - controls.reservoir_rate;
        break;
    }
    case Well::InjectorCMode::THP: {
        control_eq = bhp - bhp_from_thp();
        break;
    }
    case Well::InjectorCMode::BHP: {
        control_eq = bhp - controls.bhp_limit;
        break;
    }
    case Well::InjectorCMode::GRUP: {
        assert(baseif_.wellEcl().isAvailableForGroupControl());
        const auto& group = schedule.getGroup(baseif_.wellEcl().groupName(), baseif_.currentStep());
        this->getGroupInjectionControl(group,
                                       well_state,
                                       group_state,
                                       schedule,
                                       summaryState,
                                       injectorType,
                                       bhp,
                                       injection_rate,
                                       control_eq,
                                       efficiencyFactor,
                                       deferred_logger);
        break;
    }
    case Well::InjectorCMode::CMODE_UNDEFINED: {
        OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + baseif_.name(), deferred_logger);
    }
    }
}

template<class FluidSystem>
template<class EvalWell>
EvalWell
WellInterfaceEval<FluidSystem>::
calculateBhpFromThp(const WellState& well_state,
                    const std::vector<EvalWell>& rates,
                    const Well& well,
                    const SummaryState& summaryState,
                    const double rho,
                    DeferredLogger& deferred_logger) const
{
    // TODO: when well is under THP control, the BHP is dependent on the rates,
    // the well rates is also dependent on the BHP, so it might need to do some iteration.
    // However, when group control is involved, change of the rates might impacts other wells
    // so iterations on a higher level will be required. Some investigation might be needed when
    // we face problems under THP control.

    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    const EvalWell aqua = rates[Water];
    const EvalWell liquid = rates[Oil];
    const EvalWell vapour = rates[Gas];

    // pick the reference density
    // typically the reference in the top layer
    if (baseif_.isInjector() )
    {
        const auto& controls = well.injectionControls(summaryState);
        const double vfp_ref_depth = baseif_.vfpProperties()->getInj()->getTable(controls.vfp_table_number).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, rho, baseif_.gravity());
        return baseif_.vfpProperties()->getInj()->bhp(controls.vfp_table_number, aqua, liquid, vapour, baseif_.getTHPConstraint(summaryState)) - dp;
     }
     else if (baseif_.isProducer()) {
         const auto& controls = well.productionControls(summaryState);
         const double vfp_ref_depth = baseif_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
         const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, rho, baseif_.gravity());
         const auto& wfr =  baseif_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, baseif_.indexOfWell());
         const auto& gfr = baseif_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, baseif_.indexOfWell());
         const bool use_vfpexplicit = baseif_.useVfpExplicit();
         return baseif_.vfpProperties()->getProd()->bhp(controls.vfp_table_number, aqua, liquid, vapour, baseif_.getTHPConstraint(summaryState), baseif_.getALQ(well_state), wfr, gfr, use_vfpexplicit) - dp;
     }
     else {
         OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER for well " + baseif_.name(), deferred_logger);
     }
}

#define INSTANCE_METHODS(A,...) \
template void WellInterfaceEval<A>:: \
assembleControlEqProd_<__VA_ARGS__>(const WellState&, \
                                    const GroupState&, \
                                    const Schedule&, \
                                    const SummaryState&, \
                                    const Well::ProductionControls&, \
                                    const __VA_ARGS__&, \
                                    const std::vector<__VA_ARGS__>&, \
                                    const std::function<__VA_ARGS__()>&, \
                                    __VA_ARGS__&, \
                                    DeferredLogger&) const; \
template void WellInterfaceEval<A>:: \
assembleControlEqInj_<__VA_ARGS__>(const WellState&, \
                                   const GroupState&, \
                                   const Schedule&, \
                                   const SummaryState&, \
                                   const Well::InjectionControls&, \
                                   const __VA_ARGS__&, \
                                   const __VA_ARGS__&, \
                                   const std::function<__VA_ARGS__()>&, \
                                   __VA_ARGS__&, \
                                   DeferredLogger&) const; \
template __VA_ARGS__ WellInterfaceEval<A>:: \
calculateBhpFromThp<__VA_ARGS__>(const WellState&, \
                                 const std::vector<__VA_ARGS__>&, \
                                 const Well&, \
                                 const SummaryState&, \
                                 const double, \
                                 DeferredLogger&) const;

using FluidSys = BlackOilFluidSystem<double, BlackOilDefaultIndexTraits>;

template class WellInterfaceEval<FluidSys>;

INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,3,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,4,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,5,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,6,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,7,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,8,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,9,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,4u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,5u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,6u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,7u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,8u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,9u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,10u>)

#define INSTANCE_BHP(...) \
template double WellInterfaceEval<__VA_ARGS__>:: \
calculateBhpFromThp<double>(const WellState&, \
                            const std::vector<double>&, \
                            const Well&, \
                            const SummaryState&, \
                            const double, \
                            DeferredLogger&) const;

INSTANCE_BHP(FluidSys)

} // namespace Opm
