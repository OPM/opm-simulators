/*
  Copyright 2020 Equinor ASA.

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
#ifndef OPM_TARGETCALCULATOR_CPP_INCLUDED
#define OPM_TARGETCALCULATOR_CPP_INCLUDED

#include <config.h>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>

#include <cassert>
#include <stdexcept>

namespace Opm::WGHelpers {

template<typename FluidSystem, typename Indices>
TargetCalculator<FluidSystem, Indices>::
TargetCalculator(const Group::ProductionCMode cmode,
                 const std::vector<Scalar>& resv_coeff,
                 const Scalar group_grat_target_from_sales,
                 const std::string& group_name,
                 const GroupState<Scalar>& group_state,
                 const bool use_gpmaint)
    : cmode_(cmode)
    , resv_coeff_(resv_coeff)
    , group_grat_target_from_sales_(group_grat_target_from_sales)
    , group_name_(group_name)
    , group_state_(group_state)
    , use_gpmaint_(use_gpmaint)

{
}

template<typename FluidSystem, typename Indices>
template <typename RateType>
RateType TargetCalculator<FluidSystem, Indices>::calcModeRateFromRates(const RateType* rates) const
{
    switch (cmode_) {
    case Group::ProductionCMode::ORAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::WRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::GRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
        const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::LRAT: {
        // TODO: do need both oil and water activate to use LRAT?
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const int opos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
        const int wpos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
        return rates[opos] + rates[wpos];
    }
    case Group::ProductionCMode::RESV: {
        auto mode_rate = rates[0] * resv_coeff_[0];
        const int num_phases = Indices::numPhases;
        for (int phase = 1; phase < num_phases; ++phase) {
            mode_rate += rates[phase] * resv_coeff_[phase];
        }
        return mode_rate;
    }
    default:
        // Should never be here.
        assert(false);
        return rates[0];
    }
}

template<typename FluidSystem, typename Indices>
typename FluidSystem::Scalar
TargetCalculator<FluidSystem, Indices>::
groupTarget(const std::optional<Group::ProductionControls>& ctrl,
            DeferredLogger& deferred_logger) const
{
    if (!ctrl && !use_gpmaint_) {
        OPM_DEFLOG_THROW(std::logic_error,
                         "Production group " + this->group_name_
                         + "must either have a valid control or use GPMAINT",
                         deferred_logger);
    }
    switch (cmode_) {
    case Group::ProductionCMode::ORAT:
        return ctrl->oil_target;
    case Group::ProductionCMode::WRAT:
        return ctrl->water_target;
    case Group::ProductionCMode::GRAT:
    {
        // gas target may have been adjusted by GCONSALE
        if (group_grat_target_from_sales_ > 0)
            return group_grat_target_from_sales_;

        return ctrl->gas_target;
    }
    case Group::ProductionCMode::LRAT:
        return ctrl->liquid_target;
    case Group::ProductionCMode::RESV:
    {
        if (use_gpmaint_ && this->group_state_.has_gpmaint_target(this->group_name_))
            return this->group_state_.gpmaint_target(this->group_name_);

        return ctrl->resv_target;
    }
    default:
        // Should never be here.
        assert(false);
        return 0.0;
    }
}

template<typename FluidSystem, typename Indices>
GuideRateModel::Target
TargetCalculator<FluidSystem, Indices>::guideTargetMode() const
{
    switch (cmode_) {
    case Group::ProductionCMode::ORAT:
        return GuideRateModel::Target::OIL;
    case Group::ProductionCMode::WRAT:
        return GuideRateModel::Target::WAT;
    case Group::ProductionCMode::GRAT:
        return GuideRateModel::Target::GAS;
    case Group::ProductionCMode::LRAT:
        return GuideRateModel::Target::LIQ;
    case Group::ProductionCMode::RESV:
        return GuideRateModel::Target::RES;
    default:
        // Should never be here.
        assert(false);
        return GuideRateModel::Target::NONE;
    }
}

template<typename FluidSystem, typename Indices>
InjectionTargetCalculator<FluidSystem, Indices>::
InjectionTargetCalculator(const Group::InjectionCMode& cmode,
                          const std::vector<Scalar>& resv_coeff,
                          const std::string& group_name,
                          const Scalar sales_target,
                          const GroupState<Scalar>& group_state,
                          const Phase& injection_phase,
                          const bool use_gpmaint,
                          DeferredLogger& deferred_logger)
    : cmode_(cmode)
    , resv_coeff_(resv_coeff)
    , group_name_(group_name)
    , sales_target_(sales_target)
    , group_state_(group_state)
    , use_gpmaint_(use_gpmaint)

{
    // initialize to avoid warning
    pos_ =  FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    target_ = GuideRateModel::Target::WAT;

    switch (injection_phase) {
    case Phase::WATER: {
        pos_ =  FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
        target_ = GuideRateModel::Target::WAT;
        break;
    }
    case Phase::OIL: {
        pos_ =  FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
        target_ = GuideRateModel::Target::OIL;
        break;
    }
    case Phase::GAS: {
        pos_ =  FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
        target_ = GuideRateModel::Target::GAS;
        break;
    }
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid injection phase in InjectionTargetCalculator",
                         deferred_logger);
    }
}

template<typename FluidSystem, typename Indices>
typename FluidSystem::Scalar
InjectionTargetCalculator<FluidSystem, Indices>::
groupTarget(const std::optional<Group::InjectionControls>& ctrl,
            DeferredLogger& deferred_logger) const
{
    if (!ctrl && !use_gpmaint_) {
        OPM_DEFLOG_THROW(std::logic_error,
                         "Injection group " + this->group_name_
                         + "must either have a valid control or use GPMAINT",
                         deferred_logger);
    }
    switch (cmode_) {
    case Group::InjectionCMode::RATE:
        if (use_gpmaint_ && this->group_state_.has_gpmaint_target(this->group_name_))
            return this->group_state_.gpmaint_target(this->group_name_);

        return ctrl->surface_max_rate;
    case Group::InjectionCMode::RESV:
        if (use_gpmaint_ && this->group_state_.has_gpmaint_target(this->group_name_))
            return this->group_state_.gpmaint_target(this->group_name_) / resv_coeff_[pos_];

        return ctrl->resv_max_rate / resv_coeff_[pos_];
    case Group::InjectionCMode::REIN: {
        Scalar production_rate = this->group_state_.injection_rein_rates(ctrl->reinj_group)[pos_];
        return ctrl->target_reinj_fraction * production_rate;
    }
    case Group::InjectionCMode::VREP: {
        // We use the injection_reservoir_rates directly instead of the reduction rates here to account for the
        // possibility that the group in question has both a VREP control and another injection control for a different phase.
        const std::vector<Scalar>& group_injection_reservoir_rates = this->group_state_.injection_reservoir_rates(this->group_name_);
        Scalar voidage_rate = group_state_.injection_vrep_rate(ctrl->voidage_group) * ctrl->target_void_fraction;
        if (ctrl->phase != Phase::WATER) {
            const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[pos];
        }
        if (ctrl->phase != Phase::OIL) {
            const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[pos];
        }
        if (ctrl->phase != Phase::GAS) {
            const int pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
            voidage_rate -= group_injection_reservoir_rates[pos];
        }
        return voidage_rate / resv_coeff_[pos_];
    }
    case Group::InjectionCMode::SALE: {
        assert(pos_ == FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx) );
        // Gas injection rate = Total gas production rate + gas import rate - gas consumption rate - sales rate;
        // Gas import and consumption is already included in the REIN rates
        Scalar inj_rate = group_state_.injection_rein_rates(this->group_name_)[pos_];
        inj_rate -= sales_target_;
        return inj_rate;
    }
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid Group::InjectionCMode in InjectionTargetCalculator",
                         deferred_logger);
        return 0.0;
    }
}

template<typename FluidSystem, typename Indices>
GuideRateModel::Target
InjectionTargetCalculator<FluidSystem, Indices>::guideTargetMode() const
{
    return target_;
}

/*
#define INSTANTIATE_TARGET_CALCULATOR(T,...) \
    template __VA_ARGS__                     \
    TargetCalculator<T>::calcModeRateFromRates(const __VA_ARGS__* rates) const;

#define INSTANTIATE_TYPE(T)                                       \
    template class TargetCalculator<T>;                           \
    template class InjectionTargetCalculator<T>;                  \
    INSTANTIATE_TARGET_CALCULATOR(T,T)                            \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,3,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,4,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,5,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,6,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,7,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,8,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,9,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,10,0>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,4>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,5>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,6>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,7>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,8>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,9>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,10>) \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,11>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif
*/

} // namespace Opm::WGHelpers
#endif