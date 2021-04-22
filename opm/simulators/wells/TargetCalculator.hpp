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


#ifndef OPM_TARGETCALCULATOR_HEADER_INCLUDED
#define OPM_TARGETCALCULATOR_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace Opm
{

namespace WellGroupHelpers
{

    /// Based on a group control mode, extract or calculate rates, and
    /// provide other conveniences.
    class TargetCalculator
    {
    public:
        TargetCalculator(const Group::ProductionCMode cmode,
                         const PhaseUsage& pu,
                         const std::vector<double>& resv_coeff,
                         const double group_grat_target_from_sales)
            : cmode_(cmode)
            , pu_(pu)
            , resv_coeff_(resv_coeff)
            , group_grat_target_from_sales_(group_grat_target_from_sales)
        {

        }

        template <typename RateVec>
        auto calcModeRateFromRates(const RateVec& rates) const
        {
            // ElemType is just the plain element type of the rates container,
            // without any reference, const or volatile modifiers.
            using ElemType = std::remove_cv_t<std::remove_reference_t<decltype(rates[0])>>;
            switch (cmode_) {
            case Group::ProductionCMode::ORAT: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                const int pos = pu_.phase_pos[BlackoilPhases::Liquid];
                return rates[pos];
            }
            case Group::ProductionCMode::WRAT: {
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                const int pos = pu_.phase_pos[BlackoilPhases::Aqua];
                return rates[pos];
            }
            case Group::ProductionCMode::GRAT: {
                assert(pu_.phase_used[BlackoilPhases::Vapour]);
                const int pos = pu_.phase_pos[BlackoilPhases::Vapour];
                return rates[pos];
            }
            case Group::ProductionCMode::LRAT: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                const int opos = pu_.phase_pos[BlackoilPhases::Liquid];
                const int wpos = pu_.phase_pos[BlackoilPhases::Aqua];
                return rates[opos] + rates[wpos];
            }
            case Group::ProductionCMode::RESV: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                assert(pu_.phase_used[BlackoilPhases::Vapour]);
                ElemType mode_rate = zero<ElemType>();
                for (int phase = 0; phase < pu_.num_phases; ++phase) {
                    mode_rate += rates[phase] * resv_coeff_[phase];
                }
                return mode_rate;
            }
            default:
                // Should never be here.
                assert(false);
                return zero<ElemType>();
            }
        }

        double groupTarget(const Group::ProductionControls ctrl) const
        {
            switch (cmode_) {
            case Group::ProductionCMode::ORAT:
                return ctrl.oil_target;
            case Group::ProductionCMode::WRAT:
                return ctrl.water_target;
            case Group::ProductionCMode::GRAT:
            {
                // gas target may have been adjusted by GCONSALE
                if ( group_grat_target_from_sales_ > 0)
                    return group_grat_target_from_sales_;

                return ctrl.gas_target;
            }
            case Group::ProductionCMode::LRAT:
                return ctrl.liquid_target;
            case Group::ProductionCMode::RESV:
                return ctrl.resv_target;
            default:
                // Should never be here.
                assert(false);
                return 0.0;
            }
        }

        GuideRateModel::Target guideTargetMode() const
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

    private:
        template <typename ElemType>
        static ElemType zero()
        {
            // This is for Evaluation types.
            ElemType x;
            x = 0.0;
            return x;
        }
        Group::ProductionCMode cmode_;
        const PhaseUsage& pu_;
        const std::vector<double>& resv_coeff_;
        const double group_grat_target_from_sales_;
    };


    /// Based on a group control mode, extract or calculate rates, and
    /// provide other conveniences.
    class InjectionTargetCalculator
    {
    public:
        InjectionTargetCalculator(const Group::InjectionCMode& cmode,
                                  const PhaseUsage& pu,
                                  const std::vector<double>& resv_coeff,
                                  const std::string& group_name,
                                  const double sales_target,
                                  const WellStateFullyImplicitBlackoil& well_state,
                                  const GroupState& group_state,
                                  const Phase& injection_phase,
                                  DeferredLogger& deferred_logger)
            : cmode_(cmode)
            , pu_(pu)
            , resv_coeff_(resv_coeff)
            , group_name_(group_name)
            , sales_target_(sales_target)
            , well_state_(well_state)
            , group_state_(group_state)
        {
            // initialize to avoid warning
            pos_ = pu.phase_pos[BlackoilPhases::Aqua];
            target_ = GuideRateModel::Target::WAT;

            switch (injection_phase) {
            case Phase::WATER: {
                pos_ = pu.phase_pos[BlackoilPhases::Aqua];
                target_ = GuideRateModel::Target::WAT;
                break;
            }
            case Phase::OIL: {
                pos_ = pu.phase_pos[BlackoilPhases::Liquid];
                target_ = GuideRateModel::Target::OIL;
                break;
            }
            case Phase::GAS: {
                pos_ = pu.phase_pos[BlackoilPhases::Vapour];
                target_ = GuideRateModel::Target::GAS;
                break;
            }
            default:
                OPM_DEFLOG_THROW(std::logic_error,
                                 "Invalid injection phase in InjectionTargetCalculator",
                                 deferred_logger);
            }
        }


        template <typename RateVec>
        auto calcModeRateFromRates(const RateVec& rates) const
        {
            return rates[pos_];
        }

        double groupTarget(const Group::InjectionControls& ctrl, Opm::DeferredLogger& deferred_logger) const
        {
            switch (cmode_) {
            case Group::InjectionCMode::RATE:
                return ctrl.surface_max_rate;
            case Group::InjectionCMode::RESV:
                return ctrl.resv_max_rate;
            case Group::InjectionCMode::REIN: {
                double production_rate = this->group_state_.injection_rein_rates(ctrl.reinj_group)[pos_];
                return ctrl.target_reinj_fraction * production_rate;
            }
            case Group::InjectionCMode::VREP: {
                const std::vector<double>& group_injection_reductions = this->group_state_.injection_reduction_rates(this->group_name_);
                double voidage_rate = group_state_.injection_vrep_rate(ctrl.voidage_group) * ctrl.target_void_fraction;
                double inj_reduction = 0.0;
                if (ctrl.phase != Phase::WATER)
                    inj_reduction += group_injection_reductions[pu_.phase_pos[BlackoilPhases::Aqua]]
                            * resv_coeff_[pu_.phase_pos[BlackoilPhases::Aqua]];
                if (ctrl.phase != Phase::OIL)
                    inj_reduction += group_injection_reductions[pu_.phase_pos[BlackoilPhases::Liquid]]
                            * resv_coeff_[pu_.phase_pos[BlackoilPhases::Liquid]];
                if (ctrl.phase != Phase::GAS)
                    inj_reduction += group_injection_reductions[pu_.phase_pos[BlackoilPhases::Vapour]]
                            * resv_coeff_[pu_.phase_pos[BlackoilPhases::Vapour]];
                voidage_rate -= inj_reduction;
                return voidage_rate / resv_coeff_[pos_];
            }
            case Group::InjectionCMode::SALE: {
                assert(pos_ == pu_.phase_pos[BlackoilPhases::Vapour]);
                // Gas injection rate = Total gas production rate + gas import rate - gas consumption rate - sales rate;
                // Gas import and consumption is already included in the REIN rates
                double inj_rate = group_state_.injection_rein_rates(this->group_name_)[pos_];
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

        GuideRateModel::Target guideTargetMode() const
        {
            return target_;
        }

    private:
        template <typename ElemType>
        static ElemType zero()
        {
            // This is for Evaluation types.
            ElemType x;
            x = 0.0;
            return x;
        }
        Group::InjectionCMode cmode_;
        const PhaseUsage& pu_;
        const std::vector<double>& resv_coeff_;
        const std::string& group_name_;
        double sales_target_;
        const WellStateFullyImplicitBlackoil& well_state_;
        const GroupState& group_state_;
        int pos_;
        GuideRateModel::Target target_;

    };
} // namespace WellGroupHelpers

} // namespace Opm

#endif
