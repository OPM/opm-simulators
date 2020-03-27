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
                         const std::vector<double>& resv_coeff)
            : cmode_(cmode)
            , pu_(pu)
            , resv_coeff_(resv_coeff)
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
                return ctrl.gas_target;
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
    };

} // namespace WellGroupHelpers

} // namespace Opm

#endif
