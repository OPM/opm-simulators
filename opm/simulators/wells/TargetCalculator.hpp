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

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>

#include <optional>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GroupState;
struct PhaseUsage;

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
                         const double group_grat_target_from_sales,
                         const std::string& group_name,
                         const GroupState& group_state,
                         const bool use_gpmaint);

        template <typename RateType>
        RateType calcModeRateFromRates(const std::vector<RateType>& rates) const
        {
          return calcModeRateFromRates(rates.data());
        }

        template <typename RateType>
        RateType calcModeRateFromRates(const RateType* rates) const;

        double groupTarget(const std::optional<Group::ProductionControls>& ctrl, Opm::DeferredLogger& deferred_logger) const;

        GuideRateModel::Target guideTargetMode() const;

    private:
        Group::ProductionCMode cmode_;
        const PhaseUsage& pu_;
        const std::vector<double>& resv_coeff_;
        const double group_grat_target_from_sales_;
        const std::string& group_name_;
        const GroupState& group_state_;
        bool use_gpmaint_;
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
                                  const GroupState& group_state,
                                  const Phase& injection_phase,
                                  const bool use_gpmaint,
                                  DeferredLogger& deferred_logger);

        template <typename RateVec>
        auto calcModeRateFromRates(const RateVec& rates) const
        {
            return rates[pos_];
        }

        double groupTarget(const std::optional<Group::InjectionControls>& ctrl, Opm::DeferredLogger& deferred_logger) const;

        GuideRateModel::Target guideTargetMode() const;

    private:
        Group::InjectionCMode cmode_;
        const PhaseUsage& pu_;
        const std::vector<double>& resv_coeff_;
        const std::string& group_name_;
        double sales_target_;
        const GroupState& group_state_;
        bool use_gpmaint_;
        int pos_;
        GuideRateModel::Target target_;
    };
} // namespace WellGroupHelpers

} // namespace Opm

#endif
