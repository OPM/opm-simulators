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

namespace Opm {

class DeferredLogger;
template<class Scalar> class GroupState;
template<typename IndexTraits> class PhaseUsageInfo;

namespace WGHelpers {

/// Based on a group control mode, extract or calculate rates, and
/// provide other conveniences.
template<typename Scalar, typename IndexTraits>
class TargetCalculator
{
public:
    TargetCalculator(const Group::ProductionCMode cmode,
                     const PhaseUsageInfo<IndexTraits>& pu,
                     const std::vector<Scalar>& resv_coeff,
                     const Scalar group_grat_target_from_sales,
                     const std::string& group_name,
                     const GroupState<Scalar>& group_state,
                     const bool use_gpmaint);

    template <typename RateType>
    RateType calcModeRateFromRates(const std::vector<RateType>& rates) const
    {
      return calcModeRateFromRates(rates.data());
    }

    template <typename RateType>
    RateType calcModeRateFromRates(const RateType* rates) const;

    Scalar groupTarget(const std::optional<Group::ProductionControls>& ctrl,
                       DeferredLogger& deferred_logger) const;

    GuideRateModel::Target guideTargetMode() const;

private:
    Group::ProductionCMode cmode_;
    const PhaseUsageInfo<IndexTraits>& pu_;
    const std::vector<Scalar>& resv_coeff_;
    const Scalar group_grat_target_from_sales_;
    const std::string& group_name_;
    const GroupState<Scalar>& group_state_;
    bool use_gpmaint_;
};

/// Based on a group control mode, extract or calculate rates, and
/// provide other conveniences.
template<typename Scalar, typename IndexTraits>
class InjectionTargetCalculator
{
public:
    InjectionTargetCalculator(const Group::InjectionCMode& cmode,
                              const PhaseUsageInfo<IndexTraits>& pu,
                              const std::vector<Scalar>& resv_coeff,
                              const std::string& group_name,
                              const Scalar sales_target,
                              const GroupState<Scalar>& group_state,
                              const Phase& injection_phase,
                              const bool use_gpmaint,
                              DeferredLogger& deferred_logger);

    template <typename RateVec>
    auto calcModeRateFromRates(const RateVec& rates) const
    {
        return rates[pos_];
    }

    Scalar groupTarget(const std::optional<Group::InjectionControls>& ctrl,
                       DeferredLogger& deferred_logger) const;

    GuideRateModel::Target guideTargetMode() const;

private:
    Group::InjectionCMode cmode_;
    const PhaseUsageInfo<IndexTraits>& pu_;
    const std::vector<Scalar>& resv_coeff_;
    const std::string& group_name_;
    Scalar sales_target_;
    const GroupState<Scalar>& group_state_;
    bool use_gpmaint_;
    int pos_;
    GuideRateModel::Target target_;
};

} // namespace WGHelpers

} // namespace Opm

#endif
