/*
  Copyright 2019 Norce.
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

#include <config.h>
#include <opm/simulators/wells/FractionCalculator.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>

namespace Opm::WGHelpers {

template<class Scalar>
FractionCalculator<Scalar>::
FractionCalculator(const Schedule& schedule,
                   const WellState<Scalar>& well_state,
                   const GroupState<Scalar>& group_state,
                   const SummaryState& summary_state,
                   const int report_step,
                   const GuideRate* guide_rate,
                   const GuideRateModel::Target target,
                   const PhaseUsage& pu,
                   const bool is_producer,
                   const Phase injection_phase)
    : schedule_(schedule)
    , well_state_(well_state)
    , group_state_(group_state)
    , summary_state_(summary_state)
    , report_step_(report_step)
    , guide_rate_(guide_rate)
    , target_(target)
    , pu_(pu)
    , is_producer_(is_producer)
    , injection_phase_(injection_phase)
{
}

template<class Scalar>
Scalar FractionCalculator<Scalar>::
fraction(const std::string& name,
         const std::string& control_group_name,
         const bool always_include_this)
{
    Scalar fraction = 1.0;
    std::string current = name;
    while (current != control_group_name) {
        fraction *= localFraction(current, always_include_this ? name : "");
        current = parent(current);
    }
    return fraction;
}

template<class Scalar>
Scalar FractionCalculator<Scalar>::
localFraction(const std::string& name,
              const std::string& always_included_child)
{
    bool always_use_potentials = false;
    const Scalar my_guide_rate = guideRate(name, always_included_child, always_use_potentials);

    const Group& parent_group = schedule_.getGroup(parent(name), report_step_);
    const auto [total_guide_rate, num_active_groups] = guideRateSum(parent_group, always_included_child, always_use_potentials);

    // the group/well "name" is the only active group/well we therefore return 1 as the fraction
    // even though my_guide_rate may be zero
    if (num_active_groups == 1)
        return 1.0;

    if (total_guide_rate == 0 ) {
        // if the total guide rate is zero (for instance due to netv = 0) we use the potentials
        // to distribute the group rate
        always_use_potentials = true;
        const Scalar my_pot = guideRate(name, always_included_child, always_use_potentials);
        const Scalar my_total_pot = guideRateSum(parent_group, always_included_child, always_use_potentials).first;
        return my_pot / my_total_pot;
    }
    return my_guide_rate / total_guide_rate;
}

template<class Scalar>
std::string FractionCalculator<Scalar>::
parent(const std::string& name)
{
    if (schedule_.hasWell(name)) {
        return schedule_.getWell(name, report_step_).groupName();
    } else {
        return schedule_.getGroup(name, report_step_).parent();
    }
}

template<class Scalar>
std::pair<Scalar, int> FractionCalculator<Scalar>::
guideRateSum(const Group& group,
             const std::string& always_included_child,
             const bool always_use_potentials)
{
    Scalar total_guide_rate = 0.0;
    int number_of_included_well_or_groups = 0;
    for (const std::string& child_group : group.groups()) {
        bool included = (child_group == always_included_child);
        if (is_producer_) {
            const auto ctrl = this->group_state_.production_control(child_group);
            included |= (ctrl == Group::ProductionCMode::FLD) ||
                        (ctrl == Group::ProductionCMode::NONE);
        } else {
            const auto ctrl = this->group_state_.injection_control(child_group,
                                                                   this->injection_phase_);
            included |= (ctrl == Group::InjectionCMode::FLD) ||
                        (ctrl == Group::InjectionCMode::NONE);
        }
        if (included) {
            if (groupControlledWells(child_group, always_included_child) > 0) {
                number_of_included_well_or_groups++;
                total_guide_rate += guideRate(child_group, always_included_child, always_use_potentials);
            }
        }
    }

    for (const std::string& child_well : group.wells()) {
        bool included = (child_well == always_included_child);
        if (is_producer_) {
            included |= well_state_.isProductionGrup(child_well);
        } else {
            included |= well_state_.isInjectionGrup(child_well);
        }
        if (included) {
            number_of_included_well_or_groups++;
            total_guide_rate += guideRate(child_well, always_included_child, always_use_potentials);
        }
    }
    return {total_guide_rate, number_of_included_well_or_groups};
}

template<class Scalar>
Scalar FractionCalculator<Scalar>::
guideRate(const std::string& name,
          const std::string& always_included_child,
          const bool always_use_potentials)
{
    if (schedule_.hasWell(name, report_step_)) {
        return WellGroupHelpers<Scalar>::getGuideRate(name, schedule_, well_state_, group_state_,
                                                      report_step_, guide_rate_, target_, pu_);
    } else {
        if (groupControlledWells(name, always_included_child) > 0) {
            if (is_producer_ && guide_rate_->has(name) && !always_use_potentials) {
                return guide_rate_->get(name, target_, getGroupRateVector(name));
            } else if (!is_producer_ && guide_rate_->has(name, injection_phase_) && !always_use_potentials) {
                return guide_rate_->get(name, injection_phase_);
            } else {
                // We are a group, with default guide rate.
                // Compute guide rate by accumulating our children's guide rates.
                const Group& group = schedule_.getGroup(name, report_step_);
                const double eff = group.getGroupEfficiencyFactor();
                return eff * guideRateSum(group, always_included_child, always_use_potentials).first;
            }
        } else {
            // No group-controlled subordinate wells.
            return 0.0;
        }
    }
}

template<class Scalar>
int FractionCalculator<Scalar>::
groupControlledWells(const std::string& group_name,
                     const std::string& always_included_child)
{
    return WellGroupHelpers<Scalar>::groupControlledWells(schedule_,
                                                          well_state_,
                                                          this->group_state_,
                                                          this->summary_state_,
                                                          this->guide_rate_,
                                                          report_step_,
                                                          group_name,
                                                          always_included_child,
                                                          is_producer_,
                                                          injection_phase_);
}

template<class Scalar>
GuideRate::RateVector FractionCalculator<Scalar>::
getGroupRateVector(const std::string& group_name)
{
    assert(is_producer_);
    return WellGroupHelpers<Scalar>::getProductionGroupRateVector(this->group_state_,
                                                                  this->pu_,
                                                                  group_name);
}

template class FractionCalculator<double>;

#if FLOW_INSTANTIATE_FLOAT
template class FractionCalculator<float>;
#endif

} // namespace Opm::WGHelpers
