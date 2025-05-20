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

    //std::cout << name << " " << always_included_child << " " << parent_group.name() << " " << my_guide_rate << " " << total_guide_rate << " " << num_active_groups << std::endl;
    if (num_active_groups == 1)
        return 1.0;

    if (total_guide_rate == 0 ) {
        std::cout << name << " total_guide_rate is zero " << my_guide_rate << " " << num_active_groups << " " << always_included_child <<  std::endl;
        // if the total guide rate is zero (for instance due to netv = 0) we use the potentials
        // to distribute the group rate
        always_use_potentials = true;

        const Scalar my_pot = guideRate(name, always_included_child, always_use_potentials);
        const Scalar my_total_pot = guideRateSum(parent_group, always_included_child, always_use_potentials).first;
        std::cout << my_pot << " " << my_total_pot << std::endl;
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

    if (true)
    {
        const auto& gs_rates = is_producer_? this->group_state_.prod_guide_rates(group.name()) : this->group_state_.inj_guide_rates(group.name(), injection_phase_);
        auto group_cont_wells = is_producer_? this->group_state_.number_of_wells_under_this_control(group.name()) : this->group_state_.number_of_wells_under_this_inj_control(group.name(), injection_phase_) ;
        const auto& next_subgroup_with_guiderate = is_producer_? this->group_state_.sub_group_with_guiderate(group.name()) : this->group_state_.sub_group_inj_with_guiderate(group.name(), injection_phase_);
        double child_rate = 0.0;        
        bool isWell = schedule_.hasWell(always_included_child, report_step_);
        //std::cout << is_producer_ << " " << gs_rates << " " << group_cont_wells << " " << always_included_child << " " << group.name() << " " << isWell << " " << always_use_potentials << std::endl;
        if (isWell) {
            bool already_included = is_producer_? well_state_.isProductionGrup(always_included_child) : well_state_.isInjectionGrup(always_included_child);
            if (!already_included) {
                group_cont_wells++;
                //const auto chain = WellGroupHelpers<Scalar>::groupChainTopBot(always_included_child, group.name(), schedule_, report_step_);
                if (group_cont_wells > 0)
                    child_rate = guideRate(always_included_child, always_included_child, always_use_potentials);
                else {
                    const auto& well = schedule_.getWell(always_included_child, report_step_);
                    bool stop = false;
                    auto gr_name = well.groupName();
                    while (!stop) {
                        auto it = std::find(next_subgroup_with_guiderate.begin(), next_subgroup_with_guiderate.end(), gr_name);
                        if (it != next_subgroup_with_guiderate.end()) {
                            child_rate = guideRate(gr_name, always_included_child, always_use_potentials);
                            break;
                        }
                        const auto& group2 = schedule_.getGroup(gr_name, report_step_);
                        if (gr_name == group.name())
                            stop = true;
                        
                        gr_name = group2.parent();
                    }
                }
            }
        }
        return {gs_rates + child_rate, group_cont_wells};
    }




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

    double child_rate = 0.0;
    for (const std::string& child_well : group.wells()) {
        bool included = (child_well == always_included_child);
        if (is_producer_) {
            included |= well_state_.isProductionGrup(child_well);
        } else {
            included |= well_state_.isInjectionGrup(child_well);
        }
        if (included) {
            number_of_included_well_or_groups++;
            std::cout << child_well << " " << always_included_child << " " << always_use_potentials << " " << guideRate(child_well, always_included_child, always_use_potentials) << std::endl;
            total_guide_rate += guideRate(child_well, always_included_child, always_use_potentials);
        }
        //if (child_well == always_included_child && !well_state_.isProductionGrup(child_well))
        //child_rate = guideRate(always_included_child, always_included_child, always_use_potentials);

    }
    const auto& gs_rates = is_producer_? this->group_state_.prod_guide_rates(group.name()) : this->group_state_.inj_guide_rates(group.name(), injection_phase_);
    auto group_cont_wells = is_producer_? this->group_state_.number_of_wells_under_this_control(group.name()) : this->group_state_.number_of_wells_under_this_inj_control(group.name(), injection_phase_) ;
    const auto& next_subgroup_with_guiderate = is_producer_? this->group_state_.sub_group_with_guiderate(group.name()) : this->group_state_.sub_group_inj_with_guiderate(group.name(), injection_phase_);
    child_rate = 0.0;        
    bool isWell = schedule_.hasWell(always_included_child, report_step_);
    //std::cout << is_producer_ << " " << gs_rates << " " << group_cont_wells << " " << always_included_child << " " << group.name() << " " << isWell << " " << always_use_potentials << std::endl;
    if (isWell) {
        bool already_included = is_producer_? well_state_.isProductionGrup(always_included_child) : well_state_.isInjectionGrup(always_included_child);
        if (!already_included) {
            //const auto chain = WellGroupHelpers<Scalar>::groupChainTopBot(always_included_child, group.name(), schedule_, report_step_);
            if (group_cont_wells > 0)
                child_rate = guideRate(always_included_child, always_included_child, always_use_potentials);
            else {
                const auto& well = schedule_.getWell(always_included_child, report_step_);
                bool stop = false;
                auto gr_name = well.groupName();
                while (!stop) {
                    auto it = std::find(next_subgroup_with_guiderate.begin(), next_subgroup_with_guiderate.end(), gr_name);
                    if (it != next_subgroup_with_guiderate.end()) {
                        child_rate = guideRate(gr_name, always_included_child, always_use_potentials);
                        break;
                    }
                    const auto& group2 = schedule_.getGroup(gr_name, report_step_);
                    if (gr_name == group.name())
                        stop = true;
                    
                    gr_name = group2.parent();
                }
            }
        }
    }

   //else if (this->group_state_.has_prod_guide_rates(always_included_child) && always_included_child != group.name()) 
            //child_rate = this->group_state_.prod_guide_rates(always_included_child);

    //if (hasWell && !well_state_.isProductionGrup(always_included_child)) {
    //gs_rates += child_rate;
    //}
    //if (group_cont_wells == 0)
    //    gs_rates = 0.0;

    if (is_producer_ && std::abs(total_guide_rate - gs_rates - child_rate) > 0.1)
        std::cout << "BUG " << group.name() << " " << group_cont_wells << " " << gs_rates << " " << child_rate << " " << total_guide_rate << " " << always_included_child << std::endl;

    return {total_guide_rate, number_of_included_well_or_groups};
}

template<class Scalar>
Scalar FractionCalculator<Scalar>::
guideRate(const std::string& name,
          const std::string& always_included_child,
          const bool always_use_potentials)
{
    if (schedule_.hasWell(name, report_step_)) {
        if (guide_rate_->has(name) || guide_rate_->hasPotentials(name)) {
            return guide_rate_->get(name, target_, WellGroupHelpers<Scalar>::getWellRateVector(well_state_, pu_, name));
        } else {
            return 0.0;
        }
    } else {
        auto group_cont_wells = is_producer_? this->group_state_.number_of_wells_under_this_control(name) : this->group_state_.number_of_wells_under_this_inj_control(name, injection_phase_) ;
        if (schedule_.hasWell(always_included_child, report_step_)) {
            bool already_included = is_producer_? well_state_.isProductionGrup(always_included_child) : well_state_.isInjectionGrup(always_included_child);
            if (!already_included) 
                group_cont_wells++;
        }
        if (false ) {
            group_cont_wells = groupControlledWells(name, always_included_child) > 0;
        }
        //if (group_cont_wells != groupControlledWells(name, always_included_child) )
        //    std::cout << "group_cont_wells "<< always_included_child << " " << name << " " << group_cont_wells <<  " " << groupControlledWells(name, always_included_child) << std::endl;

        if (group_cont_wells) {
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
