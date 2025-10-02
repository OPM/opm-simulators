/*
  Copyright 2021 Equinor

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <iterator>

#include <opm/json/JsonObject.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/wells/GroupState.hpp>


namespace Opm {

template<class Scalar>
GroupState<Scalar>::GroupState(std::size_t np) :
    num_phases(np)
{}

template<class Scalar>
GroupState<Scalar> GroupState<Scalar>::serializationTestObject()
{
    GroupState result(3);
    result.m_production_rates = {{"test1", {1.0, 2.0}}};
    result.m_network_leaf_node_production_rates={{"test1", {1.0, 20}}};
    result.production_controls = {{"test2", Group::ProductionCMode::LRAT}};
    result.prod_red_rates = {{"test3", {3.0, 4.0, 5.0}}};
    result.inj_red_rates = {{"test4", {6.0, 7.0}}};
    result.inj_surface_rates = {{"test5", {8.0}}};
    result.inj_resv_rates = {{"test6", {9.0, 10.0}}};
    result.inj_rein_rates = {{"test7", {11.0}}};
    result.inj_vrep_rate = {{"test8", 12.0}, {"test9", 13.0}};
    result.m_grat_sales_target = {{"test10", 14.0}};
    result.m_gpmaint_target = {{"test11", 15.0}};
    result.injection_controls = {{{Phase::FOAM, "test12"}, Group::InjectionCMode::REIN}};
    result.gpmaint_state.add("foo", GPMaint::State::serializationTestObject());
    result.m_gconsump_rates = {{"testA", {0.2, 0.1}}};
    result.m_number_of_wells_under_group_control = {{"test13", 3}};
    result.m_number_of_wells_under_inj_group_control = {{ {Phase::WATER, "test14"}, 3}};

    return result;
}

template<class Scalar>
bool GroupState<Scalar>::operator==(const GroupState& other) const
{
    return this->m_production_rates == other.m_production_rates &&
           this->m_network_leaf_node_production_rates == other.m_network_leaf_node_production_rates &&
           this->production_controls == other.production_controls &&
           this->prod_red_rates == other.prod_red_rates &&
           this->inj_red_rates == other.inj_red_rates &&
           this->inj_resv_rates == other.inj_resv_rates &&
           this->inj_rein_rates == other.inj_rein_rates &&
           this->inj_vrep_rate == other.inj_vrep_rate &&
           this->inj_surface_rates == other.inj_surface_rates &&
           this->m_grat_sales_target == other.m_grat_sales_target &&
           this->injection_controls == other.injection_controls &&
           this->m_number_of_wells_under_group_control == m_number_of_wells_under_group_control &&
           this->m_number_of_wells_under_inj_group_control == m_number_of_wells_under_inj_group_control &&
           this->gpmaint_state == other.gpmaint_state &&
           this->m_gconsump_rates == other.m_gconsump_rates;
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::has_production_rates(const std::string& gname) const
{
    auto group_iter = this->m_production_rates.find(gname);
    return (group_iter != this->m_production_rates.end());
}

template<class Scalar>
void GroupState<Scalar>::update_production_rates(const std::string& gname,
                                                 const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->m_production_rates[gname] = rates;
}

template<class Scalar>
void GroupState<Scalar>::update_network_leaf_node_production_rates(const std::string& gname,
                                                 const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->m_network_leaf_node_production_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::production_rates(const std::string& gname) const
{
    auto group_iter = this->m_production_rates.find(gname);
    if (group_iter == this->m_production_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::network_leaf_node_production_rates(const std::string& gname) const
{
    auto group_iter = this->m_network_leaf_node_production_rates.find(gname);
    if (group_iter == this->m_network_leaf_node_production_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
GroupState::update_well_group_thp(const std::string& gname, const double& thp)
{
    this->group_thp[gname] = thp;
}

template<class Scalar>
Scalar GroupState<Scalar>::
GroupState::well_group_thp(const std::string& gname) const
{
    auto group_iter = this->group_thp.find(gname);
    if (group_iter == this->group_thp.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

template<class Scalar>
bool GroupState<Scalar>::
GroupState::is_autochoke_group(const std::string& gname) const
{
    return (this->group_thp.count(gname) > 0);
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_production_reduction_rates(const std::string& gname) const
{
    auto group_iter = this->prod_red_rates.find(gname);
    return (group_iter != this->prod_red_rates.end());
}

template<class Scalar>
void GroupState<Scalar>::
update_production_reduction_rates(const std::string& gname,
                                  const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->prod_red_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::production_reduction_rates(const std::string& gname) const
{
    auto group_iter = this->prod_red_rates.find(gname);
    if (group_iter == this->prod_red_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_injection_reduction_rates(const std::string& gname) const
{
    auto group_iter = this->inj_red_rates.find(gname);
    return (group_iter != this->inj_red_rates.end());
}

template<class Scalar>
void GroupState<Scalar>::
update_injection_reduction_rates(const std::string& gname,
                                 const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_red_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::injection_reduction_rates(const std::string& gname) const
{
    auto group_iter = this->inj_red_rates.find(gname);
    if (group_iter == this->inj_red_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}
//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_injection_surface_rates(const std::string& gname) const
{
    auto group_iter = this->inj_surface_rates.find(gname);
    return (group_iter != this->inj_surface_rates.end());
}

template<class Scalar>
void GroupState<Scalar>::
update_injection_surface_rates(const std::string& gname,
                               const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_surface_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::
injection_surface_rates(const std::string& gname) const
{
    auto group_iter = this->inj_surface_rates.find(gname);
    if (group_iter == this->inj_surface_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_injection_reservoir_rates(const std::string& gname) const
{
    auto group_iter = this->inj_resv_rates.find(gname);
    return (group_iter != this->inj_resv_rates.end());
}

template<class Scalar>
void GroupState<Scalar>::
update_injection_reservoir_rates(const std::string& gname,
                                 const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_resv_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::injection_reservoir_rates(const std::string& gname) const
{
    auto group_iter = this->inj_resv_rates.find(gname);
    if (group_iter == this->inj_resv_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
update_injection_rein_rates(const std::string& gname,
                            const std::vector<Scalar>& rates)
{
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_rein_rates[gname] = rates;
}

template<class Scalar>
const std::vector<Scalar>&
GroupState<Scalar>::
injection_rein_rates(const std::string& gname) const
{
    auto group_iter = this->inj_rein_rates.find(gname);
    if (group_iter == this->inj_rein_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
update_injection_vrep_rate(const std::string& gname, Scalar rate)
{
    this->inj_vrep_rate[gname] = rate;
}

template<class Scalar>
Scalar GroupState<Scalar>::
injection_vrep_rate(const std::string& gname) const
{
    auto group_iter = this->inj_vrep_rate.find(gname);
    if (group_iter == this->inj_vrep_rate.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
update_grat_sales_target(const std::string& gname, Scalar target)
{
    this->m_grat_sales_target[gname] = target;
}

template<class Scalar>
Scalar GroupState<Scalar>::
grat_sales_target(const std::string& gname) const
{
    auto group_iter = this->m_grat_sales_target.find(gname);
    if (group_iter == this->m_grat_sales_target.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

template<class Scalar>
bool GroupState<Scalar>::
has_grat_sales_target(const std::string& gname) const
{
    return (this->m_grat_sales_target.count(gname) > 0);
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_production_control(const std::string& gname) const
{
    auto group_iter = this->production_controls.find(gname);
    if (group_iter == this->production_controls.end())
        return false;

    return true;
}

template<class Scalar>
void GroupState<Scalar>::
production_control(const std::string& gname,
                   Group::ProductionCMode cmode)
{
    this->production_controls[gname] = cmode;
}

template<class Scalar>
Group::ProductionCMode
GroupState<Scalar>::production_control(const std::string& gname) const
{
    auto group_iter = this->production_controls.find(gname);
    if (group_iter == this->production_controls.end())
        throw std::logic_error("Could not find any control for production group: " + gname);

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
bool GroupState<Scalar>::
has_injection_control(const std::string& gname, Phase phase) const
{
    return this->injection_controls.count(std::make_pair(phase, gname)) > 0;
}

template<class Scalar>
void GroupState<Scalar>::
injection_control(const std::string& gname,
                  Phase phase, Group::InjectionCMode cmode)
{
    this->injection_controls[ std::make_pair(phase, gname) ] = cmode;
}

template<class Scalar>
Group::InjectionCMode
GroupState<Scalar>::
injection_control(const std::string& gname, Phase phase) const
{
    auto key = std::make_pair(phase, gname);
    auto group_iter = this->injection_controls.find( key );
    if (group_iter == this->injection_controls.end())
        throw std::logic_error("Could not find ontrol for injection group: " + gname);

    return group_iter->second;
}



//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
update_number_of_wells_under_group_control(const std::string& gname, int number)
{
    this->m_number_of_wells_under_group_control[gname] = number;
}
template<class Scalar>
int GroupState<Scalar>::
number_of_wells_under_group_control(const std::string& gname) const
{
    auto group_iter = this->m_number_of_wells_under_group_control.find(gname);
    if (group_iter == this->m_number_of_wells_under_group_control.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------


template<class Scalar>
void GroupState<Scalar>::
update_number_of_wells_under_inj_group_control(const std::string& gname, Phase phase, int number)
{
    this->m_number_of_wells_under_inj_group_control[{phase, gname}] = number;
}

template<class Scalar>
int GroupState<Scalar>::
number_of_wells_under_inj_group_control(const std::string& gname, Phase phase) const
{
    auto group_iter = this->m_number_of_wells_under_inj_group_control.find({phase, gname});
    if (group_iter == this->m_number_of_wells_under_inj_group_control.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

template<class Scalar>
GPMaint::State& GroupState<Scalar>::gpmaint(const std::string& gname)
{
    if (!this->gpmaint_state.has(gname))
        this->gpmaint_state.add(gname, GPMaint::State{});
    return this->gpmaint_state[gname];
}


//-------------------------------------------------------------------------

template<class Scalar>
void GroupState<Scalar>::
update_gpmaint_target(const std::string& gname, Scalar target)
{
    this->m_gpmaint_target[gname] = target;
}

template<class Scalar>
Scalar GroupState<Scalar>::gpmaint_target(const std::string& gname) const
{
    auto group_iter = this->m_gpmaint_target.find(gname);
    if (group_iter == this->m_gpmaint_target.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

template<class Scalar>
bool GroupState<Scalar>::
has_gpmaint_target(const std::string& gname) const
{
    return (this->m_gpmaint_target.count(gname) > 0);
}

template<class Scalar>
void GroupState<Scalar>::
update_gconsump(const Schedule& schedule, const int report_step, const SummaryState& summary_state) {
    this->m_gconsump_rates.clear();
    const auto& sched_state = schedule[report_step];
    const auto& gconsump = sched_state.gconsump();
    auto gcr_recursive =
        [this, &sched_state, &gconsump, &summary_state](auto self,
                                                        const std::string& group_name,
                                                        std::pair<Scalar, Scalar>& rates,
                                                        const double parent_gefac = 1.0) -> bool
        {
            // If group already has been computed, update parent rates and return true
            const auto it = this->m_gconsump_rates.find(group_name);
            if (it != this->m_gconsump_rates.end()) {
                rates.first += static_cast<Scalar>(it->second.first * parent_gefac);
                rates.second += static_cast<Scalar>(it->second.second * parent_gefac);
                return true;
            }

            // Accumulate from sub-groups and keep track of any updates in 'has_values'
            bool has_values = false;
            if (sched_state.groups.has(group_name)) {
                for (const auto& child_gname : sched_state.groups(group_name).groups()) {
                    const auto gefac = sched_state.groups(child_gname).getGroupEfficiencyFactor();
                    has_values = self(self, child_gname, rates, gefac) || has_values;
                }
            }

            // Add consumption/import rates at current level
            if (gconsump.has(group_name)) {
                const auto& group_gc = gconsump.get(group_name, summary_state);
                rates.first += static_cast<Scalar>(group_gc.consumption_rate);
                rates.second += static_cast<Scalar>(group_gc.import_rate);
                has_values = true;
            }

            // Update map if values are set
            if (has_values) this->m_gconsump_rates.insert_or_assign(group_name, rates);
            return has_values;
        };

    auto rates = std::pair { Scalar{0}, Scalar{0} };
    gcr_recursive(gcr_recursive, "FIELD", rates);
}

template<class Scalar>
const std::pair<Scalar, Scalar>& GroupState<Scalar>::
gconsump_rates(const std::string& gname) const {
    const auto it = this->m_gconsump_rates.find(gname);
    if (it != this->m_gconsump_rates.end()) {
        return it->second;
    }
    return zero_pair;
}

template<class Scalar>
void GroupState<Scalar>::
update_group_production_potential(const std::string& gname, Scalar oil_rate,
                                  Scalar gas_rate, Scalar water_rate)
{
    auto [it, inserted] = production_group_potentials.try_emplace(
        gname, oil_rate, gas_rate, water_rate
    );
    if (!inserted) {
        it->second = GroupPotential{oil_rate, gas_rate, water_rate};
    }
}

template<class Scalar>
const typename GroupState<Scalar>::GroupPotential&
GroupState<Scalar>::
get_production_group_potential(const std::string& gname) const
{
    return this->production_group_potentials.at(gname);
}

template class GroupState<double>;

#if FLOW_INSTANTIATE_FLOAT
template class GroupState<float>;
#endif

}
