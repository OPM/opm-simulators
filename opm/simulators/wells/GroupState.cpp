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

    return result;
}

template<class Scalar>
bool GroupState<Scalar>::operator==(const GroupState& other) const
{
    return this->m_production_rates == other.m_production_rates &&
           this->production_controls == other.production_controls &&
           this->prod_red_rates == other.prod_red_rates &&
           this->inj_red_rates == other.inj_red_rates &&
           this->inj_resv_rates == other.inj_resv_rates &&
           this->inj_rein_rates == other.inj_rein_rates &&
           this->inj_vrep_rate == other.inj_vrep_rate &&
           this->inj_surface_rates == other.inj_surface_rates &&
           this->m_grat_sales_target == other.m_grat_sales_target &&
           this->injection_controls == other.injection_controls &&
           this->gpmaint_state == other.gpmaint_state;
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
const std::vector<Scalar>&
GroupState<Scalar>::production_rates(const std::string& gname) const
{
    auto group_iter = this->m_production_rates.find(gname);
    if (group_iter == this->m_production_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
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

template class GroupState<double>;

}
