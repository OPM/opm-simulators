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

#include <iterator>

#include <opm/json/JsonObject.hpp>

#include <opm/simulators/wells/GroupState.hpp>

namespace Opm {

GroupState::GroupState(std::size_t np) :
    num_phases(np)
{}

bool GroupState::operator==(const GroupState& other) const {
    return this->m_production_rates == other.m_production_rates &&
           this->production_controls == other.production_controls &&
           this->prod_red_rates == other.prod_red_rates &&
           this->inj_red_rates == other.inj_red_rates &&
           this->inj_resv_rates == other.inj_resv_rates &&
           this->inj_potentials == other.inj_potentials &&
           this->inj_rein_rates == other.inj_rein_rates &&
           this->inj_vrep_rate == other.inj_vrep_rate &&
           this->m_grat_sales_target == other.m_grat_sales_target &&
           this->injection_controls == other.injection_controls;
}

//-------------------------------------------------------------------------

bool GroupState::has_production_rates(const std::string& gname) const {
    auto group_iter = this->m_production_rates.find(gname);
    return (group_iter != this->m_production_rates.end());
}

void GroupState::update_production_rates(const std::string& gname, const std::vector<double>& rates) {
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->m_production_rates[gname] = rates;
}

const std::vector<double>& GroupState::production_rates(const std::string& gname) const {
    auto group_iter = this->m_production_rates.find(gname);
    if (group_iter == this->m_production_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

bool GroupState::has_production_reduction_rates(const std::string& gname) const {
    auto group_iter = this->prod_red_rates.find(gname);
    return (group_iter != this->prod_red_rates.end());
}

void GroupState::update_production_reduction_rates(const std::string& gname, const std::vector<double>& rates) {
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->prod_red_rates[gname] = rates;
}

const std::vector<double>& GroupState::production_reduction_rates(const std::string& gname) const {
    auto group_iter = this->prod_red_rates.find(gname);
    if (group_iter == this->prod_red_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

bool GroupState::has_injection_reduction_rates(const std::string& gname) const {
    auto group_iter = this->inj_red_rates.find(gname);
    return (group_iter != this->inj_red_rates.end());
}

void GroupState::update_injection_reduction_rates(const std::string& gname, const std::vector<double>& rates) {
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_red_rates[gname] = rates;
}

const std::vector<double>& GroupState::injection_reduction_rates(const std::string& gname) const {
    auto group_iter = this->inj_red_rates.find(gname);
    if (group_iter == this->inj_red_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

bool GroupState::has_injection_reservoir_rates(const std::string& gname) const {
    auto group_iter = this->inj_resv_rates.find(gname);
    return (group_iter != this->inj_resv_rates.end());
}

void GroupState::update_injection_reservoir_rates(const std::string& gname, const std::vector<double>& rates) {
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_resv_rates[gname] = rates;
}

const std::vector<double>& GroupState::injection_reservoir_rates(const std::string& gname) const {
    auto group_iter = this->inj_resv_rates.find(gname);
    if (group_iter == this->inj_resv_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

void GroupState::update_injection_rein_rates(const std::string& gname, const std::vector<double>& rates) {
    if (rates.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_rein_rates[gname] = rates;
}

const std::vector<double>& GroupState::injection_rein_rates(const std::string& gname) const {
    auto group_iter = this->inj_rein_rates.find(gname);
    if (group_iter == this->inj_rein_rates.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

void GroupState::update_injection_vrep_rate(const std::string& gname, double rate) {
    this->inj_vrep_rate[gname] = rate;
}

double GroupState::injection_vrep_rate(const std::string& gname) const {
    auto group_iter = this->inj_vrep_rate.find(gname);
    if (group_iter == this->inj_vrep_rate.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

void GroupState::update_grat_sales_target(const std::string& gname, double target) {
    this->m_grat_sales_target[gname] = target;
}

double GroupState::grat_sales_target(const std::string& gname) const {
    auto group_iter = this->m_grat_sales_target.find(gname);
    if (group_iter == this->m_grat_sales_target.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

bool GroupState::has_grat_sales_target(const std::string& gname) const {
    return (this->m_grat_sales_target.count(gname) > 0);
}

//-------------------------------------------------------------------------

void GroupState::update_injection_potentials(const std::string& gname, const std::vector<double>& potentials) {
    if (potentials.size() != this->num_phases)
        throw std::logic_error("Wrong number of phases");

    this->inj_potentials[gname] = potentials;
}

const std::vector<double>& GroupState::injection_potentials(const std::string& gname) const {
    auto group_iter = this->inj_potentials.find(gname);
    if (group_iter == this->inj_potentials.end())
        throw std::logic_error("No such group");

    return group_iter->second;
}

//-------------------------------------------------------------------------

bool GroupState::has_production_control(const std::string& gname) const {
    auto group_iter = this->production_controls.find(gname);
    if (group_iter == this->production_controls.end())
        return false;

    return true;
}

void GroupState::production_control(const std::string& gname, Group::ProductionCMode cmode) {
    this->production_controls[gname] = cmode;
}

Group::ProductionCMode GroupState::production_control(const std::string& gname) const {
    auto group_iter = this->production_controls.find(gname);
    if (group_iter == this->production_controls.end())
        throw std::logic_error("Could not find any control for production group: " + gname);

    return group_iter->second;
}

//-------------------------------------------------------------------------

bool GroupState::has_injection_control(const std::string& gname, Phase phase) const {
    return this->injection_controls.count(std::make_pair(phase, gname)) > 0;
}

void GroupState::injection_control(const std::string& gname, Phase phase, Group::InjectionCMode cmode) {
    this->injection_controls[ std::make_pair(phase, gname) ] = cmode;
}

Group::InjectionCMode GroupState::injection_control(const std::string& gname, Phase phase) const {
    auto key = std::make_pair(phase, gname);
    auto group_iter = this->injection_controls.find( key );
    if (group_iter == this->injection_controls.end())
        throw std::logic_error("Could not find ontrol for injection group: " + gname);

    return group_iter->second;
}

//-------------------------------------------------------------------------

GPMaint::State& GroupState::gpmaint(const std::string& gname) {
    if (!this->gpmaint_state.has(gname))
        this->gpmaint_state.add(gname, GPMaint::State{});
    return this->gpmaint_state[gname];
}


//-------------------------------------------------------------------------

namespace {

template <typename T>
void dump_vector(Json::JsonObject& root, const std::string& key, const std::vector<T>& data) {
    auto data_obj = root.add_array(key);
    for (const auto& rate : data)
        data_obj.add(rate);
}


template <typename T>
void dump_groupmap(Json::JsonObject& root, const std::string& key, const std::map<std::string, std::vector<T>>& group_data) {
    auto map_obj = root.add_object(key);
    for (const auto& [gname, data] : group_data)
        dump_vector(map_obj, gname, data);
}


template <typename T>
void dump_groupmap(Json::JsonObject& root, const std::string& key, const std::map<std::string, T>& group_data) {
    auto map_obj = root.add_object(key);
    for (const auto& [gname, value] : group_data)
        map_obj.add_item(gname, value);
}

}

std::string GroupState::dump() const
{
    Json::JsonObject root;
    dump_groupmap(root, "production_rates", this->m_production_rates);
    dump_groupmap(root, "prod_red_rates", this->prod_red_rates);
    dump_groupmap(root, "inj_red_rates", this->inj_red_rates);
    dump_groupmap(root, "inj_resv_rates", this->inj_resv_rates);
    dump_groupmap(root, "inj_potentials", this->inj_potentials);
    dump_groupmap(root, "inj_rein_rates", this->inj_rein_rates);
    dump_groupmap(root, "vrep_rate", this->inj_vrep_rate);
    dump_groupmap(root, "grat_sales_target", this->m_grat_sales_target);
    {
        std::map<std::string, int> int_controls;
        for (const auto& [gname, control] : this->production_controls)
            int_controls.insert({gname, static_cast<int>(control)});
        dump_groupmap(root, "production_controls", int_controls);
    }
    {
        std::map<std::string, std::vector<std::pair<Opm::Phase, Group::InjectionCMode>>> inj_cntrl;
        for (const auto& [phase_name, cmode] : this->injection_controls) {
            const auto& [phase, gname] = phase_name;
            inj_cntrl[gname].emplace_back(phase, cmode);
        }

        auto map_obj = root.add_object("injection_controls");
        for (const auto& [gname, phase_cmode] : inj_cntrl) {
            auto group_array = map_obj.add_array(gname);
            for (const auto& [phase, cmode] : phase_cmode) {
                auto control_pair = group_array.add_array();
                control_pair.add(static_cast<int>(phase));
                control_pair.add(static_cast<int>(cmode));
            }
        }
    }
    return root.to_string();
}


}
