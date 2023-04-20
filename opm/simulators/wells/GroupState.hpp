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

#ifndef OPM_GROUPSTATE_HEADER_INCLUDED
#define OPM_GROUPSTATE_HEADER_INCLUDED

#include <map>
#include <vector>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/simulators/wells/WellContainer.hpp>

namespace Opm {

class GroupState {
public:
    GroupState() = default;
    explicit GroupState(std::size_t num_phases);

    static GroupState serializationTestObject();

    bool operator==(const GroupState& other) const;

    bool has_production_rates(const std::string& gname) const;
    void update_production_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& production_rates(const std::string& gname) const;

    bool has_production_reduction_rates(const std::string& gname) const;
    void update_production_reduction_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& production_reduction_rates(const std::string& gname) const;

    bool has_injection_reduction_rates(const std::string& gname) const;
    void update_injection_reduction_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& injection_reduction_rates(const std::string& gname) const;

    bool has_injection_reservoir_rates(const std::string& gname) const;
    void update_injection_reservoir_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& injection_reservoir_rates(const std::string& gname) const;

    bool has_injection_surface_rates(const std::string& gname) const;
    void update_injection_surface_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& injection_surface_rates(const std::string& gname) const;


    void update_injection_rein_rates(const std::string& gname, const std::vector<double>& rates);
    const std::vector<double>& injection_rein_rates(const std::string& gname) const;

    void update_injection_vrep_rate(const std::string& gname, double rate);
    double injection_vrep_rate(const std::string& gname) const;

    void update_grat_sales_target(const std::string& gname, double target);
    double grat_sales_target(const std::string& gname) const;
    bool has_grat_sales_target(const std::string& gname) const;

    void update_gpmaint_target(const std::string& gname, double target);
    double gpmaint_target(const std::string& gname) const;
    bool has_gpmaint_target(const std::string& gname) const;

    bool has_production_control(const std::string& gname) const;
    void production_control(const std::string& gname, Group::ProductionCMode cmode);
    Group::ProductionCMode production_control(const std::string& gname) const;

    bool has_injection_control(const std::string& gname, Phase phase) const;
    void injection_control(const std::string& gname, Phase phase, Group::InjectionCMode cmode);
    Group::InjectionCMode injection_control(const std::string& gname, Phase phase) const;

    std::size_t data_size() const;
    std::size_t collect(double * data) const;
    std::size_t distribute(const double * data);

    GPMaint::State& gpmaint(const std::string& gname);


    template<class Comm>
    void communicate_rates(const Comm& comm)
    {
        // Note that injection_group_vrep_rates is handled separate from
        // the forAllGroupData() function, since it contains single doubles,
        // not vectors.

        // Create a function that calls some function
        // for all the individual data items to simplify
        // the further code.
        auto iterateContainer = [](auto& container, auto& func) {
            for (auto& x : container) {
                func(x.second);
            }
        };


        auto forAllGroupData = [&](auto& func) {
            iterateContainer(m_production_rates, func);
            iterateContainer(prod_red_rates, func);
            iterateContainer(inj_red_rates, func);
            iterateContainer(inj_resv_rates, func);
            iterateContainer(inj_rein_rates, func);
            iterateContainer(inj_surface_rates, func);
        };

        // Compute the size of the data.
        std::size_t sz = 0;
        auto computeSize = [&sz](const auto& v) {
            sz += v.size();
        };
        forAllGroupData(computeSize);
        sz += this->inj_vrep_rate.size();

        // Make a vector and collect all data into it.
        std::vector<double> data(sz);
        std::size_t pos = 0;



        // That the collect function mutates the vector v is an artifact for
        // testing.
        auto collect = [&data, &pos](auto& v) {
            for (auto& x : v) {
                data[pos++] = x;
                x = -1;
            }
        };
        forAllGroupData(collect);
        for (const auto& x : this->inj_vrep_rate) {
            data[pos++] = x.second;
        }
        if (pos != sz)
            throw std::logic_error("Internal size mismatch when collecting groupData");

        // Communicate it with a single sum() call.
        comm.sum(data.data(), data.size());

        // Distribute the summed vector to the data items.
        pos = 0;
        auto distribute = [&data, &pos](auto& v) {
            for (auto& x : v) {
                x = data[pos++];
            }
        };
        forAllGroupData(distribute);
        for (auto& x : this->inj_vrep_rate) {
            x.second = data[pos++];
        }
        if (pos != sz)
            throw std::logic_error("Internal size mismatch when distributing groupData");
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(num_phases);
        serializer(m_production_rates);
        serializer(production_controls);
        serializer(prod_red_rates);
        serializer(inj_red_rates);
        serializer(inj_surface_rates);
        serializer(inj_resv_rates);
        serializer(inj_rein_rates);
        serializer(inj_vrep_rate);
        serializer(m_grat_sales_target);
        serializer(m_gpmaint_target);
        serializer(injection_controls);
        serializer(gpmaint_state);
    }

private:
    std::size_t num_phases{};
    std::map<std::string, std::vector<double>> m_production_rates;
    std::map<std::string, Group::ProductionCMode> production_controls;
    std::map<std::string, std::vector<double>> prod_red_rates;
    std::map<std::string, std::vector<double>> inj_red_rates;
    std::map<std::string, std::vector<double>> inj_surface_rates;
    std::map<std::string, std::vector<double>> inj_resv_rates;
    std::map<std::string, std::vector<double>> inj_rein_rates;
    std::map<std::string, double> inj_vrep_rate;
    std::map<std::string, double> m_grat_sales_target;
    std::map<std::string, double> m_gpmaint_target;


    std::map<std::pair<Phase, std::string>, Group::InjectionCMode> injection_controls;
    WellContainer<GPMaint::State> gpmaint_state;
};

}

#endif
