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

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <opm/simulators/wells/WellContainer.hpp>

#include <map>
#include <vector>
#include <utility>

namespace Opm {

    class GConSump;
    class Schedule;
    class SummaryState;

template<class Scalar>
class GroupState {
public:
    /// Flags indicating which rate vectors have been modified since last communication.
    /// Used by communicate_rates() to selectively communicate only modified data.
    enum class Modification : uint32_t {
        NONE                    = 0,
        PRODUCTION_RATES        = 1 << 0,  // m_production_rates
        NETWORK_RATES           = 1 << 1,  // m_network_leaf_node_production_rates
        PROD_REDUCTION_RATES    = 1 << 2,  // prod_red_rates
        INJ_REDUCTION_RATES     = 1 << 3,  // inj_red_rates
        INJ_RESV_RATES          = 1 << 4,  // inj_resv_rates
        INJ_REIN_RATES          = 1 << 5,  // inj_rein_rates
        INJ_SURFACE_RATES       = 1 << 6,  // inj_surface_rates
        INJ_VREP_RATE           = 1 << 7,  // inj_vrep_rate
        ALL                     = 0xFF
    };

    GroupState() = default;
    explicit GroupState(std::size_t num_phases);

    static GroupState serializationTestObject();

    bool operator==(const GroupState& other) const;

    bool has_production_rates(const std::string& gname) const;
    void update_production_rates(const std::string& gname,
                                 const std::vector<Scalar>& rates);
    void update_network_leaf_node_production_rates(const std::string& gname,
                                 const std::vector<Scalar>& rates);
    const std::vector<Scalar>& production_rates(const std::string& gname) const;
    const std::vector<Scalar>& network_leaf_node_production_rates(const std::string& gname) const;

    void update_well_group_thp(const std::string& gname, const double& thp);
    Scalar well_group_thp(const std::string& gname) const;
    bool is_autochoke_group(const std::string& gname) const;

    bool has_production_reduction_rates(const std::string& gname) const;
    void update_production_reduction_rates(const std::string& gname,
                                           const std::vector<Scalar>& rates);
    const std::vector<Scalar>& production_reduction_rates(const std::string& gname) const;

    bool has_injection_reduction_rates(const std::string& gname) const;
    void update_injection_reduction_rates(const std::string& gname,
                                          const std::vector<Scalar>& rates);
    const std::vector<Scalar>& injection_reduction_rates(const std::string& gname) const;

    bool has_injection_reservoir_rates(const std::string& gname) const;
    void update_injection_reservoir_rates(const std::string& gname,
                                          const std::vector<Scalar>& rates);
    const std::vector<Scalar>& injection_reservoir_rates(const std::string& gname) const;

    bool has_injection_surface_rates(const std::string& gname) const;
    void update_injection_surface_rates(const std::string& gname,
                                        const std::vector<Scalar>& rates);
    const std::vector<Scalar>& injection_surface_rates(const std::string& gname) const;

    void update_injection_rein_rates(const std::string& gname,
                                     const std::vector<Scalar>& rates);
    const std::vector<Scalar>& injection_rein_rates(const std::string& gname) const;

    void update_injection_vrep_rate(const std::string& gname, Scalar rate);
    Scalar injection_vrep_rate(const std::string& gname) const;

    void update_grat_sales_target(const std::string& gname, Scalar target);
    Scalar grat_sales_target(const std::string& gname) const;
    bool has_grat_sales_target(const std::string& gname) const;

    void update_gpmaint_target(const std::string& gname, Scalar target);
    Scalar gpmaint_target(const std::string& gname) const;
    bool has_gpmaint_target(const std::string& gname) const;

    bool has_field_or_none_control(const std::string& gname) const;
    bool has_field_or_none_control(const std::string& gname, Phase injection_phase) const;
    bool has_production_control(const std::string& gname) const;
    void production_control(const std::string& gname, Group::ProductionCMode cmode);
    Group::ProductionCMode production_control(const std::string& gname) const;
    const std::map<std::string, Group::ProductionCMode>& get_production_controls() const;

    bool has_injection_control(const std::string& gname, Phase phase) const;
    void injection_control(const std::string& gname, Phase phase, Group::InjectionCMode cmode);
    Group::InjectionCMode injection_control(const std::string& gname, Phase phase) const;

    void update_number_of_wells_under_group_control(const std::string& gname, int number);
    int number_of_wells_under_group_control(const std::string& gname) const;

    void update_number_of_wells_under_inj_group_control(const std::string& gname, Phase phase, int number);
    int number_of_wells_under_inj_group_control(const std::string& gname, Phase phase) const;


    void update_gconsump(const Schedule& schedule, const int report_step, const SummaryState& summary_state);
    const std::pair<Scalar, Scalar>& gconsump_rates(const std::string& gname) const;

    struct GroupPotential {
        Scalar oil_rate;
        Scalar gas_rate;
        Scalar water_rate;

        GroupPotential(Scalar oil = 0.0, Scalar gas = 0.0, Scalar water = 0.0)
            : oil_rate(oil), gas_rate(gas), water_rate(water) {}
    };
    bool has_production_group_potential(const std::string& gname) const;
    void update_group_production_potential(
        const std::string& gname, Scalar oil_rate, Scalar gas_rate, Scalar water_rate
    );
    const GroupPotential& get_production_group_potential(const std::string& gname) const;

    std::size_t data_size() const;
    std::size_t collect(Scalar* data) const;
    std::size_t distribute(const Scalar* data);

    GPMaint::State& gpmaint(const std::string& gname);

    /// Check if a modification flag is set (for testing/debugging)
    bool is_modified(Modification mod) const {
        return (modifications_ & static_cast<uint32_t>(mod)) != 0;
    }

    template<class Comm>
    void communicate_rates(const Comm& comm)
    {
        // Note that injection_group_vrep_rates is handled separate from
        // the forModifiedData() function, since it contains single doubles,
        // not vectors.

        // Create a function that calls some function
        // for all the individual data items to simplify
        // the further code.
        auto iterateContainer = [](auto& container, const auto& func) {
            for (auto& x : container) {
                func(x.second);
            }
        };

        // Only include modified vectors in communication
        auto forModifiedData = [&](auto& func) {
            if (is_modified(Modification::PRODUCTION_RATES))
                iterateContainer(m_production_rates, func);
            if (is_modified(Modification::NETWORK_RATES))
                iterateContainer(m_network_leaf_node_production_rates, func);
            if (is_modified(Modification::PROD_REDUCTION_RATES))
                iterateContainer(prod_red_rates, func);
            if (is_modified(Modification::INJ_REDUCTION_RATES))
                iterateContainer(inj_red_rates, func);
            if (is_modified(Modification::INJ_RESV_RATES))
                iterateContainer(inj_resv_rates, func);
            if (is_modified(Modification::INJ_REIN_RATES))
                iterateContainer(inj_rein_rates, func);
            if (is_modified(Modification::INJ_SURFACE_RATES))
                iterateContainer(inj_surface_rates, func);
        };

        // Compute the size of the modified data.
        std::size_t sz = 0;
        auto computeSize = [&sz](const auto& v) {
            sz += v.size();
        };
        forModifiedData(computeSize);
        if (is_modified(Modification::INJ_VREP_RATE))
            sz += this->inj_vrep_rate.size();

        // Early exit if nothing to communicate
        if (sz == 0) {
            return;
        }

        // Make a vector and collect all data into it.
        std::vector<Scalar> data(sz);
        std::size_t pos = 0;

        // That the collect function mutates the vector v is an artifact for
        // testing.
        auto doCollect = [&data, &pos](auto& v) {
            for (auto& x : v) {
                data[pos++] = x;
                x = -1;
            }
        };
        forModifiedData(doCollect);
        if (is_modified(Modification::INJ_VREP_RATE)) {
            for (auto& x : this->inj_vrep_rate) {
                data[pos++] = x.second;
                x.second = -1;
            }
        }
        if (pos != sz)
            throw std::logic_error("Internal size mismatch when collecting groupData");

        // Communicate it with a single sum() call.
        comm.sum(data.data(), data.size());

        // Distribute the summed vector to the data items.
        pos = 0;
        auto doDistribute = [&data, &pos](auto& v) {
            std::copy_n(data.begin() + pos, v.size(), v.begin());
            pos += v.size();
        };
        forModifiedData(doDistribute);
        if (is_modified(Modification::INJ_VREP_RATE)) {
            for (auto& x : this->inj_vrep_rate) {
                x.second = data[pos++];
            }
        }
        if (pos != sz)
            throw std::logic_error("Internal size mismatch when distributing groupData");

        // Clear modification flags after successful communication
        clear_modifications();
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(num_phases);
        serializer(m_production_rates);
        serializer(m_network_leaf_node_production_rates);
        serializer(production_controls);
        serializer(group_thp);
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
        serializer(m_gconsump_rates);
        serializer(m_number_of_wells_under_group_control);
        serializer(m_number_of_wells_under_inj_group_control);
    }

private:
    std::size_t num_phases{};
    std::map<std::string, std::vector<Scalar>> m_production_rates;
    std::map<std::string, std::vector<Scalar>> m_network_leaf_node_production_rates;
    std::map<std::string, Group::ProductionCMode> production_controls;
    std::map<std::string, std::vector<Scalar>> prod_red_rates;
    std::map<std::string, std::vector<Scalar>> inj_red_rates;
    std::map<std::string, std::vector<Scalar>> inj_surface_rates;
    std::map<std::string, std::vector<Scalar>> inj_resv_rates;
    std::map<std::string, std::vector<Scalar>> inj_rein_rates;
    std::map<std::string, Scalar> inj_vrep_rate;
    std::map<std::string, Scalar> m_grat_sales_target;
    std::map<std::string, Scalar> m_gpmaint_target;
    std::map<std::string, Scalar> group_thp;
    std::map<std::string, GroupPotential> production_group_potentials;
    std::map<std::string, int> m_number_of_wells_under_group_control;
    std::map<std::pair<Phase, std::string>, int> m_number_of_wells_under_inj_group_control;


    std::map<std::pair<Phase, std::string>, Group::InjectionCMode> injection_controls;
    WellContainer<GPMaint::State> gpmaint_state;
    std::map<std::string, std::pair<Scalar, Scalar>> m_gconsump_rates; // Pair with {consumption_rate, import_rate} for each group
    static constexpr std::pair<Scalar, Scalar> zero_pair = {0.0, 0.0};

    /// Bitmask tracking which rate vectors have been modified since last communication
    uint32_t modifications_ = 0;

    /// Mark a vector type as modified (called by update_*() methods)
    void mark_modified(Modification mod) {
        modifications_ |= static_cast<uint32_t>(mod);
    }

    /// Clear all modification flags (called after successful communication)
    void clear_modifications() {
        modifications_ = 0;
    }
};

}

#endif
