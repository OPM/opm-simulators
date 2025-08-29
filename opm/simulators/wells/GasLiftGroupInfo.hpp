/*
  Copyright 2021 Equinor ASA.

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

#ifndef OPM_GASLIFT_GROUP_INFO_HEADER_INCLUDED
#define OPM_GASLIFT_GROUP_INFO_HEADER_INCLUDED

#include <opm/simulators/wells/GasLiftCommon.hpp>

#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace Opm {

class DeferredLogger;
class GasLiftOpt;
class Group;
template<class Scalar> class GroupState;
template<typename IndexTraits> class PhaseUsageInfo;
class Schedule;
class SummaryState;
class Well;
template<typename Scalar, typename IndexTraits> class WellState;


template<typename Scalar, typename IndexTraits>
class GasLiftGroupInfo : public GasLiftCommon<Scalar, IndexTraits>
{
protected:
    class GroupRates;
    // NOTE: In the Well2GroupMap below, in the std::vector value we store
    //    pairs of group names and efficiency factors. The efficiency factors
    //    are the product of the wells efficiency factor and all the efficiency
    //    factors of the child groups of the group all the way down
    //    to the well group.
    using Well2GroupMap =
        std::map<std::string, std::vector<std::pair<std::string, Scalar>>>;
    using GroupRateMap =
        std::map<std::string, GroupRates>;
    using GroupIdxMap = std::map<std::string, int>;
    using Communication = Dune::Communication<Dune::MPIHelper::MPICommunicator>;

public:
    enum class Rate {oil, gas, water, liquid};

    using GLiftEclWells = std::map<std::string,std::pair<const Well *,int>>;
    GasLiftGroupInfo(GLiftEclWells& ecl_wells,
                     const Schedule& schedule,
                     const SummaryState& summary_state,
                     const int report_step_idx,
                     const int iteration_idx,
                     DeferredLogger& deferred_logger,
                     WellState<Scalar, IndexTraits>& well_state,
                     const GroupState<Scalar>& group_state,
                     const Parallel::Communication& comm,
                     bool glift_debug);

    std::vector<std::pair<std::string,Scalar>>&
    getWellGroups(const std::string& well_name);

    Scalar alqRate(const std::string& group_name);
    Scalar gasRate(const std::string& group_name) const;
    Scalar gasPotential(const std::string& group_name) const;
    Scalar waterPotential(const std::string& group_name) const;
    Scalar oilPotential(const std::string& group_name) const;
    int getGroupIdx(const std::string& group_name);
    Scalar getRate(Rate rate_type, const std::string& group_name) const;
    Scalar getPotential(Rate rate_type, const std::string& group_name) const;
    std::tuple<Scalar,Scalar,Scalar,Scalar> getRates(const int group_idx) const;
    std::optional<Scalar> gasTarget(const std::string& group_name) const;
    std::optional<Scalar> getTarget(Rate rate_type, const std::string& group_name) const;
    const std::string& groupIdxToName(int group_idx) const;
    bool hasAnyTarget(const std::string& group_name) const;
    bool hasWell(const std::string& well_name);
    void initialize();
    std::optional<Scalar> liquidTarget(const std::string& group_name) const;
    std::optional<Scalar> maxAlq(const std::string& group_name);
    std::optional<Scalar> maxTotalGasRate(const std::string& group_name);
    Scalar oilRate(const std::string& group_name) const;
    std::optional<Scalar> oilTarget(const std::string& group_name) const;
    static const std::string rateToString(Rate rate);
    Scalar waterRate(const std::string& group_name) const;
    std::optional<Scalar> waterTarget(const std::string& group_name) const;
    void update(const std::string& well_name,
                Scalar delta_oil,
                Scalar delta_gas,
                Scalar delta_water,
                Scalar delta_alq);
    void updateRate(int idx,
                    Scalar oil_rate,
                    Scalar gas_rate,
                    Scalar water_rate,
                    Scalar alq);
    const Well2GroupMap& wellGroupMap() { return well_group_map_; }

protected:
    bool checkDoGasLiftOptimization_(const std::string& well_name);
    bool checkNewtonIterationIdxOk_(const std::string& well_name);
    void debugDisplayWellContribution_(const std::string& gr_name,
                                       const std::string& well_name,
                                       Scalar eff_factor,
                                       Scalar well_oil_rate,
                                       Scalar well_gas_rate,
                                       Scalar well_water_rate,
                                       Scalar well_alq,
                                       Scalar oil_rate,
                                       Scalar gas_rate,
                                       Scalar water_rate,
                                       Scalar alq) const;
    void debugDisplayUpdatedGroupRates(const std::string& name,
                                       Scalar oil_rate,
                                       Scalar gas_rate,
                                       Scalar water_rate,
                                       Scalar alq) const;
    void debugEndInitializeGroup(const std::string& name) const;
    void debugStartInitializeGroup(const std::string& name) const;
    void displayDebugMessage_(const std::string& msg) const override;
    void displayDebugMessage_(const std::string& msg, const std::string& well_name);

    std::tuple<Scalar, Scalar, Scalar, Scalar, Scalar, Scalar>
    getProducerWellRates_(const Well* well, const int index);

    std::tuple<Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar>
    initializeGroupRatesRecursive_(const Group &group);

    void initializeWell2GroupMapRecursive_(const Group& group,
                                           std::vector<std::string>& group_names,
                                           std::vector<Scalar>& group_efficiency,
                                           Scalar cur_efficiency);
    void updateGroupIdxMap_(const std::string& group_name);

    class GroupRates {
    public:
        GroupRates(Scalar oil_rate,
                   Scalar gas_rate,
                   Scalar water_rate,
                   Scalar alq,
                   Scalar oil_potential,
                   Scalar gas_potential,
                   Scalar water_potential,
                   std::optional<Scalar> oil_target,
                   std::optional<Scalar> gas_target,
                   std::optional<Scalar> water_target,
                   std::optional<Scalar> liquid_target,
                   std::optional<Scalar> total_gas,
                   std::optional<Scalar> max_alq)
            : oil_rate_{oil_rate}
            , gas_rate_{gas_rate}
            , water_rate_{water_rate}
            , alq_{alq}
            , oil_potential_{oil_potential}
            , gas_potential_{gas_potential}
            , water_potential_{water_potential}
            , oil_target_{oil_target}
            , gas_target_{gas_target}
            , water_target_{water_target}
            , liquid_target_{liquid_target}
            , total_gas_{total_gas}
            , max_alq_{max_alq}
        {}

        Scalar alq() const { return alq_; }
        void assign(Scalar oil_rate,
                    Scalar gas_rate,
                    Scalar water_rate,
                    Scalar alq)
        {
            oil_rate_ = oil_rate;
            gas_rate_ = gas_rate;
            water_rate_ = water_rate;
            alq_ = alq;
        }
        Scalar gasRate() const { return gas_rate_; }
        Scalar waterRate() const { return water_rate_; }
        std::optional<Scalar> gasTarget() const { return gas_target_; }
        std::optional<Scalar> waterTarget() const { return water_target_; }
        std::optional<Scalar> maxAlq() const { return max_alq_; }
        std::optional<Scalar > maxTotalGasRate() const { return total_gas_; }
        Scalar oilRate() const { return oil_rate_; }
        std::optional<Scalar> oilTarget() const { return oil_target_; }
        std::optional<Scalar> liquidTarget() const { return liquid_target_; }
        Scalar oilPotential() const { return oil_potential_; }
        Scalar gasPotential() const { return gas_potential_; }
        Scalar waterPotential() const { return water_potential_; }

        void update(Scalar delta_oil,
                    Scalar delta_gas,
                    Scalar delta_water,
                    Scalar delta_alq)
        {
            oil_rate_ += delta_oil;
            gas_rate_ += delta_gas;
            water_rate_ += delta_water;
            alq_ += delta_alq;
            // Note. We don't updata the potentials at this point. They
            // are only needed initially.
        }

    private:
        Scalar oil_rate_;
        Scalar gas_rate_;
        Scalar water_rate_;
        Scalar alq_;
        Scalar oil_potential_;
        Scalar gas_potential_;
        Scalar water_potential_;
        std::optional<Scalar> oil_target_;
        std::optional<Scalar> gas_target_;
        std::optional<Scalar> water_target_;
        std::optional<Scalar> liquid_target_;
        std::optional<Scalar> total_gas_;
        std::optional<Scalar> max_alq_;
    };

    GLiftEclWells& ecl_wells_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const int report_step_idx_;
    const int iteration_idx_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
    const GasLiftOpt& glo_;
    GroupRateMap group_rate_map_;
    Well2GroupMap well_group_map_;
    GroupIdxMap group_idx_;
    int next_group_idx_ = 0;
    // Optimize only wells under THP control
    bool optimize_only_thp_wells_ = false;
};

} // namespace Opm

#endif // OPM_GASLIFT_GROUP_INFO_INCLUDED
