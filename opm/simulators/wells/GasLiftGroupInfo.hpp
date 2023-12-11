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

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/GasLiftCommon.hpp>

#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GasLiftOpt;
class Group;
class GroupState;
class Schedule;
class SummaryState;
class Well;
class WellState;

class GasLiftGroupInfo : public GasLiftCommon
{
protected:
    class GroupRates;
    // NOTE: In the Well2GroupMap below, in the std::vector value we store
    //    pairs of group names and efficiency factors. The efficiency factors
    //    are the product of the wells efficiency factor and all the efficiency
    //    factors of the child groups of the group all the way down
    //    to the well group.
    using Well2GroupMap =
        std::map<std::string, std::vector<std::pair<std::string,double>>>;
    using GroupRateMap =
        std::map<std::string, GroupRates>;
    using GroupIdxMap = std::map<std::string, int>;
    using Communication = Dune::Communication<Dune::MPIHelper::MPICommunicator>;

    // TODO: same definition with WellInterface, and
    //   WellState eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;
public:
    enum class Rate {oil, gas, water, liquid};

    using GLiftEclWells = std::map<std::string,std::pair<const Well *,int>>;
    GasLiftGroupInfo(
        GLiftEclWells& ecl_wells,
        const Schedule& schedule,
        const SummaryState& summary_state,
        const int report_step_idx,
        const int iteration_idx,
        const PhaseUsage& phase_usage,
        DeferredLogger& deferred_logger,
        WellState& well_state,
        const GroupState& group_state,
        const Parallel::Communication& comm,
        bool glift_debug
    );
    std::vector<std::pair<std::string,double>>& getWellGroups(
        const std::string& well_name);

    double alqRate(const std::string& group_name);
    double gasRate(const std::string& group_name) const;
    double gasPotential(const std::string& group_name) const;
    double waterPotential(const std::string& group_name) const;
    double oilPotential(const std::string& group_name) const;
    int getGroupIdx(const std::string& group_name);
    double getRate(Rate rate_type, const std::string& group_name) const;
    double getPotential(Rate rate_type, const std::string& group_name) const;
    std::tuple<double,double,double,double> getRates(const int group_idx) const;
    std::optional<double> gasTarget(const std::string& group_name) const;
    std::optional<double> getTarget(
        Rate rate_type, const std::string& group_name) const;
    const std::string& groupIdxToName(int group_idx) const;
    bool hasAnyTarget(const std::string& group_name) const;
    bool hasWell(const std::string& well_name);
    void initialize();
    std::optional<double> liquidTarget(const std::string& group_name) const;
    std::optional<double> maxAlq(const std::string& group_name);
    std::optional<double> maxTotalGasRate(const std::string& group_name);
    double oilRate(const std::string& group_name) const;
    std::optional<double> oilTarget(const std::string& group_name) const;
    static const std::string rateToString(Rate rate);
    double waterRate(const std::string& group_name) const;
    std::optional<double> waterTarget(const std::string& group_name) const;
    void update(const std::string& well_name,
        double delta_oil, double delta_gas, double delta_water, double delta_alq);
    void updateRate(int idx, double oil_rate, double gas_rate, double water_rate, double alq);
    const Well2GroupMap& wellGroupMap() { return well_group_map_; }
protected:
    bool checkDoGasLiftOptimization_(const std::string& well_name);
    bool checkNewtonIterationIdxOk_(const std::string& well_name);
    void debugDisplayWellContribution_(
        const std::string& gr_name, const std::string& well_name,
        double eff_factor,
        double well_oil_rate, double well_gas_rate, double well_water_rate,
        double well_alq,
        double oil_rate, double gas_rate, double water_rate,
        double alq
    ) const;
    void debugDisplayUpdatedGroupRates(const std::string& name,
        double oil_rate, double gas_rate, double water_rate, double alq) const;
    void debugEndInitializeGroup(const std::string& name) const;
    void debugStartInitializeGroup(const std::string& name) const;
    void displayDebugMessage_(const std::string& msg) const override;
    void displayDebugMessage_(const std::string& msg, const std::string& well_name);
    std::tuple<double, double, double, double, double, double>
        getProducerWellRates_(const Well* well, const int index);
    std::tuple<double, double, double, double, double, double, double>
        initializeGroupRatesRecursive_(const Group &group);
    void initializeWell2GroupMapRecursive_(
        const Group& group, std::vector<std::string>& group_names,
        std::vector<double>& group_efficiency, double cur_efficiency);
    void updateGroupIdxMap_(const std::string& group_name);


    class GroupRates {
    public:
        GroupRates( double oil_rate, double gas_rate, double water_rate, double alq,
            double oil_potential, double gas_potential, double water_potential,
            std::optional<double> oil_target,
            std::optional<double> gas_target,
            std::optional<double> water_target,
            std::optional<double> liquid_target,
            std::optional<double> total_gas,
            std::optional<double> max_alq
                  ) :
            oil_rate_{oil_rate},
            gas_rate_{gas_rate},
            water_rate_{water_rate},
            alq_{alq},
            oil_potential_{oil_potential},
            gas_potential_{gas_potential},
            water_potential_{water_potential},
            oil_target_{oil_target},
            gas_target_{gas_target},
            water_target_{water_target},
            liquid_target_{liquid_target},
            total_gas_{total_gas},
            max_alq_{max_alq}
        {}
        double alq() const { return alq_; }
        void assign(double oil_rate, double gas_rate, double water_rate, double alq)
        {
            oil_rate_ = oil_rate;
            gas_rate_ = gas_rate;
            water_rate_ = water_rate;
            alq_ = alq;
        }
        double gasRate() const { return gas_rate_; }
        double waterRate() const { return water_rate_; }
        std::optional<double> gasTarget() const { return gas_target_; }
        std::optional<double> waterTarget() const { return water_target_; }
        std::optional<double> maxAlq() const { return max_alq_; }
        std::optional<double> maxTotalGasRate() const { return total_gas_; }
        double oilRate() const { return oil_rate_; }
        std::optional<double> oilTarget() const { return oil_target_; }
        std::optional<double> liquidTarget() const { return liquid_target_; }
        double oilPotential() const { return oil_potential_; }
        double gasPotential() const { return gas_potential_; }
        double waterPotential() const { return water_potential_; }

        void update(double delta_oil, double delta_gas, double delta_water, double delta_alq)
        {
            oil_rate_ += delta_oil;
            gas_rate_ += delta_gas;
            water_rate_ += delta_water;
            alq_ += delta_alq;
            // Note. We don't updata the potentials at this point. They
            // are only needed initially.
        }
    private:
        double oil_rate_;
        double gas_rate_;
        double water_rate_;
        double alq_;
        double oil_potential_;
        double gas_potential_;
        double water_potential_;
        std::optional<double> oil_target_;
        std::optional<double> gas_target_;
        std::optional<double> water_target_;
        std::optional<double> liquid_target_;
        std::optional<double> total_gas_;
        std::optional<double> max_alq_;
    };

    GLiftEclWells &ecl_wells_;
    const Schedule &schedule_;
    const SummaryState &summary_state_;
    const int report_step_idx_;
    const int iteration_idx_;
    const PhaseUsage &phase_usage_;
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
