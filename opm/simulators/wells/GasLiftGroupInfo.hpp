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

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <fmt/format.h>

namespace Opm
{
class GasLiftGroupInfo
{
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
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>;
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif

    // TODO: same definition with WellInterface, and
    //   WellState eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;
public:
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
        const Communication& comm);
    std::vector<std::pair<std::string,double>>& getWellGroups(
        const std::string& well_name);

    double alqRate(const std::string& group_name);
    double gasRate(const std::string& group_name);
    int getGroupIdx(const std::string& group_name);
    std::tuple<double,double,double> getRates(int group_idx);
    std::optional<double> gasTarget(const std::string& group_name);
    const std::string& groupIdxToName(int group_idx);
    bool hasWell(const std::string& well_name);
    void initialize();
    std::optional<double> maxAlq(const std::string& group_name);
    double oilRate(const std::string& group_name);
    std::optional<double> oilTarget(const std::string& group_name);
    void update(const std::string& well_name,
        double delta_oil, double delta_gas, double delta_alq);
    void updateRate(int idx, double oil_rate, double gas_rate, double alq);
    const Well2GroupMap& wellGroupMap() { return well_group_map_; }
private:
    bool checkDoGasLiftOptimization_(const std::string& well_name);
    bool checkNewtonIterationIdxOk_(const std::string& well_name);
    void displayDebugMessage_(const std::string& msg);
    void displayDebugMessage_(const std::string& msg, const std::string& well_name);
    std::pair<double, double> getProducerWellRates_(const int index);
    std::tuple<double, double, double>
        initializeGroupRatesRecursive_(const Group &group);
    void initializeWell2GroupMapRecursive_(
        const Group& group, std::vector<std::string>& group_names,
        std::vector<double>& group_efficiency, double cur_efficiency);
    void updateGroupIdxMap_(const std::string& group_name);


    class GroupRates {
    public:
        GroupRates( double oil_rate, double gas_rate, double alq,
            std::optional<double> oil_target,
            std::optional<double> gas_target,
            std::optional<double> total_gas,
            std::optional<double> max_alq
                  ) :
            oil_rate_{oil_rate},
            gas_rate_{gas_rate},
            alq_{alq},
            oil_target_{oil_target},
            gas_target_{gas_target},
            total_gas_{total_gas},
            max_alq_{max_alq}
        {}
        double alq() const { return alq_; }
        void assign(double oil_rate, double gas_rate, double alq)
        {
            oil_rate_ = oil_rate;
            gas_rate_ = gas_rate;
            alq_ = alq;
        }
        double gasRate() const { return gas_rate_; }
        std::optional<double> gasTarget() const { return gas_target_; }
        std::optional<double> maxAlq() const { return max_alq_; }
        double oilRate() const { return oil_rate_; }
        std::optional<double> oilTarget() const { return oil_target_; }
        void update(double delta_oil, double delta_gas, double delta_alq)
        {
            oil_rate_ += delta_oil;
            gas_rate_ += delta_gas;
            alq_ += delta_alq;
        }
    private:
        double oil_rate_;
        double gas_rate_;
        double alq_;
        std::optional<double> oil_target_;
        std::optional<double> gas_target_;
        std::optional<double> total_gas_;
        std::optional<double> max_alq_;
    };

    GLiftEclWells &ecl_wells_;
    const Schedule &schedule_;
    const SummaryState &summary_state_;
    const int report_step_idx_;
    const int iteration_idx_;
    const PhaseUsage &phase_usage_;
    DeferredLogger &deferred_logger_;
    WellState &well_state_;
    const Communication &comm_;
    const GasLiftOpt& glo_;
    GroupRateMap group_rate_map_;
    Well2GroupMap well_group_map_;
    bool debug;
    GroupIdxMap group_idx_;
    int next_group_idx_ = 0;
    // Optimize only wells under THP control
    bool optimize_only_thp_wells_ = false;

};

} // namespace Opm

#endif // OPM_GASLIFT_GROUP_INFO_INCLUDED
