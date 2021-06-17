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
        WellState& well_state);
    std::vector<std::pair<std::string,double>>& getWellGroups(
        const std::string& well_name);

    // TODO: See comment below for initializeGroupRatesRecursive_() for why
    //   the implementation of initialize() is kept here in the header file instead
    //   of in the .cpp file...
    template<class Comm>
    void
    initialize(const Comm& comm)
    {
        const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);
        initializeGroupRatesRecursive_(comm, group);
        std::vector<std::string> group_names;
        std::vector<double> group_efficiency;
        initializeWell2GroupMapRecursive_(
            group, group_names, group_efficiency, /*current efficiency=*/1.0);
    }

    double alqRate(const std::string& group_name);
    double gasRate(const std::string& group_name);
    int getGroupIdx(const std::string& group_name);
    std::tuple<double,double,double> getRates(int group_idx);
    std::optional<double> gasTarget(const std::string& group_name);
    const std::string& groupIdxToName(int group_idx);
    bool hasWell(const std::string& well_name);
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
    void initializeWell2GroupMapRecursive_(
        const Group& group, std::vector<std::string>& group_names,
        std::vector<double>& group_efficiency, double cur_efficiency);
    void updateGroupIdxMap_(const std::string& group_name);

    // TODO: I first tried to pass the MPI Communicator as a constructor argument
    //   to this class (GasLiftGroupInfo) similar to what is done for
    //   GasLiftStage2 (see GasLiftStage2.hpp), hower this did not work for this
    //   class since we are also constructing a GasLiftGroupInfo object in the
    //   test file test1_glift.cpp and when the linker tries to find a definition
    //   of the GasLiftGroupInfo(...) constructor in libopmsimulators.a,
    //   the template type of the MPI communicator (Dune::Communication<..>)
    //   is not of the same type as the one needed by the test case.
    //   The test case needs Dune::Communication<ompi_communicator_t*>, whereas
    //   the one in libopmsimulators.a is Dune::Communication<Dune::No_Comm>.
    //
    //   The error I got from the linker is:
    //
    //   /bin/ld: CMakeFiles/test_glift1.dir/tests/test_glift1.cpp.o:
    //      in function `G1::test_method()':
    //        test_glift1.cpp:(.text+0x15b36): undefined reference to
    //          `Opm::GasLiftGroupInfo::GasLiftGroupInfo(....)
    //
    //   to work around this issue this function is templetized in terms of Comm
    //   here in the header file instead of having it in the .cpp file.
    //   (thanks to Tor Harald S. for the suggestion)
    template<class Comm>
    std::tuple<double, double, double>
    initializeGroupRatesRecursive_(const Comm& comm, const Group &group)
    {
        double oil_rate = 0.0;
        double gas_rate = 0.0;
        double alq = 0.0;
        if (group.wellgroup()) {
            for (const std::string& well_name : group.wells()) {
                // NOTE: we cannot simply use:
                //
                //  const auto &well =
                //    this->schedule_.getWell(well_name, this->report_step_idx_);
                //
                // since the well may not be active (present in the well container)
                auto itr = this->ecl_wells_.find(well_name);
                if (itr != this->ecl_wells_.end()) {
                    const Well *well = (itr->second).first;
                    assert(well); // Should never be nullptr
                    const int index = (itr->second).second;
                    if (well->isProducer()) {
                        auto [sw_oil_rate, sw_gas_rate] = getProducerWellRates_(index);
                        auto sw_alq = this->well_state_.getALQ(well_name);
                        double factor = well->getEfficiencyFactor();
                        oil_rate += (factor * sw_oil_rate);
                        gas_rate += (factor * sw_gas_rate);
                        alq += (factor * sw_alq);
                    }
                }
            }
            // These sums needs to be communictated
            oil_rate = comm.sum(oil_rate);
            gas_rate = comm.sum(gas_rate);
            alq = comm.sum(alq);
        }
        else {
            for (const std::string& group_name : group.groups()) {
                if (!this->schedule_.back().groups.has(group_name))
                    continue;
                const Group& sub_group = this->schedule_.getGroup(
                    group_name, this->report_step_idx_);
                auto [sg_oil_rate, sg_gas_rate, sg_alq]
                    = initializeGroupRatesRecursive_(comm, sub_group);
                const auto gefac = sub_group.getGroupEfficiencyFactor();
                oil_rate += (gefac * sg_oil_rate);
                gas_rate += (gefac * sg_gas_rate);
                alq += (gefac * sg_alq);
            }
        }
        std::optional<double> oil_target, gas_target, max_total_gas, max_alq;
        const auto controls = group.productionControls(this->summary_state_);
        if (group.has_control(Group::ProductionCMode::ORAT)) {
            oil_target = controls.oil_target;
        }
        if (group.has_control(Group::ProductionCMode::GRAT)) {
            gas_target = controls.gas_target;
        }
        if (this->glo_.has_group(group.name())) {
            const auto &gl_group = this->glo_.group(group.name());
            max_alq = gl_group.max_lift_gas();
            max_total_gas = gl_group.max_total_gas();
        }
        if (oil_target || gas_target || max_total_gas || max_alq) {
            updateGroupIdxMap_(group.name());
            this->group_rate_map_.try_emplace(group.name(),
                oil_rate, gas_rate, alq, oil_target, gas_target, max_total_gas, max_alq);
        }
        return std::make_tuple(oil_rate, gas_rate, alq);
    }

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
    const GasLiftOpt& glo_;
    GroupRateMap group_rate_map_;
    Well2GroupMap well_group_map_;
    bool debug;
    GroupIdxMap group_idx_;
    int next_group_idx_ = 0;
    // Optimize only wells under THP control
    bool optimize_only_thp_wells_ = true;

};

} // namespace Opm

#endif // OPM_GASLIFT_GROUP_INFO_INCLUDED
