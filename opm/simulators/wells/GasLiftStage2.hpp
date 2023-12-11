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

#ifndef OPM_GASLIFT_STAGE2_HEADER_INCLUDED
#define OPM_GASLIFT_STAGE2_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GasLiftOpt;
class GasLiftWellState;
class Group;
class GroupState;
class Schedule;
class WellInterfaceGeneric;
class WellState;

class GasLiftStage2 : public GasLiftCommon {
    using GasLiftSingleWell = GasLiftSingleWellGeneric;
    using GLiftOptWells = std::map<std::string,std::unique_ptr<GasLiftSingleWell>>;
    using GLiftProdWells = std::map<std::string,const WellInterfaceGeneric*>;
    using GLiftWellStateMap = std::map<std::string,std::unique_ptr<GasLiftWellState>>;
    using GradPair = std::pair<std::string, double>;
    using GradPairItr = std::vector<GradPair>::iterator;
    using GradInfo = typename GasLiftSingleWellGeneric::GradInfo;
    using GradMap = std::map<std::string, GradInfo>;
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;
public:
    GasLiftStage2(
        const int report_step_idx,
        const Parallel::Communication& comm,
        const Schedule& schedule,
        const SummaryState& summary_state,
        DeferredLogger& deferred_logger,
        WellState& well_state,
        const GroupState& group_state,
        GLiftProdWells& prod_wells,
        GLiftOptWells& glift_wells,
        GasLiftGroupInfo& group_info,
        GLiftWellStateMap& state_map,
        bool glift_debug
    );
    void runOptimize();
protected:
    void addOrRemoveALQincrement_(
        GradMap& grad_map, const std::string& well_name, bool add);
    std::optional<GradInfo> calcIncOrDecGrad_(
        const std::string name, const GasLiftSingleWell& gs_well, const std::string& gr_name_dont_limit, bool increase);
    bool checkRateAlreadyLimited_(const std::string& well_name, GasLiftWellState& state, bool increase);
    GradInfo deleteDecGradItem_(const std::string& name);
    GradInfo deleteIncGradItem_(const std::string& name);
    GradInfo deleteGrad_(const std::string& name, bool increase);
    void displayDebugMessage_(const std::string& msg) const override;
    void displayDebugMessage2B_(const std::string& msg);
    void displayDebugMessage_(const std::string& msg, const std::string& group_name);
    void displayWarning_(const std::string& msg, const std::string& group_name);
    void displayWarning_(const std::string& msg);
    std::tuple<double, double, double, double> getCurrentGroupRates_(const Group& group);
    std::optional<double> getGroupMaxALQ_(const Group &group);
    std::optional<double> getGroupMaxTotalGas_(const Group &group);
    std::vector<GasLiftSingleWell *> getGroupGliftWells_(
        const Group& group);
    void getGroupGliftWellsRecursive_(
        const Group& group, std::vector<GasLiftSingleWell *>& wells);
    void optimizeGroup_(const Group& group);
    void optimizeGroupsRecursive_(const Group& group);
    void recalculateGradientAndUpdateData_(
        GradPairItr& grad_itr, const std::string& gr_name_dont_limit, bool increase,
        std::vector<GradPair>& grads, std::vector<GradPair>& other_grads);
    void redistributeALQ_(
        std::vector<GasLiftSingleWell *>& wells,  const Group& group,
        std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
    void removeSurplusALQ_(
        const Group& group,
        std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
    void saveGrad_(GradMap& map, const std::string& name, GradInfo& grad);
    void saveDecGrad_(const std::string& name, GradInfo& grad);
    void saveIncGrad_(const std::string& name, GradInfo& grad);
    void sortGradients_(std::vector<GradPair>& grads);
    std::optional<GradInfo> updateGrad_(
        const std::string& name, GradInfo& grad, bool increase);
    void updateGradVector_(
        const std::string& name, std::vector<GradPair>& grads, double grad);
    void mpiSyncGlobalGradVector_(std::vector<GradPair>& grads_global) const;
    void mpiSyncLocalToGlobalGradVector_(
        const std::vector<GradPair>& grads_local,
        std::vector<GradPair>& grads_global) const;


    GLiftProdWells& prod_wells_;
    GLiftOptWells& stage1_wells_;
    GasLiftGroupInfo& group_info_;
    GLiftWellStateMap& well_state_map_;

    int report_step_idx_;
    const SummaryState& summary_state_;
    const Schedule& schedule_;
    const GasLiftOpt& glo_;
    GradMap inc_grads_;
    GradMap dec_grads_;
    int max_iterations_ = 1000;
    //int time_step_idx_;

    struct OptimizeState {
        OptimizeState(GasLiftStage2& parent_, const Group& group_ ) :
            parent{parent_},
            group{group_},
            it{0}
        {}
        GasLiftStage2& parent;
        const Group& group;
        int it;

        using GradInfo = typename GasLiftStage2::GradInfo;
        using GradPair = typename GasLiftStage2::GradPair;
        using GradPairItr = typename GasLiftStage2::GradPairItr;
        using GradMap = typename GasLiftStage2::GradMap;
        void calculateEcoGradients(std::vector<GasLiftSingleWell *>& wells,
            std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
        bool checkAtLeastTwoWells(std::vector<GasLiftSingleWell *>& wells);
        void debugShowIterationInfo();
        std::pair<std::optional<GradPairItr>,std::optional<GradPairItr>>
        getEcoGradients(
            std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
        void recalculateGradients(
            std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads,
            GradPairItr& min_dec_grad_itr, GradPairItr &max_inc_grad_itr);
        void redistributeALQ( GradPairItr& min_dec_grad, GradPairItr& max_inc_grad);

    private:
        void displayDebugMessage_(const std::string& msg);
        void displayWarning_(const std::string& msg);

    };

    struct SurplusState {
        SurplusState( GasLiftStage2& parent_, const Group& group_,
              double oil_rate_, double gas_rate_, double water_rate_, double alq_,
              double min_eco_grad_,
              double oil_target_, double gas_target_, double water_target_, double liquid_target_,
              std::optional<double> max_glift_, std::optional<double> max_total_gas_) :
            parent{parent_},
            group{group_},
            oil_rate{oil_rate_},
            gas_rate{gas_rate_},
            water_rate{water_rate_},
            alq{alq_},
            min_eco_grad{min_eco_grad_},
            oil_target{oil_target_},
            gas_target{gas_target_},
            water_target(water_target_),
            liquid_target{liquid_target_},
            max_glift{max_glift_},
            max_total_gas{max_total_gas_},
            it{0}
        {}
        GasLiftStage2 &parent;
        const Group &group;
        double oil_rate;
        double gas_rate;
        double water_rate;
        double alq;
        const double min_eco_grad;
        const double oil_target;
        const double gas_target;
        const double water_target;
        const double liquid_target;
        std::optional<double> max_glift;
        std::optional<double> max_total_gas;
        int it;

        void addOrRemoveALQincrement(
            GradMap &grad_map, const std::string& well_name, bool add);
        bool checkALQlimit(double delta_alq, double delta_gas);
        bool checkEcoGradient(const std::string& well_name, double eco_grad);
        bool checkGasTarget(double delta_gas);
        bool checkLiquidTarget(double delta_liquid);
        bool checkOilTarget(double delta_oil);
        bool checkWaterTarget(double delta_water);
        std::array<double, 4> computeDelta(const std::string& name);
        void updateRates(const std::array<double, 4>& delta);
    };
};

} // namespace Opm

#endif // OPM_GASLIFT_STAGE2_HEADER_INCLUDED
