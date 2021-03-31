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

#include <ebos/eclproblem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
// NOTE: BlackoilWellModel.hpp includes ourself (GasLiftStage2.hpp), so we need
//   to forward declare BlackoilWellModel for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class BlackoilWellModel;
}
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftStage2 {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using BlackoilWellModel = Opm::BlackoilWellModel<TypeTag>;
        using GasLiftSingleWell = Opm::GasLiftSingleWell<TypeTag>;
        using GLiftWellState = Opm::GasLiftWellState<TypeTag>;
        using GLiftOptWells = typename BlackoilWellModel::GLiftOptWells;
        using GLiftProdWells = typename BlackoilWellModel::GLiftProdWells;
        using GLiftWellStateMap = typename BlackoilWellModel::GLiftWellStateMap;
        using GradPair = std::pair<std::string, double>;
        using GradPairItr = std::vector<GradPair>::iterator;
        using GradInfo = typename GasLiftSingleWell::GradInfo;
        using GradMap = std::map<std::string, GradInfo>;
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
    public:
        GasLiftStage2(
            const BlackoilWellModel &well_model,
            const Simulator &ebos_simulator,
            DeferredLogger &deferred_logger,
            WellState &well_state,
            GLiftProdWells &prod_wells,
            GLiftOptWells &glift_wells,
            GLiftWellStateMap &state_map
        );
        void runOptimize();
    private:
        void addOrRemoveALQincrement_(
            GradMap &grad_map, const std::string well_name, bool add);
        std::optional<GradInfo> calcIncOrDecGrad_(
            const std::string name, const GasLiftSingleWell &gs_well, bool increase);
        bool checkRateAlreadyLimited_(GLiftWellState &state, bool increase);
        GradInfo deleteDecGradItem_(const std::string &name);
        GradInfo deleteIncGradItem_(const std::string &name);
        GradInfo deleteGrad_(const std::string &name, bool increase);
        void displayDebugMessage_(const std::string &msg);
        void displayDebugMessage2B_(const std::string &msg);
        void displayDebugMessage_(const std::string &msg, const std::string &group_name);
        void displayWarning_(const std::string &msg, const std::string &group_name);
        void displayWarning_(const std::string &msg);
        std::tuple<double, double, double> getCurrentGroupRates_(
            const Opm::Group &group);
        std::tuple<double, double, double> getCurrentGroupRatesRecursive_(
            const Opm::Group &group);
        std::tuple<double, double, double> getCurrentWellRates_(
            const std::string &well_name, const std::string &group_name);
        std::vector<GasLiftSingleWell *> getGroupGliftWells_(
            const Opm::Group &group);
        void getGroupGliftWellsRecursive_(
            const Opm::Group &group, std::vector<GasLiftSingleWell *> &wells);
        std::pair<double, double> getStdWellRates_(const WellInterface<TypeTag> &well);
        void optimizeGroup_(const Opm::Group &group);
        void optimizeGroupsRecursive_(const Opm::Group &group);
        void recalculateGradientAndUpdateData_(
            GradPairItr &grad_itr, bool increase,
            std::vector<GradPair> &grads, std::vector<GradPair> &other_grads);
        void redistributeALQ_(
            std::vector<GasLiftSingleWell *> &wells,  const Opm::Group &group,
            std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads);
        void removeSurplusALQ_(
            const Opm::Group &group,
            std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads);
        void saveGrad_(GradMap &map, const std::string &name, GradInfo &grad);
        void saveDecGrad_(const std::string &name, GradInfo &grad);
        void saveIncGrad_(const std::string &name, GradInfo &grad);
        void sortGradients_(std::vector<GradPair> &grads);
        std::optional<GradInfo> updateGrad_(
            const std::string &name, GradInfo &grad, bool increase);
        void updateGradVector_(
            const std::string &name, std::vector<GradPair> &grads, double grad);

        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        const BlackoilWellModel &well_model_;
        WellState &well_state_;
        GLiftProdWells &prod_wells_;
        GLiftOptWells &stage1_wells_;
        GLiftWellStateMap &well_state_map_;

        int report_step_idx_;
        const SummaryState &summary_state_;
        const Schedule &schedule_;
        const PhaseUsage &phase_usage_;
        const GasLiftOpt& glo_;
        GradMap inc_grads_;
        GradMap dec_grads_;
        bool debug_;
        int max_iterations_ = 1000;
        //int time_step_idx_;
        int nonlinear_iteration_idx_;

        struct OptimizeState {
            OptimizeState( GasLiftStage2 &parent_, const Opm::Group &group_ ) :
                parent{parent_},
                group{group_},
                it{0}
            {}
            GasLiftStage2 &parent;
            const Opm::Group &group;
            int it;

            using GradInfo = typename GasLiftStage2::GradInfo;
            using GradPair = typename GasLiftStage2::GradPair;
            using GradPairItr = typename GasLiftStage2::GradPairItr;
            using GradMap = typename GasLiftStage2::GradMap;
            void calculateEcoGradients(std::vector<GasLiftSingleWell *> &wells,
                std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads);
            bool checkAtLeastTwoWells(std::vector<GasLiftSingleWell *> &wells);
            void debugShowIterationInfo();
            std::pair<std::optional<GradPairItr>,std::optional<GradPairItr>>
               getEcoGradients(
                   std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads);
            void recalculateGradients(
                std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads,
                GradPairItr &min_dec_grad_itr, GradPairItr &max_inc_grad_itr);
            void redistributeALQ( GradPairItr &min_dec_grad, GradPairItr &max_inc_grad);

        private:
            void displayDebugMessage_(const std::string &msg);
            void displayWarning_(const std::string &msg);

        };

        struct SurplusState {
            SurplusState( GasLiftStage2 &parent_, const Opm::Group &group_,
                double oil_rate_, double gas_rate_, double alq_, double min_eco_grad_,
                double oil_target_, double gas_target_,
                std::optional<double> max_glift_) :
                parent{parent_},
                group{group_},
                oil_rate{oil_rate_},
                gas_rate{gas_rate_},
                alq{alq_},
                min_eco_grad{min_eco_grad_},
                oil_target{oil_target_},
                gas_target{gas_target_},
                max_glift{max_glift_},
                it{0}
            {}
            GasLiftStage2 &parent;
            const Opm::Group &group;
            double oil_rate;
            double gas_rate;
            double alq;
            const double min_eco_grad;
            const double oil_target;
            const double gas_target;
            std::optional<double> max_glift;
            int it;

            void addOrRemoveALQincrement(
                GradMap &grad_map, const std::string well_name, bool add);
            bool checkALQlimit();
            bool checkEcoGradient(const std::string &well_name, double eco_grad);
            bool checkGasTarget();
            bool checkOilTarget();
            void updateRates(const std::string &name, const GradInfo &gi);
        private:
        };
    };

#include "GasLiftStage2_impl.hpp"

} // namespace Opm

#endif // OPM_GASLIFT_STAGE2_HEADER_INCLUDED
