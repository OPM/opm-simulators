/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
// NOTE: StandardWell.hpp includes ourself (GasLiftSingleWell.hpp), so we need
//   to forward declare StandardWell for it to be defined in this file.
namespace Opm {
    template<typename TypeTag> class StandardWell;
}
#include <opm/simulators/wells/StandardWell.hpp>

#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <fmt/format.h>

namespace Opm
{
    template<class TypeTag>
    class GasLiftSingleWell
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using WellState = WellStateFullyImplicitBlackoil;
        using StdWell = StandardWell<TypeTag>;
        using GLiftWellState = GasLiftWellState;
        // TODO: same definition with WellInterface, and
        //  WellStateFullyImplicitBlackoil eventually they should go
        //  to a common header file.
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;
        static constexpr double ALQ_EPSILON = 1e-8;
        struct OptimizeState;
        class Stage2State;
    public:
        GasLiftSingleWell(
            const StdWell &std_well,
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            WellState &well_state
        );
        struct GradInfo;
        std::optional<GradInfo> calcIncOrDecGradient(
            double oil_rate, double gas_rate, double alq, bool increase) const;
        std::pair<double, double> getStage2Rates();
        const WellInterface<TypeTag> &getStdWell() const { return std_well_; }
        bool hasStage2Rates();
        std::unique_ptr<GLiftWellState> runOptimize();
        const std::string& name() const {return well_name_; }
        void updateStage2State(const GradInfo &gi, bool increase);

        struct GradInfo
        {
            GradInfo() { }

            GradInfo(double grad_, double new_oil_rate_, bool oil_is_limited_,
                     double new_gas_rate_, bool gas_is_limited_,
                     double alq_, bool alq_is_limited_) :
                grad{grad_},
                new_oil_rate{new_oil_rate_},
                oil_is_limited{oil_is_limited_},
                new_gas_rate{new_gas_rate_},
                gas_is_limited{gas_is_limited_},
                alq{alq_},
                alq_is_limited{alq_is_limited_} {}
            double grad;
            double new_oil_rate;
            bool oil_is_limited;
            double new_gas_rate;
            bool gas_is_limited;
            double alq;
            bool alq_is_limited;
        };

    private:
        std::pair<std::optional<double>, bool> addOrSubtractAlqIncrement_(
            double alq, bool increase) const;
        double calcEcoGradient_(double oil_rate, double new_oil_rate,
            double gas_rate, double new_gas_rate, bool increase) const;
        bool checkALQequal_(double alq1, double alq2) const;
        bool checkInitialALQmodified_(double alq, double initial_alq) const;
        bool checkWellRatesViolated_(
            std::vector<double> &potentials,
            const std::function<bool(double, double, const std::string &)> &callback,
            bool increase
        );
        std::optional<double> computeBhpAtThpLimit_(double alq) const;
        bool computeInitialWellRates_(std::vector<double> &potentials);
        void computeWellRates_(
            double bhp, std::vector<double> &potentials, bool debug_output=true) const;
        void debugCheckNegativeGradient_(double grad, double alq, double new_alq,
            double oil_rate, double new_oil_rate, double gas_rate,
            double new_gas_rate, bool increase) const;
        void debugShowAlqIncreaseDecreaseCounts_();
        void debugShowBhpAlqTable_();
        void debugShowStartIteration_(double alq, bool increase, double oil_rate);
        void debugShowTargets_();
        void displayDebugMessage_(const std::string &msg) const;
        void displayWarning_(std::string warning);
        std::pair<double, bool> getBhpWithLimit_(double bhp) const;
        std::pair<double, bool> getGasRateWithLimit_(
            const std::vector<double> &potentials) const;
        std::tuple<double,double,bool,bool> getInitialRatesWithLimit_(
            const std::vector<double> &potentials);
        std::pair<double, bool> getOilRateWithLimit_(
            const std::vector<double> &potentials) const;
        std::tuple<double,double,bool,bool,double>
            increaseALQtoPositiveOilRate_(double alq, double oil_rate, double gas_rate,
              bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials);
        std::tuple<double,double,bool,bool,double>
            increaseALQtoMinALQ_(double alq, double oil_rate, double gas_rate,
              bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials);
        void logSuccess_(double alq);
        std::tuple<double,double,bool,bool,double>
        reduceALQtoOilTarget_(double alq, double oil_rate, double gas_rate,
            bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials);

        std::unique_ptr<GLiftWellState> runOptimize1_();
        std::unique_ptr<GLiftWellState> runOptimize2_();
        std::unique_ptr<GLiftWellState> runOptimizeLoop_(bool increase);
        void setAlqMaxRate_(const GasLiftOpt::Well &well);
        void setAlqMinRate_(const GasLiftOpt::Well &well);
        std::unique_ptr<GLiftWellState> tryIncreaseLiftGas_();
        std::unique_ptr<GLiftWellState> tryDecreaseLiftGas_();
        void updateWellStateAlqFixedValue_(const GasLiftOpt::Well &well);
        bool useFixedAlq_(const GasLiftOpt::Well &well);
        void warnMaxIterationsExceeded_();

        DeferredLogger &deferred_logger_;
        const Simulator &ebos_simulator_;
        const StdWell &std_well_;
        const SummaryState &summary_state_;
        WellState &well_state_;
        std::string well_name_;
        const Well &ecl_well_;
        const Well::ProductionControls controls_;
        int num_phases_;
        bool debug;  // extra debug output
        bool debug_limit_increase_decrease_;
        bool debug_abort_if_decrease_and_oil_is_limited_ = false;
        bool debug_abort_if_increase_and_gas_is_limited_ = false;

        double alpha_w_;
        double alpha_g_;
        double eco_grad_;
        int gas_pos_;
        bool has_run_init_ = false;
        double increment_;
        double max_alq_;
        int max_iterations_;
        double min_alq_;
        int oil_pos_;
        bool optimize_;
        double orig_alq_;
        int water_pos_;

        struct OptimizeState
        {
            OptimizeState( GasLiftSingleWell &parent_, bool increase_ ) :
                parent{parent_},
                increase{increase_},
                it{0},
                stop_iteration{false},
                bhp{-1}
            {}

            GasLiftSingleWell &parent;
            bool increase;
            int it;
            bool stop_iteration;
            double bhp;

            std::pair<std::optional<double>,bool> addOrSubtractAlqIncrement(double alq);
            double calcEcoGradient(double oil_rate, double new_oil_rate,
                double gas_rate, double new_gas_rate);
            bool checkAlqOutsideLimits(double alq, double oil_rate);
            bool checkEcoGradient(double gradient);
            bool checkNegativeOilRate(double oil_rate);
            bool checkOilRateExceedsTarget(double oil_rate);
            bool checkRate(double rate, double limit, const std::string &rate_str) const;
            bool checkWellRatesViolated(std::vector<double> &potentials);
            bool computeBhpAtThpLimit(double alq);
            void debugShowIterationInfo(double alq);
            double getBhpWithLimit();
            void warn_(std::string msg) {parent.displayWarning_(msg);}
        };

    };

#include "GasLiftSingleWell_impl.hpp"

} // namespace Opm


#endif // OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
