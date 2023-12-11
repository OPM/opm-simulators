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

#ifndef OPM_GASLIFT_SINGLE_WELL_GENERIC_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_GENERIC_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/Well/WellProductionControls.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftCommon.hpp>

#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace Opm
{

class DeferredLogger;
class GasLiftWell;
class GasLiftWellState;
class Schedule;
class SummaryState;
class WellInterfaceGeneric;
class WellState;
class GroupState;

class GasLiftSingleWellGeneric : public GasLiftCommon
{
protected:
    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int NUM_PHASES = 3;
    static constexpr double ALQ_EPSILON = 1e-8;

public:
    using GLiftSyncGroups = std::set<int>;
    using Rate = GasLiftGroupInfo::Rate;
    struct GradInfo
    {
        GradInfo() { }

        GradInfo(double grad_, double new_oil_rate_, bool oil_is_limited_,
                 double new_gas_rate_, bool gas_is_limited_,
                 double new_water_rate_, bool water_is_limited_,
                 double alq_, bool alq_is_limited_) :
            grad{grad_},
            new_oil_rate{new_oil_rate_},
            oil_is_limited{oil_is_limited_},
            new_gas_rate{new_gas_rate_},
            gas_is_limited{gas_is_limited_},
            new_water_rate{new_water_rate_},
            water_is_limited{water_is_limited_},
            alq{alq_},
            alq_is_limited{alq_is_limited_} {}
        double grad;
        double new_oil_rate;
        bool oil_is_limited;
        double new_gas_rate;
        bool gas_is_limited;
        double new_water_rate;
        bool water_is_limited;
        double alq;
        bool alq_is_limited;
    };

    virtual ~GasLiftSingleWellGeneric() = default;

    const std::string& name() const { return well_name_; }

    std::optional<GradInfo> calcIncOrDecGradient(double oil_rate, double gas_rate,
                                                 double water_rate,
                                                 double alq,
                                                 const std::string& gr_name_dont_limit,
                                                 bool increase,
                                                 bool debug_output = true
                                                 ) const;

    std::unique_ptr<GasLiftWellState> runOptimize(const int iteration_idx);

    virtual const WellInterfaceGeneric& getWell() const = 0;

protected:
    GasLiftSingleWellGeneric(
        DeferredLogger& deferred_logger,
        WellState& well_state,
        const GroupState& group_state,
        const Well& ecl_well,
        const SummaryState& summary_state,
        GasLiftGroupInfo& group_info,
        const PhaseUsage& phase_usage,
        const Schedule& schedule,
        const int report_step_idx,
        GLiftSyncGroups& sync_groups,
        const Parallel::Communication& comm,
        bool glift_debug
    );

    struct LimitedRates;
    struct BasicRates
    {
        BasicRates(const BasicRates& rates) :
            oil{rates.oil},
            gas{rates.gas},
            water{rates.water},
            bhp_is_limited{rates.bhp_is_limited}
        {}
        BasicRates(double oil_, double gas_, double water_, bool bhp_is_limited_) :
            oil{oil_},
            gas{gas_},
            water{water_},
            bhp_is_limited{bhp_is_limited_}
        {}
        BasicRates& operator=(const BasicRates& rates) {
            oil = rates.oil;
            gas = rates.gas;
            water = rates.water;
            bhp_is_limited = rates.bhp_is_limited;
            return *this;
        }
        // This copy constructor cannot be defined inline here since LimitedRates
        //   has not been defined yet (it is defined below). Instead it is defined in
        //   in the .cpp file
        BasicRates(const LimitedRates& rates);
        double operator[](Rate rate_type) const {
            switch (rate_type) {
            case Rate::oil:
                return this->oil;
            case Rate::gas:
                return this->gas;
            case Rate::water:
                return this->water;
            case Rate::liquid:
                return this->oil + this->water;
            default:
                throw std::runtime_error("This should not happen");
            }
        }

        double oil, gas, water;
        bool bhp_is_limited;
    };

    struct LimitedRates : public BasicRates
    {
        enum class LimitType {well, group, none};
        LimitedRates(
              double oil_, double gas_, double water_,
              bool oil_is_limited_, bool gas_is_limited_,
              bool water_is_limited_, bool bhp_is_limited_,
              std::optional<Rate> oil_limiting_target_,
              std::optional<Rate> water_limiting_target_
        ) :
            BasicRates(oil_, gas_, water_, bhp_is_limited_),
            oil_is_limited{oil_is_limited_},
            gas_is_limited{gas_is_limited_},
            water_is_limited{water_is_limited_},
            oil_limiting_target{oil_limiting_target_},
            water_limiting_target{water_limiting_target_}
        {
            set_initial_limit_type_();
        }

        LimitedRates(
              const BasicRates& rates,
              bool oil_is_limited_, bool gas_is_limited_,
              bool water_is_limited_
        ) :
            BasicRates(rates),
            oil_is_limited{oil_is_limited_},
            gas_is_limited{gas_is_limited_},
            water_is_limited{water_is_limited_}
        {
            set_initial_limit_type_();
        }

        bool limited() const {
            return oil_is_limited || gas_is_limited || water_is_limited;
        }
        // For a given ALQ value, were the rates limited due to group targets
        //   or due to well targets?
        LimitType limit_type;
        bool oil_is_limited;
        bool gas_is_limited;
        bool water_is_limited;
        std::optional<Rate> oil_limiting_target;
        std::optional<Rate> water_limiting_target;
    private:
        void set_initial_limit_type_() {
            limit_type = limited() ? LimitType::well : LimitType::none;
        }
    };

    struct OptimizeState
    {
        OptimizeState( GasLiftSingleWellGeneric& parent_, bool increase_ ) :
            parent{parent_},
            increase{increase_},
            it{0},
            stop_iteration{false},
            bhp{-1}
        {}

        GasLiftSingleWellGeneric& parent;
        bool increase;
        int it;
        bool stop_iteration;
        double bhp;

        std::pair<std::optional<double>,bool> addOrSubtractAlqIncrement(double alq);
        double calcEcoGradient(double oil_rate, double new_oil_rate,
                               double gas_rate, double new_gas_rate);
        bool checkAlqOutsideLimits(double alq, double oil_rate);
        bool checkEcoGradient(double gradient);
        bool checkOilRateExceedsTarget(double oil_rate);
        bool checkRatesViolated(const LimitedRates& rates) const;
        void debugShowIterationInfo(double alq);
        double getBhpWithLimit();
        void warn_(std::string msg) {parent.displayWarning_(msg);}
    };
    bool checkGroupALQrateExceeded(double delta_alq, const std::string& gr_name_dont_limit = "") const;
    bool checkGroupTotalRateExceeded(double delta_alq, double delta_gas_rate) const;

    std::pair<std::optional<double>, bool> addOrSubtractAlqIncrement_(
                            double alq, bool increase) const;
    double calcEcoGradient_(double oil_rate, double new_oil_rate,
                            double gas_rate, double new_gas_rate, bool increase) const;
    bool checkALQequal_(double alq1, double alq2) const;
    bool checkGroupTargetsViolated(
                      const BasicRates& rates, const BasicRates& new_rates) const;
    bool checkInitialALQmodified_(double alq, double initial_alq) const;
    virtual bool checkThpControl_() const = 0;
    virtual std::optional<double> computeBhpAtThpLimit_(double alq, bool debug_output = true) const = 0;
    std::pair<std::optional<double>,double> computeConvergedBhpAtThpLimitByMaybeIncreasingALQ_() const;
    std::pair<std::optional<BasicRates>,double> computeInitialWellRates_() const;
    std::optional<LimitedRates> computeLimitedWellRatesWithALQ_(double alq) const;
    virtual BasicRates computeWellRates_(double bhp, bool bhp_is_limited, bool debug_output = true) const = 0;
    std::optional<BasicRates> computeWellRatesWithALQ_(double alq) const;
    void debugCheckNegativeGradient_(double grad, double alq, double new_alq,
                                     double oil_rate, double new_oil_rate,
                                     double gas_rate, double new_gas_rate,
                                     bool increase) const;
    void debugPrintWellStateRates() const;
    void debugShowAlqIncreaseDecreaseCounts_();
    void debugShowBhpAlqTable_();
    void debugShowLimitingTargets_(const LimitedRates& rates) const;
    void debugShowProducerControlMode() const;
    void debugShowStartIteration_(double alq, bool increase, double oil_rate);
    void debugShowTargets_();
    void displayDebugMessage_(const std::string& msg) const override;
    void displayWarning_(const std::string& warning);
    std::pair<double, bool> getBhpWithLimit_(double bhp) const;
    std::pair<double, bool> getGasRateWithLimit_(
                           const BasicRates& rates) const;
    std::pair<double, bool> getGasRateWithGroupLimit_(
                           double new_gas_rate, double gas_rate, const std::string& gr_name_dont_limit) const;
    std::pair<std::optional<LimitedRates>,double> getInitialRatesWithLimit_() const;
    LimitedRates getLimitedRatesFromRates_(const BasicRates& rates) const;
    std::tuple<double,double,bool,bool> getLiquidRateWithGroupLimit_(
                           const double new_oil_rate, const double oil_rate,
                           const double new_water_rate, const double water_rate, const std::string& gr_name_dont_limit) const;
    std::pair<double, bool> getOilRateWithGroupLimit_(
                           double new_oil_rate, double oil_rate, const std::string& gr_name_dont_limit) const;
    std::pair<double, bool> getOilRateWithLimit_(const BasicRates& rates) const;
    std::pair<double, std::optional<Rate>> getOilRateWithLimit2_(
                           const BasicRates& rates) const;
    double getProductionTarget_(Rate rate) const;
    double getRate_(Rate rate_type, const BasicRates& rates) const;
    std::pair<double, std::optional<Rate>> getRateWithLimit_(
                    Rate rate_type, const BasicRates& rates) const;
    std::tuple<double, const std::string*, double> getRateWithGroupLimit_(
                 Rate rate_type, const double new_rate, const double old_rate, const std::string& gr_name_dont_limit) const;
    std::pair<double, bool> getWaterRateWithGroupLimit_(
                           double new_water_rate, double water_rate, const std::string& gr_name_dont_limit) const;
    std::pair<double, bool> getWaterRateWithLimit_(const BasicRates& rates) const;
    std::pair<double, std::optional<Rate>> getWaterRateWithLimit2_(
                           const BasicRates& rates) const;
    BasicRates getWellStateRates_() const;
    bool hasProductionControl_(Rate rate) const;
    std::pair<LimitedRates, double> increaseALQtoPositiveOilRate_(
                           double alq, const LimitedRates& orig_rates) const;
    std::pair<LimitedRates, double> increaseALQtoMinALQ_(
                           double alq, const LimitedRates& orig_rates) const;
    void logSuccess_(double alq,
                     const int iteration_idx);
    std::pair<LimitedRates, double> maybeAdjustALQbeforeOptimizeLoop_(
                        const LimitedRates& rates, double alq, bool increase) const;
    std::pair<LimitedRates, double> reduceALQtoGroupAlqLimits_(
                        double alq, const LimitedRates& rates) const;
    std::pair<LimitedRates, double> reduceALQtoGroupTarget(
                        double alq, const LimitedRates& rates) const;
    std::pair<LimitedRates, double> reduceALQtoWellTarget_(
                        double alq, const LimitedRates& rates) const;
    std::unique_ptr<GasLiftWellState> runOptimize1_();
    std::unique_ptr<GasLiftWellState> runOptimize2_();
    std::unique_ptr<GasLiftWellState> runOptimizeLoop_(bool increase);
    void setAlqMinRate_(const GasLiftWell& well);
    std::unique_ptr<GasLiftWellState> tryIncreaseLiftGas_();
    std::unique_ptr<GasLiftWellState> tryDecreaseLiftGas_();
    void updateGroupRates_(
        const LimitedRates& rates,
        const LimitedRates& new_rates,
        double delta_alq) const;
    LimitedRates updateRatesToGroupLimits_(
        const BasicRates& rates, const LimitedRates& new_rates, const std::string& gr_name = "") const;
    void updateWellStateAlqFixedValue_(const GasLiftWell& well);
    bool useFixedAlq_(const GasLiftWell& well);
    void debugInfoGroupRatesExceedTarget(
        Rate rate_type, const std::string& gr_name, double rate, double target) const;
    void warnMaxIterationsExceeded_();

    const Well& ecl_well_;
    const SummaryState& summary_state_;
    GasLiftGroupInfo& group_info_;
    const PhaseUsage& phase_usage_;
    GLiftSyncGroups& sync_groups_;
    const WellProductionControls controls_;

    double increment_;
    double max_alq_;
    double min_alq_;
    double orig_alq_;

    double alpha_w_;
    double alpha_g_;
    double eco_grad_;

    int gas_pos_;
    int oil_pos_;
    int water_pos_;

    int max_iterations_;

    std::string well_name_;

    const GasLiftWell* gl_well_;

    bool optimize_;
    bool debug_limit_increase_decrease_;
    bool debug_abort_if_decrease_and_oil_is_limited_ = false;
    bool debug_abort_if_increase_and_gas_is_limited_ = false;
};

} // namespace Opm

#endif // OPM_GASLIFT_SINGLE_WELL_GENERIC_HEADER_INCLUDED
