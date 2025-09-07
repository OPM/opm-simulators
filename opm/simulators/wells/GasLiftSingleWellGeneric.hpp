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

#include <opm/input/eclipse/Schedule/Well/WellProductionControls.hpp>

#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftCommon.hpp>

#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace Opm {

class DeferredLogger;
class GasLiftWell;
template<class Scalar> class GasLiftWellState;
class Schedule;
class SummaryState;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
template<typename Scalar, typename IndexTraits> class WellState;
template<class Scalar> class GroupState;

template<typename Scalar, typename IndexTraits>
class GasLiftSingleWellGeneric : public GasLiftCommon<Scalar, IndexTraits>
{
protected:
    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;
    static constexpr int NUM_PHASES = 3;
    static constexpr Scalar ALQ_EPSILON = 1e-8;

public:
    using GLiftSyncGroups = std::set<int>;
    using Rate = typename GasLiftGroupInfo<Scalar, IndexTraits>::Rate;
    using MessageType = typename GasLiftCommon<Scalar, IndexTraits>::MessageType;

    struct GradInfo
    {
        GradInfo() = default;
        GradInfo(Scalar grad_,
                 Scalar new_oil_rate_,
                 bool oil_is_limited_,
                 Scalar new_gas_rate_,
                 bool gas_is_limited_,
                 Scalar new_water_rate_,
                 bool water_is_limited_,
                 Scalar alq_,
                 bool alq_is_limited_)
            : grad{grad_}
            , new_oil_rate{new_oil_rate_}
            , oil_is_limited{oil_is_limited_}
            , new_gas_rate{new_gas_rate_}
            , gas_is_limited{gas_is_limited_}
            , new_water_rate{new_water_rate_}
            , water_is_limited{water_is_limited_}
            , alq{alq_}
            , alq_is_limited{alq_is_limited_}
        {}

        Scalar grad;
        Scalar new_oil_rate;
        bool oil_is_limited;
        Scalar new_gas_rate;
        bool gas_is_limited;
        Scalar new_water_rate;
        bool water_is_limited;
        Scalar alq;
        bool alq_is_limited;
    };

    const std::string& name() const { return well_name_; }

    std::optional<GradInfo> calcIncOrDecGradient(Scalar oil_rate,
                                                 Scalar gas_rate,
                                                 Scalar water_rate,
                                                 Scalar alq,
                                                 const std::string& gr_name_dont_limit,
                                                 bool increase,
                                                 bool debug_output = true) const;

    std::unique_ptr<GasLiftWellState<Scalar>> runOptimize(const int iteration_idx);

    std::pair<Scalar, bool> wellTestALQ();

    virtual const WellInterfaceGeneric<Scalar, IndexTraits>& getWell() const = 0;

protected:
    GasLiftSingleWellGeneric(DeferredLogger& deferred_logger,
                             WellState<Scalar, IndexTraits>& well_state,
                             const GroupState<Scalar>& group_state,
                             const Well& ecl_well,
                             const SummaryState& summary_state,
                             GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                             const Schedule& schedule,
                             const int report_step_idx,
                             GLiftSyncGroups& sync_groups,
                             const Parallel::Communication& comm,
                             bool glift_debug);

    struct LimitedRates;
    struct BasicRates
    {
        BasicRates(const BasicRates& rates) :
            oil{rates.oil},
            gas{rates.gas},
            water{rates.water},
            bhp_is_limited{rates.bhp_is_limited}
        {}

        BasicRates(Scalar oil_,
                   Scalar gas_,
                   Scalar water_,
                   bool bhp_is_limited_)
            : oil{oil_}
            , gas{gas_}
            , water{water_}
            , bhp_is_limited{bhp_is_limited_}
        {}

        BasicRates& operator=(const BasicRates& rates)
        {
            oil = rates.oil;
            gas = rates.gas;
            water = rates.water;
            bhp_is_limited = rates.bhp_is_limited;
            return *this;
        }

        // This copy constructor cannot be defined inline here since LimitedRates
        //   has not been defined yet (it is defined below). Instead it is defined in
        //   in the .cpp file
        explicit BasicRates(const LimitedRates& rates);

        Scalar operator[](Rate rate_type) const
        {
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

        Scalar oil, gas, water;
        bool bhp_is_limited;
    };

    struct LimitedRates : public BasicRates
    {
        enum class LimitType {well, group, none};
        LimitedRates(Scalar oil_,
                     Scalar gas_,
                     Scalar water_,
                     bool oil_is_limited_,
                     bool gas_is_limited_,
                     bool water_is_limited_,
                     bool bhp_is_limited_,
                     std::optional<Rate> oil_limiting_target_,
                     std ::optional<Rate> water_limiting_target_)
            :  BasicRates(oil_, gas_, water_, bhp_is_limited_)
            , oil_is_limited{oil_is_limited_}
            , gas_is_limited{gas_is_limited_}
            , water_is_limited{water_is_limited_}
            , oil_limiting_target{oil_limiting_target_}
            , water_limiting_target{water_limiting_target_}
        {
            set_initial_limit_type_();
        }

        LimitedRates(const BasicRates& rates,
                     bool oil_is_limited_,
                     bool gas_is_limited_,
                     bool water_is_limited_)
            : BasicRates(rates)
            , oil_is_limited{oil_is_limited_}
            , gas_is_limited{gas_is_limited_}
            , water_is_limited{water_is_limited_}
        {
            set_initial_limit_type_();
        }

        bool limited() const
        {
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
        void set_initial_limit_type_()
        {
            limit_type = limited() ? LimitType::well : LimitType::none;
        }
    };

    struct OptimizeState
    {
        OptimizeState( GasLiftSingleWellGeneric& parent_, bool increase_ )
            : parent{parent_}
            , increase{increase_}
            , it{0}
            , stop_iteration{false}
            , bhp{-1}
        {}

        GasLiftSingleWellGeneric& parent;
        bool increase;
        int it;
        bool stop_iteration;
        Scalar bhp;

        std::pair<std::optional<Scalar>,bool> addOrSubtractAlqIncrement(Scalar alq);
        Scalar calcEcoGradient(Scalar oil_rate,
                               Scalar new_oil_rate,
                               Scalar gas_rate,
                               Scalar new_gas_rate);

        bool checkAlqOutsideLimits(Scalar alq, Scalar oil_rate);
        bool checkEcoGradient(Scalar gradient);
        bool checkOilRateExceedsTarget(Scalar oil_rate);
        bool checkRatesViolated(const LimitedRates& rates) const;

        void debugShowIterationInfo(Scalar alq);

        Scalar getBhpWithLimit();

        void warn_(const std::string& msg) { parent.displayWarning_(msg); }
    };

    bool checkGroupALQrateExceeded(Scalar delta_alq,
                                   const std::string& gr_name_dont_limit = "") const;
    bool checkGroupTotalRateExceeded(Scalar delta_alq,
                                     Scalar delta_gas_rate) const;

    std::pair<std::optional<Scalar>, bool>
    addOrSubtractAlqIncrement_(Scalar alq, bool increase) const;

    Scalar calcEcoGradient_(Scalar oil_rate, Scalar new_oil_rate,
                            Scalar gas_rate, Scalar new_gas_rate, bool increase) const;

    bool checkALQequal_(Scalar alq1, Scalar alq2) const;

    bool checkGroupTargetsViolated(const BasicRates& rates,
                                   const BasicRates& new_rates) const;
    bool checkInitialALQmodified_(Scalar alq, Scalar initial_alq) const;

    virtual bool checkThpControl_() const = 0;
    virtual std::optional<Scalar > computeBhpAtThpLimit_(Scalar alq,
                                                        bool debug_output = true) const = 0;

    std::pair<std::optional<Scalar>,Scalar>
    computeConvergedBhpAtThpLimitByMaybeIncreasingALQ_() const;

    std::pair<std::optional<BasicRates>,Scalar>
    computeInitialWellRates_() const;

    std::optional<LimitedRates>
    computeLimitedWellRatesWithALQ_(Scalar alq) const;

    virtual BasicRates computeWellRates_(Scalar bhp,
                                         bool bhp_is_limited,
                                         bool debug_output = true) const = 0;

    std::optional<BasicRates> computeWellRatesWithALQ_(Scalar alq) const;

    void debugCheckNegativeGradient_(Scalar grad, Scalar alq, Scalar new_alq,
                                     Scalar oil_rate, Scalar new_oil_rate,
                                     Scalar gas_rate, Scalar new_gas_rate,
                                     bool increase) const;

    void debugPrintWellStateRates() const;
    void debugShowAlqIncreaseDecreaseCounts_();
    void debugShowBhpAlqTable_();
    void debugShowLimitingTargets_(const LimitedRates& rates) const;
    void debugShowProducerControlMode() const;
    void debugShowStartIteration_(Scalar alq, bool increase, Scalar oil_rate);
    void debugShowTargets_();
    void displayDebugMessage_(const std::string& msg) const override;
    void displayWarning_(const std::string& warning);

    std::pair<Scalar, bool> getBhpWithLimit_(Scalar bhp) const;
    std::pair<Scalar, bool> getGasRateWithLimit_(const BasicRates& rates) const;
    std::pair<Scalar, bool> getGasRateWithGroupLimit_(Scalar new_gas_rate,
                                                      Scalar gas_rate,
                                                      const std::string& gr_name_dont_limit) const;

    std::pair<std::optional<LimitedRates>,Scalar >
    getInitialRatesWithLimit_() const;

    LimitedRates getLimitedRatesFromRates_(const BasicRates& rates) const;

    std::tuple<Scalar,Scalar,bool,bool>
    getLiquidRateWithGroupLimit_(const Scalar new_oil_rate,
                                 const Scalar oil_rate,
                                 const Scalar new_water_rate,
                                 const Scalar water_rate,
                                 const std::string& gr_name_dont_limit) const;

    std::pair<Scalar, bool>
    getOilRateWithGroupLimit_(Scalar new_oil_rate,
                              Scalar oil_rate,
                              const std::string& gr_name_dont_limit) const;

    std::pair<Scalar, bool> getOilRateWithLimit_(const BasicRates& rates) const;

    std::pair<Scalar, std::optional<Rate>>
    getOilRateWithLimit2_(const BasicRates& rates) const;

    Scalar getProductionTarget_(Rate rate) const;
    Scalar getRate_(Rate rate_type, const BasicRates& rates) const;

    std::pair<Scalar, std::optional<Rate>>
    getRateWithLimit_(Rate rate_type, const BasicRates& rates) const;

    std::tuple<Scalar, const std::string*, Scalar>
    getRateWithGroupLimit_(Rate rate_type,
                           const Scalar new_rate,
                           const Scalar old_rate,
                           const std::string& gr_name_dont_limit) const;

    std::pair<Scalar, bool>
    getWaterRateWithGroupLimit_(Scalar new_water_rate,
                                Scalar water_rate,
                                const std::string& gr_name_dont_limit) const;

    std::pair<Scalar, bool> getWaterRateWithLimit_(const BasicRates& rates) const;

    std::pair<Scalar, std::optional<Rate>>
    getWaterRateWithLimit2_(const BasicRates& rates) const;

    BasicRates getWellStateRates_() const;
    bool hasProductionControl_(Rate rate) const;

    std::pair<LimitedRates, Scalar>
    increaseALQtoPositiveOilRate_(Scalar alq,
                                  const LimitedRates& orig_rates) const;

    std::pair<LimitedRates, Scalar>
    increaseALQtoMinALQ_(Scalar alq,
                         const LimitedRates& orig_rates) const;

    void logSuccess_(Scalar alq,
                     const int iteration_idx);

    std::pair<LimitedRates, Scalar>
    maybeAdjustALQbeforeOptimizeLoop_(const LimitedRates& rates,
                                      Scalar alq,
                                      bool increase) const;

    std::pair<LimitedRates, Scalar>
    reduceALQtoGroupAlqLimits_(Scalar alq,
                               const LimitedRates& rates) const;

    std::pair<LimitedRates, Scalar>
    reduceALQtoGroupTarget(Scalar alq,
                           const LimitedRates& rates) const;

    std::pair<LimitedRates, Scalar>
    reduceALQtoWellTarget_(Scalar alq,
                           const LimitedRates& rates) const;

    std::unique_ptr<GasLiftWellState<Scalar>> runOptimize1_();
    std::unique_ptr<GasLiftWellState<Scalar>> runOptimize2_();
    std::unique_ptr<GasLiftWellState<Scalar>> runOptimizeLoop_(bool increase);

    void setAlqMinRate_(const GasLiftWell& well);
    std::unique_ptr<GasLiftWellState<Scalar>> tryIncreaseLiftGas_();
    std::unique_ptr<GasLiftWellState<Scalar>> tryDecreaseLiftGas_();

    void updateGroupRates_(const LimitedRates& rates,
                           const LimitedRates& new_rates,
                           Scalar delta_alq) const;

    LimitedRates
    updateRatesToGroupLimits_(const BasicRates& old_rates,
                              const LimitedRates& rates,
                              const std::string& gr_name = "") const;

    void updateWellStateAlqFixedValue_(const GasLiftWell& well);
    bool useFixedAlq_(const GasLiftWell& well);

    void debugInfoGroupRatesExceedTarget(Rate rate_type,
                                         const std::string& gr_name,
                                         Scalar rate,
                                         Scalar target) const;
    void warnMaxIterationsExceeded_();

    const Well& ecl_well_;
    const SummaryState& summary_state_;
    GasLiftGroupInfo<Scalar, IndexTraits>& group_info_;
    GLiftSyncGroups& sync_groups_;
    const WellProductionControls controls_;

    Scalar increment_;
    Scalar max_alq_;
    Scalar min_alq_;
    Scalar orig_alq_;

    Scalar alpha_w_;
    Scalar alpha_g_;
    Scalar eco_grad_;

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
