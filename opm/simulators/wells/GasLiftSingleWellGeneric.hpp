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
                 Scalar new_oil_pot_,
                 bool oil_is_limited_,
                 Scalar new_gas_rate_,
                 Scalar new_gas_pot_,
                 bool gas_is_limited_,
                 Scalar new_water_rate_,
                 Scalar new_water_pot_,
                 bool water_is_limited_,
                 Scalar alq_,
                 bool alq_is_limited_,
                 Scalar bhp_)
            : grad{grad_}
            , new_oil_rate{new_oil_rate_}
            , new_oil_pot{new_oil_pot_}
            , oil_is_limited{oil_is_limited_}
            , new_gas_rate{new_gas_rate_}
            , new_gas_pot{new_gas_pot_}
            , gas_is_limited{gas_is_limited_}
            , new_water_rate{new_water_rate_}
            , new_water_pot{new_water_pot_}
            , water_is_limited{water_is_limited_}
            , alq{alq_}
            , alq_is_limited{alq_is_limited_}
            , bhp{bhp_}
        {}

        Scalar grad;
        Scalar new_oil_rate;
        Scalar new_oil_pot;
        bool oil_is_limited;
        Scalar new_gas_rate;
        Scalar new_gas_pot;
        bool gas_is_limited;
        Scalar new_water_rate;
        Scalar new_water_pot;
        bool water_is_limited;
        Scalar alq;
        bool alq_is_limited;
        Scalar bhp;
    };

    const std::string& name() const { return well_name_; }

    std::optional<GradInfo> calcIncOrDecGradient(const GasLiftWellState<Scalar>& state,
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

    struct LimitedRatesAndBhp;
    struct RatesAndBhp
    {
        RatesAndBhp(const RatesAndBhp& rates) :
            oil{rates.oil},
            gas{rates.gas},
            water{rates.water},
            bhp{rates.bhp},
            bhp_is_limited{rates.bhp_is_limited}
        {}

        RatesAndBhp(Scalar oil_,
                   Scalar gas_,
                   Scalar water_,
                   Scalar bhp_,
                   bool bhp_is_limited_)
            : oil{oil_}
            , gas{gas_}
            , water{water_}
            , bhp{bhp_}
            , bhp_is_limited{bhp_is_limited_}
        {}

        RatesAndBhp& operator=(const RatesAndBhp& rates)
        {
            oil = rates.oil;
            gas = rates.gas;
            water = rates.water;
            bhp = rates.bhp;
            bhp_is_limited = rates.bhp_is_limited;
            return *this;
        }

        // This copy constructor cannot be defined inline here since LimitedRatesAndBhp
        //   has not been defined yet (it is defined below). Instead it is defined in
        //   in the .cpp file
        explicit RatesAndBhp(const LimitedRatesAndBhp& rates);

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

        Scalar oil, gas, water, bhp;
        bool bhp_is_limited;
    };

    struct LimitedRatesAndBhp : public RatesAndBhp
    {
        enum class LimitType {well, group, none};
        LimitedRatesAndBhp(Scalar oil_,
                     Scalar oil_pot_,
                     Scalar gas_,
                     Scalar gas_pot_,
                     Scalar water_,
                     Scalar water_pot_,
                     Scalar bhp_,
                     bool oil_is_limited_,
                     bool gas_is_limited_,
                     bool water_is_limited_,
                     bool bhp_is_limited_)
            :  RatesAndBhp(oil_, gas_, water_, bhp_, bhp_is_limited_)
            , oil_pot(oil_pot_)
            , gas_pot(gas_pot_)
            , water_pot(water_pot_)
            , oil_is_limited{oil_is_limited_}
            , gas_is_limited{gas_is_limited_}
            , water_is_limited{water_is_limited_}
        {
            set_initial_limit_type_();
        }

        LimitedRatesAndBhp(const RatesAndBhp& rates,
                     Scalar oil_pot_,
                     Scalar gas_pot_,
                     Scalar water_pot_,
                     bool oil_is_limited_,
                     bool gas_is_limited_,
                     bool water_is_limited_)
            : RatesAndBhp(rates)
            , oil_pot(oil_pot_)
            , gas_pot(gas_pot_)
            , water_pot(water_pot_)
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
        Scalar oil_pot;
        Scalar gas_pot;
        Scalar water_pot;
        bool oil_is_limited;
        bool gas_is_limited;
        bool water_is_limited;

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
        bool checkRatesViolated(const LimitedRatesAndBhp& rates) const;

        void debugShowIterationInfo(Scalar alq);

        Scalar getBhpWithLimit();

        void warn_(const std::string& msg) { parent.displayWarning_(msg); }
    };

    bool checkGroupALQrateExceeded(Scalar delta_alq,
                                   const std::string& gr_name_dont_limit = "") const;
    bool checkGroupTotalRateExceeded(Scalar delta_alq,
                                     Scalar delta_gas_rate,
                                     const std::string& gr_name_dont_limit = "") const;

    std::pair<std::optional<Scalar>, bool>
    addOrSubtractAlqIncrement_(Scalar alq, bool increase) const;

    Scalar calcEcoGradient_(Scalar oil_rate, Scalar new_oil_rate,
                            Scalar gas_rate, Scalar new_gas_rate, bool increase) const;

    bool checkALQequal_(Scalar alq1, Scalar alq2) const;

    bool checkGroupTargetsViolated(const RatesAndBhp& rates,
                                   const RatesAndBhp& new_rates) const;
    bool checkInitialALQmodified_(Scalar alq, Scalar initial_alq) const;

    virtual bool checkThpControl_() const = 0;
    virtual std::optional<Scalar > computeBhpAtThpLimit_(Scalar alq, Scalar current_bhp,
                                                        bool debug_output = true) const = 0;

    std::pair<std::optional<Scalar>,Scalar>
    computeConvergedBhpAtThpLimitByMaybeIncreasingALQ_() const;

    std::pair<std::optional<RatesAndBhp>,Scalar>
    computeInitialWellRates_() const;

    std::optional<LimitedRatesAndBhp>
    computeLimitedWellRatesWithALQ_(Scalar alq, Scalar bhp) const;

    virtual RatesAndBhp computeWellRates_(Scalar bhp,
                                         bool bhp_is_limited,
                                         bool debug_output = true) const = 0;

    std::optional<RatesAndBhp> computeWellRatesWithALQ_(Scalar alq, Scalar bhp) const;

    void debugCheckNegativeGradient_(Scalar grad, Scalar alq, Scalar new_alq,
                                     Scalar oil_rate, Scalar new_oil_rate,
                                     Scalar gas_rate, Scalar new_gas_rate,
                                     bool increase) const;

    void debugPrintWellStateRates() const;
    void debugShowAlqIncreaseDecreaseCounts_();
    void debugShowBhpAlqTable_();
    void debugShowLimitingTargets_(const LimitedRatesAndBhp& rates) const;
    void debugShowProducerControlMode() const;
    void debugShowStartIteration_(Scalar alq, bool increase, Scalar oil_rate);
    void debugShowTargets_();
    void displayDebugMessage_(const std::string& msg) const override;
    void displayWarning_(const std::string& warning);

    std::pair<Scalar, bool> getBhpWithLimit_(Scalar bhp) const;
    std::pair<Scalar, bool> getGasRateWithGroupLimit_(Scalar new_gas_rate,
                                                      Scalar gas_rate,
                                                      const std::string& gr_name_dont_limit) const;

    std::pair<std::optional<LimitedRatesAndBhp>,Scalar >
    getInitialRatesWithLimit_() const;

    LimitedRatesAndBhp getLimitedRatesAndBhp_(const RatesAndBhp& rates) const;


    Scalar getProductionTarget_(Rate rate) const;
    Scalar getRate_(Rate rate_type, const RatesAndBhp& rates) const;

    std::pair<Scalar, std::optional<Rate>>
    getRateWithLimit_(Rate rate_type, const RatesAndBhp& rates) const;

    std::tuple<Scalar, const std::string*>
    getRateWithGroupLimit_(Rate rate_type,
                           const Scalar new_rate,
                           const Scalar old_rate,
                           const std::string& gr_name_dont_limit) const;

    RatesAndBhp getWellStateRates_() const;
    bool hasProductionControl_(Rate rate) const;

    std::pair<LimitedRatesAndBhp, Scalar>
    increaseALQtoPositiveOilRate_(Scalar alq,
                                  const LimitedRatesAndBhp& orig_rates) const;

    std::pair<LimitedRatesAndBhp, Scalar>
    increaseALQtoMinALQ_(Scalar alq,
                         const LimitedRatesAndBhp& orig_rates) const;

    void logSuccess_(Scalar alq,
                     const int iteration_idx);

    std::pair<LimitedRatesAndBhp, Scalar>
    maybeAdjustALQbeforeOptimizeLoop_(const LimitedRatesAndBhp& rates,
                                      Scalar alq,
                                      bool increase) const;

    std::pair<LimitedRatesAndBhp, Scalar>
    reduceALQtoGroupAlqLimits_(Scalar alq,
                               const LimitedRatesAndBhp& rates) const;

    std::pair<LimitedRatesAndBhp, Scalar>
    reduceALQtoGroupTarget(Scalar alq,
                           const LimitedRatesAndBhp& rates) const;

    std::pair<LimitedRatesAndBhp, Scalar>
    reduceALQtoWellTarget_(Scalar alq,
                           const LimitedRatesAndBhp& rates) const;

    std::unique_ptr<GasLiftWellState<Scalar>> runOptimize1_();
    std::unique_ptr<GasLiftWellState<Scalar>> runOptimize2_();
    std::unique_ptr<GasLiftWellState<Scalar>> runOptimizeLoop_(bool increase);

    void setAlqMinRate_(const GasLiftWell& well);
    std::unique_ptr<GasLiftWellState<Scalar>> tryIncreaseLiftGas_();
    std::unique_ptr<GasLiftWellState<Scalar>> tryDecreaseLiftGas_();

    void updateGroupRates_(const LimitedRatesAndBhp& rates,
                           const LimitedRatesAndBhp& new_rates,
                           Scalar delta_alq) const;

    LimitedRatesAndBhp
    updateRatesToGroupLimits_(const RatesAndBhp& old_rates,
                              const LimitedRatesAndBhp& rates,
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
