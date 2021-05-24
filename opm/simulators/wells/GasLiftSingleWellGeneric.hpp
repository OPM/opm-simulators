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

#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

#include <functional>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#include <utility>

namespace Opm
{

class DeferredLogger;
class GasLiftWellState;
class Schedule;
class SummaryState;
class WellState;

class GasLiftSingleWellGeneric
{
protected:
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;
    static constexpr double ALQ_EPSILON = 1e-8;

public:
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

    virtual ~GasLiftSingleWellGeneric() = default;

    const std::string& name() const { return well_name_; }

    std::optional<GradInfo> calcIncOrDecGradient(double oil_rate, double gas_rate,
                                                 double alq, bool increase) const;

    std::unique_ptr<GasLiftWellState> runOptimize(const int iteration_idx);

protected:
    GasLiftSingleWellGeneric(DeferredLogger &deferred_logger,
                             WellState &well_state,
                             const Well& ecl_well,
                             const SummaryState& summary_state,
                             const Schedule& schedule,
                             const int report_step_idx);

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
        bool checkNegativeOilRate(double oil_rate);
        bool checkOilRateExceedsTarget(double oil_rate);
        bool checkRate(double rate, double limit, const std::string &rate_str) const;
        bool checkWellRatesViolated(std::vector<double> &potentials);
        bool computeBhpAtThpLimit(double alq);
        void debugShowIterationInfo(double alq);
        double getBhpWithLimit();
        void warn_(std::string msg) {parent.displayWarning_(msg);}
    };

    std::pair<std::optional<double>, bool>
    addOrSubtractAlqIncrement_(double alq, bool increase) const;

    double calcEcoGradient_(double oil_rate, double new_oil_rate,
                            double gas_rate, double new_gas_rate, bool increase) const;

    bool checkALQequal_(double alq1, double alq2) const;
    bool checkInitialALQmodified_(double alq, double initial_alq) const;

    bool checkWellRatesViolated_(std::vector<double>& potentials,
                                 const std::function<bool(double, double, const std::string &)>& callback,
                                 bool increase);

    virtual std::optional<double> computeBhpAtThpLimit_(double alq) const = 0;
    virtual void computeWellRates_(double bhp,
                                   std::vector<double>& potentials,
                                   bool debug_output = true) const = 0;

    bool computeInitialWellRates_(std::vector<double>& potentials);

    void debugCheckNegativeGradient_(double grad, double alq, double new_alq,
                                     double oil_rate, double new_oil_rate, double gas_rate,
                                     double new_gas_rate, bool increase) const;

    void debugShowAlqIncreaseDecreaseCounts_();

    void debugShowBhpAlqTable_();

    void debugShowStartIteration_(double alq, bool increase, double oil_rate);

    void debugShowTargets_();

    void displayDebugMessage_(const std::string& msg) const;
    void displayWarning_(const std::string& warning);

    std::pair<double, bool> getBhpWithLimit_(double bhp) const;
    std::pair<double, bool> getGasRateWithLimit_(const std::vector<double>& potentials) const;
    std::tuple<double,double,bool,bool>
    getInitialRatesWithLimit_(const std::vector<double>& potentials);
    std::pair<double, bool> getOilRateWithLimit_(const std::vector<double>& potentials) const;

    std::tuple<double,double,bool,bool,double>
    increaseALQtoPositiveOilRate_(double alq,
                                  double oil_rate,
                                  double gas_rate,
                                  bool oil_is_limited,
                                  bool gas_is_limited,
                                  std::vector<double>& potentials);

    std::tuple<double,double,bool,bool,double>
    increaseALQtoMinALQ_(double alq,
                         double oil_rate,
                         double gas_rate,
                         bool oil_is_limited,
                         bool gas_is_limited,
                         std::vector<double>& potentials);

    void logSuccess_(double alq,
                     const int iteration_idx);

    std::tuple<double,double,bool,bool,double>
    reduceALQtoOilTarget_(double alq,
                          double oil_rate,
                          double gas_rate,
                          bool oil_is_limited,
                          bool gas_is_limited,
                          std::vector<double>& potentials);


    std::unique_ptr<GasLiftWellState> runOptimize1_();
    std::unique_ptr<GasLiftWellState> runOptimize2_();
    std::unique_ptr<GasLiftWellState> runOptimizeLoop_(bool increase);

    void setAlqMinRate_(const GasLiftOpt::Well& well);

    std::unique_ptr<GasLiftWellState> tryIncreaseLiftGas_();
    std::unique_ptr<GasLiftWellState> tryDecreaseLiftGas_();

    void updateWellStateAlqFixedValue_(const GasLiftOpt::Well& well);

    bool useFixedAlq_(const GasLiftOpt::Well& well);

    void warnMaxIterationsExceeded_();

    DeferredLogger& deferred_logger_;
    WellState& well_state_;
    const Well& ecl_well_;
    const SummaryState& summary_state_;

    const Well::ProductionControls controls_;

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
    int num_phases_;

    std::string well_name_;

    const GasLiftOpt::Well* gl_well_;

    bool optimize_;
    bool debug_;  // extra debug output
    bool debug_limit_increase_decrease_;
    bool debug_abort_if_decrease_and_oil_is_limited_ = false;
    bool debug_abort_if_increase_and_gas_is_limited_ = false;
};

} // namespace Opm

#endif // OPM_GASLIFT_SINGLE_WELL_GENERIC_HEADER_INCLUDED
