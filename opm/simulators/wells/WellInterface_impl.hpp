/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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

#ifndef OPM_WELLINTERFACE_IMPL_HEADER_INCLUDED
#define OPM_WELLINTERFACE_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/WellInterface.hpp>
#endif

#include <opm/common/Exceptions.hpp>

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>

#include <dune/common/version.hh>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

#include <fmt/format.h>

namespace Opm
{


    template<typename TypeTag>
    WellInterface<TypeTag>::
    WellInterface(const Well& well,
                  const ParallelWellInfo<Scalar>& pw_info,
                  const int time_step,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components,
                  const int num_phases,
                  const int index_of_well,
                  const std::vector<PerforationData<Scalar>>& perf_data)
      : WellInterfaceIndices<FluidSystem,Indices>(well,
                                                  pw_info,
                                                  time_step,
                                                  param,
                                                  rate_converter,
                                                  pvtRegionIdx,
                                                  num_components,
                                                  num_phases,
                                                  index_of_well,
                                                  perf_data)
    {
        connectionRates_.resize(this->number_of_perforations_);

        if constexpr (has_solvent || has_zFraction) {
            if (well.isInjector()) {
                auto injectorType = this->well_ecl_.injectorType();
                if (injectorType == InjectorType::GAS) {
                    this->wsolvent_ = this->well_ecl_.getSolventFraction();
                }
            }
        }
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<Scalar>& /* depth_arg */,
         const Scalar gravity_arg,
         const std::vector<Scalar>& B_avg,
         const bool changed_to_open_this_step)
    {
        this->phase_usage_ = phase_usage_arg;
        this->gravity_ = gravity_arg;
        B_avg_ = B_avg;
        this->changed_to_open_this_step_ = changed_to_open_this_step;
    }




    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    wpolymer() const
    {
        if constexpr (has_polymer) {
            return this->wpolymer_();
        }

        return 0.0;
    }





    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    wfoam() const
    {
        if constexpr (has_foam) {
            return this->wfoam_();
        }

        return 0.0;
    }



    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    wsalt() const
    {
        if constexpr (has_brine) {
            return this->wsalt_();
        }

        return 0.0;
    }

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    wmicrobes() const
    {
      if constexpr (has_micp) {
          return this->wmicrobes_();
      }

      return 0.0;
    }

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    woxygen() const
    {
      if constexpr (has_micp) {
          return this->woxygen_();
      }

      return 0.0;
    }

    // The urea injection concentration is scaled down by a factor of 10, since its value
    // can be much bigger than 1 (not doing this slows the simulations). The
    // corresponding values are scaled accordingly in blackoilmicpmodules.hh when computing
    // the reactions and also when writing the output files (vtk and eclipse format, i.e.,
    // vtkblackoilmicpmodule.hh and ecloutputblackoilmodel.hh respectively).

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    wurea() const
    {
      if constexpr (has_micp) {
          return this->wurea_();
      }

      return 0.0;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellControl(const Simulator& simulator,
                      const IndividualOrGroup iog,
                      WellState<Scalar>& well_state,
                      const GroupState<Scalar>& group_state,
                      DeferredLogger& deferred_logger) /* const */
    {
        OPM_TIMEFUNCTION();
        if (stoppedOrZeroRateTarget(simulator, well_state, deferred_logger)) {
            return false;
        }

        const auto& summaryState = simulator.vanguard().summaryState();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& well = this->well_ecl_;
        auto& ws = well_state.well(this->index_of_well_);
        std::string from;
        if (well.isInjector()) {
            from = WellInjectorCMode2String(ws.injection_cmode);
        } else {
            from = WellProducerCMode2String(ws.production_cmode);
        }
        bool oscillating = std::count(this->well_control_log_.begin(), this->well_control_log_.end(), from) >= this->param_.max_number_of_well_switches_;
        const int episodeIdx = simulator.episodeIndex();
        const int iterationIdx = simulator.model().newtonMethod().numIterations();
        const int nupcol = schedule[episodeIdx].nupcol();
        if (oscillating && iterationIdx > nupcol) {
            // only output frist time
            bool output = std::count(this->well_control_log_.begin(), this->well_control_log_.end(), from) == this->param_.max_number_of_well_switches_;
            if (output) {
                std::ostringstream ss;
                ss << "    The control mode for well " << this->name()
                   << " is oscillating\n"
                   << "    We don't allow for more than "
                   << this->param_.max_number_of_well_switches_
                   << " switches. The control is kept at " << from;
                deferred_logger.info(ss.str());
                // add one more to avoid outputting the same info again
                this->well_control_log_.push_back(from);
            }
            return false;
        }
        bool changed = false;
        if (iog == IndividualOrGroup::Individual) {
            changed = this->checkIndividualConstraints(ws, summaryState, deferred_logger);
        } else if (iog == IndividualOrGroup::Group) {
            changed = this->checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        } else {
            assert(iog == IndividualOrGroup::Both);
            changed = this->checkConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        }
        Parallel::Communication cc = simulator.vanguard().grid().comm();
        // checking whether control changed
        if (changed) {
            std::string to;
            if (well.isInjector()) {
                to = WellInjectorCMode2String(ws.injection_cmode);
            } else {
                to = WellProducerCMode2String(ws.production_cmode);
            }
            std::ostringstream ss;
            ss << "    Switching control mode for well " << this->name()
               << " from " << from
               << " to " <<  to;
            if (cc.size() > 1) {
               ss << " on rank " << cc.rank();
            }
            deferred_logger.debug(ss.str());

            this->well_control_log_.push_back(from);
            updateWellStateWithTarget(simulator, group_state, well_state, deferred_logger);
            updatePrimaryVariables(simulator, well_state, deferred_logger);
        }

        return changed;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellControlAndStatusLocalIteration(const Simulator& simulator,
                                             WellState<Scalar>& well_state,
                                             GroupState<Scalar>& group_state,
                                             const Well::InjectionControls& inj_controls,
                                             const Well::ProductionControls& prod_controls,
                                             const Scalar wqTotal,
                                             DeferredLogger& deferred_logger, 
                                             const bool fixed_control, 
                                             const bool fixed_status)
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();
        const auto& schedule = simulator.vanguard().schedule();
        auto& ws = well_state.well(this->index_of_well_);
        std::string from;
        if (this->isInjector()) {
            from = WellInjectorCMode2String(ws.injection_cmode);
        } else {
            from = WellProducerCMode2String(ws.production_cmode);
        }
        const bool oscillating = std::count(this->well_control_log_.begin(), this->well_control_log_.end(), from) >= this->param_.max_number_of_well_switches_;

        if (oscillating || this->wellUnderZeroRateTarget(simulator, well_state, deferred_logger) || !(this->well_ecl_.getStatus() == WellStatus::OPEN)) {
           return false;
        }

        const Scalar sgn = this->isInjector() ? 1.0 : -1.0;
        if (!this->wellIsStopped()){
            if (wqTotal*sgn <= 0.0 && !fixed_status){
                this->stopWell();
                return true;
            } else {
                bool changed = false;
                if (!fixed_control) {
                    // Changing to group controls here may lead to inconsistencies in the group handling which in turn 
                    // may result in excessive back and forth switching. However, we currently allow this by default.
                    // The switch check_group_constraints_inner_well_iterations_ is a temporary solution.
                    
                    const bool hasGroupControl = this->isInjector() ? inj_controls.hasControl(Well::InjectorCMode::GRUP) :
                                                                      prod_controls.hasControl(Well::ProducerCMode::GRUP);
                    bool isGroupControl = ws.production_cmode == Well::ProducerCMode::GRUP || ws.injection_cmode == Well::InjectorCMode::GRUP; 
                    if (! (isGroupControl && !this->param_.check_group_constraints_inner_well_iterations_)) {
                        changed = this->checkIndividualConstraints(ws, summary_state, deferred_logger, inj_controls, prod_controls);
                    }
                    if (hasGroupControl && this->param_.check_group_constraints_inner_well_iterations_) {
                        changed = changed || this->checkGroupConstraints(well_state, group_state, schedule, summary_state,deferred_logger);
                    }
               
                    if (changed) {
                        const bool thp_controlled = this->isInjector() ? ws.injection_cmode == Well::InjectorCMode::THP :
                                                                        ws.production_cmode == Well::ProducerCMode::THP;
                        if (!thp_controlled){
                            // don't call for thp since this might trigger additional local solve
                            updateWellStateWithTarget(simulator, group_state, well_state, deferred_logger);
                        } else {
                            ws.thp = this->getTHPConstraint(summary_state);
                        }
                        updatePrimaryVariables(simulator, well_state, deferred_logger);
                    }

                    bool isNowGroupControl = ws.production_cmode == Well::ProducerCMode::GRUP || ws.injection_cmode == Well::InjectorCMode::GRUP; 
                    if (isGroupControl != isNowGroupControl) {
                        well_state.updateGlobalIsGrup();
                        int reportStepIdx = simulator.episodeIndex();
                        const Group& fieldGroup = schedule.getGroup("FIELD", reportStepIdx);
                        std::vector<Scalar> groupTargetReduction(well_state.numPhases(), 0.0);
                        WellGroupHelpers<Scalar>::updateGroupTargetReduction(fieldGroup,
                                                         schedule,
                                                         reportStepIdx,
                                                         this->isInjector(),
                                                         well_state.phaseUsage(),
                                                         *this->guideRate(),
                                                         well_state,
                                                         summary_state,
                                                         group_state,
                                                         groupTargetReduction);

                        if (isNowGroupControl) {
                            updateWellStateWithTarget(simulator, group_state, well_state, deferred_logger);
                            updatePrimaryVariables(simulator, well_state, deferred_logger);
                        }
                    }    
                }
                return changed;
            }
        } else if (!fixed_status){
            // well is stopped, check if current bhp allows reopening
            const Scalar bhp = well_state.well(this->index_of_well_).bhp;
            Scalar prod_limit = prod_controls.bhp_limit;
            Scalar inj_limit = inj_controls.bhp_limit;
            const bool has_thp = this->wellHasTHPConstraints(summary_state);
            if (has_thp){
                std::vector<Scalar> rates(this->num_components_);
                if (this->isInjector()){
                    const Scalar bhp_thp = WellBhpThpCalculator(*this).
                                                calculateBhpFromThp(well_state, rates,
                                                                    this->well_ecl_,
                                                                    summary_state,
                                                                    this->getRefDensity(),
                                                                    deferred_logger);
                    inj_limit = std::min(bhp_thp, static_cast<Scalar>(inj_controls.bhp_limit));
                } else {
                    // if the well can operate, it must at least be able to produce
                    // at the lowest bhp of the bhp-curve (explicit fractions)
                    const Scalar bhp_min = WellBhpThpCalculator(*this).
                                                calculateMinimumBhpFromThp(well_state,
                                                                           this->well_ecl_,
                                                                           summary_state,
                                                                           this->getRefDensity());
                    prod_limit = std::max(bhp_min, static_cast<Scalar>(prod_controls.bhp_limit));
                }
            }
            const Scalar bhp_diff = (this->isInjector())? inj_limit - bhp: bhp - prod_limit;
            if (bhp_diff > 0){
                this->openWell();
                well_state.well(this->index_of_well_).bhp = (this->isInjector())? inj_limit : prod_limit;
                if (has_thp) {
                    well_state.well(this->index_of_well_).thp = this->getTHPConstraint(summary_state);
                }
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTesting(const Simulator& simulator,
                const double simulation_time,
                /* const */ WellState<Scalar>& well_state,
                const GroupState<Scalar>& group_state,
                WellTestState& well_test_state,
                const PhaseUsage& phase_usage,
                GLiftEclWells& ecl_well_map,
                std::map<std::string, double>& open_times,
                DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        deferred_logger.info(" well " + this->name() + " is being tested");

        WellState<Scalar> well_state_copy = well_state;
        GroupState<Scalar> group_state_copy = group_state;

        auto& ws = well_state_copy.well(this->indexOfWell());

        updateWellStateWithTarget(simulator, group_state, well_state_copy, deferred_logger);
        calculateExplicitQuantities(simulator, well_state_copy, deferred_logger);
        updatePrimaryVariables(simulator, well_state_copy, deferred_logger);
        initPrimaryVariablesEvaluation();

        if (this->isProducer()) {
            const auto& schedule = simulator.vanguard().schedule();
            const auto report_step = simulator.episodeIndex();
            const auto& glo = schedule.glo(report_step);
            if (glo.active()) {
                gliftBeginTimeStepWellTestUpdateALQ(simulator, 
                                                    well_state_copy, 
                                                    group_state,
                                                    phase_usage, 
                                                    ecl_well_map, 
                                                    deferred_logger);
            }
        }

        WellTestState welltest_state_temp;

        bool testWell = true;
        // if a well is closed because all completions are closed, we need to check each completion
        // individually. We first open all completions, then we close one by one by calling updateWellTestState
        // untill the number of closed completions do not increase anymore.
        while (testWell) {
            const std::size_t original_number_closed_completions = welltest_state_temp.num_closed_completions();
            bool converged = solveWellForTesting(simulator, well_state_copy, group_state_copy, deferred_logger);
            if (!converged) {
                const auto msg = fmt::format("WTEST: Well {} is not solvable (physical)", this->name());
                deferred_logger.debug(msg);
                return;
            }


            updateWellOperability(simulator, well_state_copy, deferred_logger);
            if ( !this->isOperableAndSolvable() ) {
                const auto msg = fmt::format("WTEST: Well {} is not operable (physical)", this->name());
                deferred_logger.debug(msg);
                return;
            }
            std::vector<Scalar> potentials;
            try {
                computeWellPotentials(simulator, well_state_copy, potentials, deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg = fmt::format("well {}: computeWellPotentials() "
                                                    "failed during testing for re-opening: ",
                                                    this->name(), e.what());
                deferred_logger.info(msg);
                return;
            }
            const int np = well_state_copy.numPhases();
            for (int p = 0; p < np; ++p) {
                ws.well_potentials[p] = std::max(Scalar{0.0}, potentials[p]);
            }
            const bool under_zero_target = this->wellUnderZeroGroupRateTarget(simulator, well_state_copy, deferred_logger);
            this->updateWellTestState(well_state_copy.well(this->indexOfWell()),
                                     simulation_time,
                                      /*writeMessageToOPMLog=*/ false,
                                      under_zero_target,
                                      welltest_state_temp,
                                      deferred_logger);
            this->closeCompletions(welltest_state_temp);

            // Stop testing if the well is closed or shut due to all completions shut
            // Also check if number of completions has increased. If the number of closed completions do not increased
            // we stop the testing.
            // TODO: it can be tricky here, if the well is shut/closed due to other reasons
            if ( welltest_state_temp.num_closed_wells() > 0 ||
                (original_number_closed_completions == welltest_state_temp.num_closed_completions()) ) {
                    testWell = false; // this terminates the while loop
            }
        }

        // update wellTestState if the well test succeeds
        if (!welltest_state_temp.well_is_closed(this->name())) {
            well_test_state.open_well(this->name());

            std::string msg = std::string("well ") + this->name() + std::string(" is re-opened");
            deferred_logger.info(msg);

            // also reopen completions
            for (const auto& completion : this->well_ecl_.getCompletions()) {
                if (!welltest_state_temp.completion_is_closed(this->name(), completion.first))
                    well_test_state.open_completion(this->name(), completion.first);
            }
            // set the status of the well_state to open
            ws.open();
            well_state = well_state_copy;
            open_times.try_emplace(this->name(), well_test_state.lastTestTime(this->name()));
        }
    }




    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    iterateWellEquations(const Simulator& simulator,
                         const double dt,
                         WellState<Scalar>& well_state,
                         GroupState<Scalar>& group_state,
                         DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();
        const auto inj_controls = this->well_ecl_.isInjector() ? this->well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = this->well_ecl_.isProducer() ? this->well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        bool converged = false;
        try {
            // TODO: the following two functions will be refactored to be one to reduce the code duplication
            if (!this->param_.local_well_solver_control_switching_){
                converged = this->iterateWellEqWithControl(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
            } else {
                if (this->param_.use_implicit_ipr_ && this->well_ecl_.isProducer() && this->wellHasTHPConstraints(summary_state) && (this->well_ecl_.getStatus() == WellStatus::OPEN)) {
                    converged = solveWellWithTHPConstraint(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
                } else {
                    converged = this->iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
                }
            }

        } catch (NumericalProblem& e ) {
            const std::string msg = "Inner well iterations failed for well " + this->name() + " Treat the well as unconverged. ";
            deferred_logger.warning("INNER_ITERATION_FAILED", msg);
            converged = false;
        }
        return converged;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellWithTHPConstraint(const Simulator& simulator,
                               const double dt,
                               const Well::InjectionControls& inj_controls,
                               const Well::ProductionControls& prod_controls,
                               WellState<Scalar>& well_state,
                               GroupState<Scalar>& group_state,
                               DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();
        bool is_operable = true;
        bool converged = true;
        auto& ws = well_state.well(this->index_of_well_);
        // if well is stopped, check if we can reopen
        if (this->wellIsStopped()) {
            this->openWell();
            auto bhp_target = estimateOperableBhp(simulator, dt, well_state, summary_state, deferred_logger);
            if (!bhp_target.has_value()) {
                // no intersection with ipr
                const auto msg = fmt::format("estimateOperableBhp: Did not find operable BHP for well {}", this->name());
                deferred_logger.debug(msg);
                is_operable = false;
                // solve with zero rates
                solveWellWithZeroRate(simulator, dt, well_state, deferred_logger);
                this->stopWell();
            } else {
                // solve well with the estimated target bhp (or limit)
                ws.thp = this->getTHPConstraint(summary_state);
                const Scalar bhp = std::max(bhp_target.value(),
                                            static_cast<Scalar>(prod_controls.bhp_limit));
                solveWellWithBhp(simulator, dt, bhp, well_state, deferred_logger);
            }
        }
        // solve well-equation
        if (is_operable) {
            converged = this->iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
        }

        const bool isThp = ws.production_cmode == Well::ProducerCMode::THP;
        // check stability of solution under thp-control
        if (converged && !stoppedOrZeroRateTarget(simulator, well_state, deferred_logger) && isThp) {
            auto rates = well_state.well(this->index_of_well_).surface_rates;
            this->adaptRatesForVFP(rates);
            this->updateIPRImplicit(simulator, well_state, deferred_logger);
            bool is_stable = WellBhpThpCalculator(*this).isStableSolution(well_state, this->well_ecl_, rates, summary_state);
            if (!is_stable) {
                // solution converged to an unstable point!
                this->operability_status_.use_vfpexplicit = true;
                auto bhp_stable = WellBhpThpCalculator(*this).estimateStableBhp(well_state, this->well_ecl_, rates, this->getRefDensity(), summary_state);
                // if we find an intersection with a sufficiently lower bhp, re-solve equations
                const Scalar reltol = 1e-3;
                const Scalar cur_bhp = ws.bhp;
                if (bhp_stable.has_value() && cur_bhp - bhp_stable.value() > cur_bhp*reltol){
                    const auto msg = fmt::format("Well {} converged to an unstable solution, re-solving", this->name());
                    deferred_logger.debug(msg);
                    solveWellWithBhp(simulator, dt, bhp_stable.value(), well_state, deferred_logger);
                    // re-solve with hopefully good initial guess
                    ws.thp = this->getTHPConstraint(summary_state);
                    converged = this->iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
                }
            }
        }

        if (!converged) {
            // Well did not converge, switch to explicit fractions
            this->operability_status_.use_vfpexplicit = true;
            this->openWell();
            auto bhp_target = estimateOperableBhp(simulator, dt, well_state, summary_state, deferred_logger);
            if (!bhp_target.has_value()) {
                // well can't operate using explicit fractions
                is_operable = false;
                // solve with zero rate
                converged = solveWellWithZeroRate(simulator, dt, well_state, deferred_logger);
                this->stopWell();
            } else {
                // solve well with the estimated target bhp (or limit)
                const Scalar bhp = std::max(bhp_target.value(),
                                            static_cast<Scalar>(prod_controls.bhp_limit));
                solveWellWithBhp(simulator, dt, bhp, well_state, deferred_logger);
                ws.thp = this->getTHPConstraint(summary_state);
                converged = this->iterateWellEqWithSwitching(simulator, dt,
                                                             inj_controls,
                                                             prod_controls,
                                                             well_state,
                                                             group_state,
                                                             deferred_logger);
            }
        }
        // update operability
        is_operable = is_operable && !this->wellIsStopped();
        this->operability_status_.can_obtain_bhp_with_thp_limit = is_operable;
        this->operability_status_.obey_thp_limit_under_bhp_limit = is_operable;
        return converged;
    }

    template<typename TypeTag>
    std::optional<typename WellInterface<TypeTag>::Scalar>
    WellInterface<TypeTag>::
    estimateOperableBhp(const Simulator& simulator,
                        const double dt,
                        WellState<Scalar>& well_state,
                        const SummaryState& summary_state,
                        DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // Given an unconverged well or closed well, estimate an operable bhp (if any)
        // Get minimal bhp from vfp-curve
        Scalar bhp_min =  WellBhpThpCalculator(*this).calculateMinimumBhpFromThp(well_state, this->well_ecl_, summary_state, this->getRefDensity());
        // Solve
        const bool converged = solveWellWithBhp(simulator, dt, bhp_min, well_state, deferred_logger);
        if (!converged || this->wellIsStopped()) {
            return std::nullopt;
        }
        this->updateIPRImplicit(simulator, well_state, deferred_logger);
        auto rates = well_state.well(this->index_of_well_).surface_rates;
        this->adaptRatesForVFP(rates);
        return WellBhpThpCalculator(*this).estimateStableBhp(well_state, this->well_ecl_, rates, this->getRefDensity(), summary_state);
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellWithBhp(const Simulator& simulator,
                     const double dt,
                     const Scalar bhp,
                     WellState<Scalar>& well_state,
                     DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // Solve a well using single bhp-constraint (but close if not operable under this)
        auto group_state = GroupState<Scalar>(); // empty group
        auto inj_controls = Well::InjectionControls(0);
        auto prod_controls = Well::ProductionControls(0);
        auto& ws = well_state.well(this->index_of_well_);
        auto cmode_inj = ws.injection_cmode;
        auto cmode_prod = ws.production_cmode;
        if (this->isInjector()) {
            inj_controls.addControl(Well::InjectorCMode::BHP);
            inj_controls.bhp_limit = bhp;
            inj_controls.cmode = Well::InjectorCMode::BHP;
            ws.injection_cmode = Well::InjectorCMode::BHP;
        } else {
            prod_controls.addControl(Well::ProducerCMode::BHP);
            prod_controls.bhp_limit = bhp;
            prod_controls.cmode = Well::ProducerCMode::BHP;
            ws.production_cmode = Well::ProducerCMode::BHP;
        }
        // update well-state
        ws.bhp = bhp;
        // solve
        const bool converged =  this->iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger, /*fixed_control*/true);
        ws.injection_cmode = cmode_inj;
        ws.production_cmode = cmode_prod;
        return converged;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellWithZeroRate(const Simulator& simulator,
                          const double dt,
                          WellState<Scalar>& well_state,
                          DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // Solve a well as stopped
        const auto well_status_orig = this->wellStatus_;
        this->stopWell();

        auto group_state = GroupState<Scalar>(); // empty group
        auto inj_controls = Well::InjectionControls(0);
        auto prod_controls = Well::ProductionControls(0);
        const bool converged =  this->iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger, /*fixed_control*/true, /*fixed_status*/ true);
        this->wellStatus_ = well_status_orig;
        return converged;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellForTesting(const Simulator& simulator,
                        WellState<Scalar>& well_state,
                        GroupState<Scalar>& group_state,
                        DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // keep a copy of the original well state
        const WellState<Scalar> well_state0 = well_state;
        const double dt = simulator.timeStepSize();
        const auto& summary_state = simulator.vanguard().summaryState();
        const bool has_thp_limit = this->wellHasTHPConstraints(summary_state);
        bool converged;
        if (has_thp_limit) {
            well_state.well(this->indexOfWell()).production_cmode = Well::ProducerCMode::THP;
        }
        else {
            well_state.well(this->indexOfWell()).production_cmode = Well::ProducerCMode::BHP;
        }
        converged = iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("WellTest: Well equation for well " + this->name() +  " converged");
            return true;
        }
        const int max_iter = this->param_.max_welleq_iter_;
        deferred_logger.debug("WellTest: Well equation for well " + this->name() + " failed converging in "
                              + std::to_string(max_iter) + " iterations");
        well_state = well_state0;
        return false;
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellEquation(const Simulator& simulator,
                      WellState<Scalar>& well_state,
                      GroupState<Scalar>& group_state,
                      DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        if (!this->isOperableAndSolvable() && !this->wellIsStopped())
            return;

        // keep a copy of the original well state
        const WellState<Scalar> well_state0 = well_state;
        const double dt = simulator.timeStepSize();
        bool converged = iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger);

        // Newly opened wells with THP control sometimes struggles to
        // converge due to bad initial guess. Or due to the simple fact
        // that the well needs to change to another control.
        // We therefore try to solve the well with BHP control to get
        // an better initial guess.
        // If the well is supposed to operate under THP control
        // "updateWellControl" will switch it back to THP later.
        if (!converged) {
            auto& ws = well_state.well(this->indexOfWell());
            bool thp_control = false;
            if (this->well_ecl_.isInjector()) {
                thp_control = ws.injection_cmode == Well::InjectorCMode::THP;
                if (thp_control) {
                    ws.injection_cmode = Well::InjectorCMode::BHP;
                    this->well_control_log_.push_back(WellInjectorCMode2String(Well::InjectorCMode::THP));
                }
            } else {
                thp_control = ws.production_cmode == Well::ProducerCMode::THP;
                if (thp_control) {
                    ws.production_cmode = Well::ProducerCMode::BHP;
                    this->well_control_log_.push_back(WellProducerCMode2String(Well::ProducerCMode::THP));
                }
            }
            if (thp_control) {
                const std::string msg = std::string("The newly opened well ") + this->name()
                    + std::string(" with THP control did not converge during inner iterations, we try again with bhp control");
                deferred_logger.debug(msg);
                converged = this->iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger);
            }
        }

        if (!converged) {
            const int max_iter = this->param_.max_welleq_iter_;
            deferred_logger.debug("Compute initial well solution for well " + this->name() + ". Failed to converge in "
                                  + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    assembleWellEq(const Simulator& simulator,
                   const double dt,
                   WellState<Scalar>& well_state,
                   GroupState<Scalar>& group_state,
                   DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        prepareWellBeforeAssembling(simulator, dt, well_state, group_state, deferred_logger);
        assembleWellEqWithoutIteration(simulator, dt, well_state, group_state, deferred_logger);
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& simulator,
                                   const double dt,
                                   WellState<Scalar>& well_state,
                                   const GroupState<Scalar>& group_state,
                                   DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();
        const auto inj_controls = this->well_ecl_.isInjector() ? this->well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = this->well_ecl_.isProducer() ? this->well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        // TODO: the reason to have inj_controls and prod_controls in the arguments, is that we want to change the control used for the well functions
        // TODO: maybe we can use std::optional or pointers to simplify here
        assembleWellEqWithoutIteration(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    prepareWellBeforeAssembling(const Simulator& simulator,
                                const double dt,
                                WellState<Scalar>& well_state,
                                GroupState<Scalar>& group_state,
                                DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const bool old_well_operable = this->operability_status_.isOperableAndSolvable();

        if (this->param_.check_well_operability_iter_)
            checkWellOperability(simulator, well_state, deferred_logger);

        // only use inner well iterations for the first newton iterations.
        const int iteration_idx = simulator.model().newtonMethod().numIterations();
        if (iteration_idx < this->param_.max_niter_inner_well_iter_ || this->well_ecl_.isMultiSegment()) {
            const auto& ws = well_state.well(this->indexOfWell());
            const auto pmode_orig = ws.production_cmode;
            const auto imode_orig = ws.injection_cmode;

            const bool nonzero_rate_original =
                std::any_of(ws.surface_rates.begin(),
                            ws.surface_rates.begin() + well_state.numPhases(),
                            [](Scalar rate) { return rate != Scalar(0.0); });

            this->operability_status_.solvable = true;
            bool converged = this->iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger);

            if (converged) {
                const bool zero_target = this->wellUnderZeroRateTarget(simulator, well_state, deferred_logger);
                if (this->wellIsStopped() && !zero_target && nonzero_rate_original) {
                    // Well had non-zero rate, but was stopped during local well-solve. We re-open the well 
                    // for the next global iteration, but if the zero rate persists, it will be stopped.
                    // This logic is introduced to prevent/ameliorate stopped/revived oscillations  
                    this->operability_status_.resetOperability();
                    this->openWell();
                    deferred_logger.debug("    " + this->name() + " is re-opened after being stopped during local solve");
                }
                // Add debug info for switched controls
                if (ws.production_cmode != pmode_orig || ws.injection_cmode != imode_orig) {
                    std::string from,to;
                    if (this->isInjector()) {
                        from = WellInjectorCMode2String(imode_orig);
                        to = WellInjectorCMode2String(ws.injection_cmode);
                    } else {
                        from = WellProducerCMode2String(pmode_orig);
                        to = WellProducerCMode2String(ws.production_cmode);
                    }
                    deferred_logger.debug("    " + this->name() + " switched from " + from + " to " + to + " during local solve");
                    }

            } else {
                // unsolvable wells are treated as not operable and will not be solved for in this iteration.
                if (this->param_.shut_unsolvable_wells_)
                    this->operability_status_.solvable = false;
            }
        }
        if (this->operability_status_.has_negative_potentials) {
            auto well_state_copy = well_state;
            std::vector<Scalar> potentials;
            try {
                computeWellPotentials(simulator, well_state_copy, potentials, deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg = fmt::format("well {}: computeWellPotentials() failed "
                                                    "during attempt to recompute potentials for well: ",
                                                    this->name(), e.what());
                deferred_logger.info(msg);
                this->operability_status_.has_negative_potentials = true;
            }
            auto& ws = well_state.well(this->indexOfWell());
            const int np = well_state.numPhases();
            for (int p = 0; p < np; ++p) {
                ws.well_potentials[p] = std::max(Scalar{0.0}, potentials[p]);
            }
        }
        this->changed_to_open_this_step_ = false;
        const bool well_operable = this->operability_status_.isOperableAndSolvable();

        if (!well_operable && old_well_operable) {
            deferred_logger.info(" well " + this->name() + " gets STOPPED during iteration ");
            this->stopWell();
            changed_to_stopped_this_step_ = true;
        } else if (well_operable && !old_well_operable) {
            deferred_logger.info(" well " + this->name() + " gets REVIVED during iteration ");
            this->openWell();
            changed_to_stopped_this_step_ = false;
            this->changed_to_open_this_step_ = true;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::addCellRates(RateVector& rates, int cellIdx) const
    {
        if(!this->isOperableAndSolvable() && !this->wellIsStopped())
            return;

        for (int perfIdx = 0; perfIdx < this->number_of_perforations_; ++perfIdx) {
            if (this->cells()[perfIdx] == cellIdx) {
                for (int i = 0; i < RateVector::dimension; ++i) {
                    rates[i] += connectionRates_[perfIdx][i];
                }
            }
        }
    }

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const
    {
        for (int perfIdx = 0; perfIdx < this->number_of_perforations_; ++perfIdx) {
            if (this->cells()[perfIdx] == cellIdx) {
                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                return connectionRates_[perfIdx][activeCompIdx].value();
            }
        }
        // this is not thread safe
        OPM_THROW(std::invalid_argument, "The well with name " + this->name()
                  + " does not perforate cell " + std::to_string(cellIdx));
        return 0.0;
    }




    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkWellOperability(const Simulator& simulator,
                         const WellState<Scalar>& well_state,
                         DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        if (!this->param_.check_well_operability_) {
            return;
        }

        if (this->wellIsStopped() && !changed_to_stopped_this_step_) {
            return;
        }

        updateWellOperability(simulator, well_state, deferred_logger);
        if (!this->operability_status_.isOperableAndSolvable()) {
            this->operability_status_.use_vfpexplicit = true;
            deferred_logger.debug("EXPLICIT_LOOKUP_VFP",
                                "well not operable, trying with explicit vfp lookup: " + this->name());
            updateWellOperability(simulator, well_state, deferred_logger);
        }
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    gliftBeginTimeStepWellTestIterateWellEquations(const Simulator& simulator,
                                                   const double dt,
                                                   WellState<Scalar>& well_state,
                                                   GroupState<Scalar>& group_state,
                                                   DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const auto& well_name = this->name();
        assert(this->wellHasTHPConstraints(simulator.vanguard().summaryState()));
        const auto& schedule = simulator.vanguard().schedule();
        auto report_step_idx = simulator.episodeIndex();
        const auto& glo = schedule.glo(report_step_idx);
        if(glo.active() && glo.has_well(well_name)) {
            const auto increment = glo.gaslift_increment();
            auto alq = well_state.getALQ(well_name);
            bool converged;
            while (alq > 0) {
                well_state.setALQ(well_name, alq);
                if ((converged =
                      iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger)))
                {
                    return converged;
                }
                alq -= increment;
            }
            return false;
        }
        else {
            return iterateWellEquations(simulator, dt, well_state, group_state, deferred_logger);
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    gliftBeginTimeStepWellTestUpdateALQ(const Simulator& simulator,
                                        WellState<Scalar>& well_state,
                                        const GroupState<Scalar>& group_state,
                                        const PhaseUsage& phase_usage,
                                        GLiftEclWells& ecl_well_map,
                                        DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();
        const auto& well_name = this->name();
        if (!this->wellHasTHPConstraints(summary_state)) {
            const std::string msg = fmt::format("GLIFT WTEST: Well {} does not have THP constraints", well_name);
            deferred_logger.info(msg);
            return;
        }
        const auto& schedule = simulator.vanguard().schedule();
        const auto report_step_idx = simulator.episodeIndex();
        const auto& glo = schedule.glo(report_step_idx);
        if (!glo.has_well(well_name)) {
            const std::string msg = fmt::format(
                "GLIFT WTEST: Well {} : Gas lift not activated: "
                "WLIFTOPT is probably missing. Skipping.", well_name);
            deferred_logger.info(msg);
            return;
        }
        const auto& gl_well = glo.well(well_name);

        // Use gas lift optimization to get ALQ for well test
        std::unique_ptr<GasLiftSingleWell> glift =
            initializeGliftWellTest_<GasLiftSingleWell>(simulator,
                                                        well_state,
                                                        group_state,
                                                        phase_usage,
                                                        ecl_well_map,
                                                        deferred_logger);
        auto [wtest_alq, success] = glift->wellTestALQ();
        std::string msg;
        const auto& unit_system = schedule.getUnits();
        if (success) {
            well_state.setALQ(well_name, wtest_alq);
            msg = fmt::format(
                "GLIFT WTEST: Well {} : Setting ALQ to optimized value = {}",
                well_name, unit_system.from_si(UnitSystem::measure::gas_surface_rate, wtest_alq));
        }
        else {
            if (!gl_well.use_glo()) {
                msg = fmt::format(
                    "GLIFT WTEST: Well {} : Gas lift optimization deactivated. Setting ALQ to WLIFTOPT item 3 = {}",
                    well_name, 
                    unit_system.from_si(UnitSystem::measure::gas_surface_rate, well_state.getALQ(well_name)));
                
            }
            else {
                msg = fmt::format(
                    "GLIFT WTEST: Well {} : Gas lift optimization failed, no ALQ set.",
                    well_name);
            }
        }
        deferred_logger.info(msg);
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellOperability(const Simulator& simulator,
                          const WellState<Scalar>& well_state,
                          DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        if (this->param_.local_well_solver_control_switching_) {
            const bool success = updateWellOperabilityFromWellEq(simulator, well_state, deferred_logger);
            if (success) {
                return;
            } else {
                deferred_logger.debug("Operability check using well equations did not converge for well "
                                      + this->name() + ", reverting to classical approach." );
            }
        }
        this->operability_status_.resetOperability();

        bool thp_controlled = this->isInjector() ? well_state.well(this->index_of_well_).injection_cmode == Well::InjectorCMode::THP:
                                              well_state.well(this->index_of_well_).production_cmode == Well::ProducerCMode::THP;
        bool bhp_controlled = this->isInjector() ? well_state.well(this->index_of_well_).injection_cmode == Well::InjectorCMode::BHP:
                                              well_state.well(this->index_of_well_).production_cmode == Well::ProducerCMode::BHP;

        // Operability checking is not free
        // Only check wells under BHP and THP control
        bool check_thp = thp_controlled || this->operability_status_.thp_limit_violated_but_not_switched;
        if (check_thp || bhp_controlled) {
            updateIPR(simulator, deferred_logger);
            checkOperabilityUnderBHPLimit(well_state, simulator, deferred_logger);
        }
        // we do some extra checking for wells under THP control.
        if (check_thp) {
            checkOperabilityUnderTHPLimit(simulator, well_state, deferred_logger);
        }
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellOperabilityFromWellEq(const Simulator& simulator,
                                    const WellState<Scalar>& well_state,
                                    DeferredLogger& deferred_logger)
    {
        OPM_TIMEFUNCTION();
        // only makes sense if we're using this parameter is true
        assert(this->param_.local_well_solver_control_switching_);
        this->operability_status_.resetOperability();
        WellState<Scalar> well_state_copy = well_state;
        auto group_state = simulator.problem().wellModel().groupState();
        const double dt = simulator.timeStepSize();
        // equations should be converged at this stage, so only one it is needed
        bool converged = iterateWellEquations(simulator, dt, well_state_copy, group_state, deferred_logger);
        return converged;
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateWithTarget(const Simulator& simulator,
                              const GroupState<Scalar>& group_state,
                              WellState<Scalar>& well_state,
                              DeferredLogger& deferred_logger) const
    {
        OPM_TIMEFUNCTION();
        // only bhp and wellRates are used to initilize the primaryvariables for standard wells
        const auto& well = this->well_ecl_;
        const int well_index = this->index_of_well_;
        auto& ws = well_state.well(well_index);
        const auto& pu = this->phaseUsage();
        const int np = well_state.numPhases();
        const auto& summaryState = simulator.vanguard().summaryState();
        const auto& schedule = simulator.vanguard().schedule();

        if (this->wellIsStopped()) {
            for (int p = 0; p<np; ++p) {
                ws.surface_rates[p] = 0;
            }
            ws.thp = 0;
            return;
        }

        if (this->isInjector() )
        {
            const auto& controls = well.injectionControls(summaryState);

            InjectorType injectorType = controls.injector_type;
            int phasePos;
            switch (injectorType) {
            case InjectorType::WATER:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Aqua];
                break;
            }
            case InjectorType::OIL:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Liquid];
                break;
            }
            case InjectorType::GAS:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Vapour];
                break;
            }
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Expected WATER, OIL or GAS as type for injectors "  + this->name(), deferred_logger );
            }

            const auto current = ws.injection_cmode;

            switch (current) {
            case Well::InjectorCMode::RATE:
            {
                ws.surface_rates[phasePos] = (1.0 - this->rsRvInj()) * controls.surface_rate;
                if(this->rsRvInj() > 0) {
                    if (injectorType == InjectorType::OIL && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        ws.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = controls.surface_rate * this->rsRvInj();
                    } else if (injectorType == InjectorType::GAS && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = controls.surface_rate * this->rsRvInj();
                    } else {
                        OPM_DEFLOG_THROW(std::runtime_error, "Expected OIL or GAS as type for injectors when RS/RV (item 10) is non-zero "  + this->name(), deferred_logger );
                    }
                }
                break;
            }

            case Well::InjectorCMode::RESV:
            {
                std::vector<Scalar> convert_coeff(this->number_of_phases_, 1.0);
                this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);
                const Scalar coeff = convert_coeff[phasePos];
                ws.surface_rates[phasePos] = controls.reservoir_rate/coeff;
                break;
            }

            case Well::InjectorCMode::THP:
            {
                auto rates = ws.surface_rates;
                Scalar bhp = WellBhpThpCalculator(*this).calculateBhpFromThp(well_state,
                                                                             rates,
                                                                             well,
                                                                             summaryState,
                                                                             this->getRefDensity(),
                                                                             deferred_logger);
                ws.bhp = bhp;
                ws.thp = this->getTHPConstraint(summaryState);

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                Scalar total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0)
                    ws.surface_rates = ws.well_potentials;

                break;
            }
            case Well::InjectorCMode::BHP:
            {
                ws.bhp = controls.bhp_limit;
                Scalar total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate += ws.surface_rates[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0)
                    ws.surface_rates = ws.well_potentials;

                break;
            }
            case Well::InjectorCMode::GRUP:
            {
                assert(well.isAvailableForGroupControl());
                const auto& group = schedule.getGroup(well.groupName(), this->currentStep());
                const Scalar efficiencyFactor = well.getEfficiencyFactor() *
                                                well_state[well.name()].efficiency_scaling_factor;
                std::optional<Scalar> target =
                        this->getGroupInjectionTargetRate(group,
                                                          well_state,
                                                          group_state,
                                                          schedule,
                                                          summaryState,
                                                          injectorType,
                                                          efficiencyFactor,
                                                          deferred_logger);
                if (target)
                    ws.surface_rates[phasePos] = *target;
                break;
            }
            case Well::InjectorCMode::CMODE_UNDEFINED:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + this->name(), deferred_logger );
            }

            }
            // for wells with zero injection rate, if we assign exactly zero rate,
            // we will have to assume some trivial composition in the wellbore.
            // here, we use some small value (about 0.01 m^3/day ~= 1.e-7) to initialize
            // the zero rate target, then we can use to retain the composition information
            // within the wellbore from the previous result, and hopefully it is a good
            // initial guess for the zero rate target.
            ws.surface_rates[phasePos] = std::max(Scalar{1.e-7}, ws.surface_rates[phasePos]);

            if (ws.bhp == 0.) {
                ws.bhp = controls.bhp_limit;
            }
        }
        //Producer
        else
        {
            const auto current = ws.production_cmode;
            const auto& controls = well.productionControls(summaryState);
            switch (current) {
            case Well::ProducerCMode::ORAT:
            {
                Scalar current_rate = -ws.surface_rates[ pu.phase_pos[Oil] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.oil_rate/current_rate;
                    }
                } else {
                    const std::vector<Scalar> fractions = initialWellRateFractions(simulator, well_state);
                    double control_fraction = fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.oil_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::WRAT:
            {
                Scalar current_rate = -ws.surface_rates[ pu.phase_pos[Water] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.water_rate/current_rate;
                    }
                } else {
                    const std::vector<Scalar> fractions = initialWellRateFractions(simulator, well_state);
                    const Scalar control_fraction = fractions[pu.phase_pos[Water]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.water_rate / control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRAT:
            {
                Scalar current_rate = -ws.surface_rates[pu.phase_pos[Gas] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.gas_rate/current_rate;
                    }
                } else {
                    const std::vector<Scalar > fractions = initialWellRateFractions(simulator, well_state);
                    const Scalar control_fraction = fractions[pu.phase_pos[Gas]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.gas_rate / control_fraction;
                        }
                    }
                }

                break;

            }
            case Well::ProducerCMode::LRAT:
            {
                Scalar current_rate = - ws.surface_rates[ pu.phase_pos[Water] ]
                                      - ws.surface_rates[ pu.phase_pos[Oil] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.liquid_rate/current_rate;
                    }
                } else {
                    const std::vector<Scalar> fractions = initialWellRateFractions(simulator, well_state);
                    const Scalar control_fraction = fractions[pu.phase_pos[Water]] + fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.liquid_rate / control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::CRAT:
            {
                OPM_DEFLOG_THROW(std::runtime_error,
                                 fmt::format("CRAT control not supported, well {}", this->name()),
                                 deferred_logger);
            }
            case Well::ProducerCMode::RESV:
            {
                std::vector<Scalar> convert_coeff(this->number_of_phases_, 1.0);
                this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, ws.surface_rates, convert_coeff);
                Scalar total_res_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_res_rate -= ws.surface_rates[p] * convert_coeff[p];
                }
                if (controls.prediction_mode) {
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] *= controls.resv_rate/total_res_rate;
                        }
                    } else {
                        const std::vector<Scalar> fractions = initialWellRateFractions(simulator, well_state);
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.resv_rate / convert_coeff[p];
                        }
                    }
                } else {
                    std::vector<Scalar> hrates(this->number_of_phases_,0.);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        hrates[pu.phase_pos[Water]] = controls.water_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        hrates[pu.phase_pos[Oil]] = controls.oil_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        hrates[pu.phase_pos[Gas]] = controls.gas_rate;
                    }
                    std::vector<Scalar> hrates_resv(this->number_of_phases_,0.);
                    this->rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, this->pvtRegionIdx_, hrates, hrates_resv);
                    Scalar target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] *= target/total_res_rate;
                        }
                    } else {
                        const std::vector<Scalar> fractions = initialWellRateFractions(simulator, well_state);
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * target / convert_coeff[p];
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::BHP:
            {
                ws.bhp = controls.bhp_limit;
                Scalar total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate -= ws.surface_rates[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] = -ws.well_potentials[p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::THP:
            {
                const bool update_success = updateWellStateWithTHPTargetProd(simulator, well_state, deferred_logger);

                if (!update_success) {
                    // the following is the original way of initializing well state with THP constraint
                    // keeping it for robust reason in case that it fails to get a bhp value with THP constraint
                    // more sophisticated design might be needed in the future
                    auto rates = ws.surface_rates;
                    this->adaptRatesForVFP(rates);
                    const Scalar bhp = WellBhpThpCalculator(*this).calculateBhpFromThp(
                        well_state, rates, well, summaryState, this->getRefDensity(), deferred_logger);
                    ws.bhp = bhp;
                    ws.thp = this->getTHPConstraint(summaryState);
                    // if the total rates are negative or zero
                    // we try to provide a better initial well rate
                    // using the well potentials
                    const Scalar total_rate = -std::accumulate(rates.begin(), rates.end(), 0.0);
                    if (total_rate <= 0.0) {
                        for (int p = 0; p < this->number_of_phases_; ++p) {
                            ws.surface_rates[p] = -ws.well_potentials[p];
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRUP:
            {
                assert(well.isAvailableForGroupControl());
                const auto& group = schedule.getGroup(well.groupName(), this->currentStep());
                const Scalar efficiencyFactor = well.getEfficiencyFactor() *
                                                well_state[well.name()].efficiency_scaling_factor;
                Scalar scale = this->getGroupProductionTargetRate(group,
                                                                  well_state,
                                                                  group_state,
                                                                  schedule,
                                                                  summaryState,
                                                                  efficiencyFactor,
                                                                  deferred_logger);

                // we don't want to scale with zero and get zero rates.
                if (scale > 0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= scale;
                    }
                    ws.trivial_target = false;
                } else {
                    ws.trivial_target = true;
                }
                break;
            }
            case Well::ProducerCMode::CMODE_UNDEFINED:
            case Well::ProducerCMode::NONE:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + this->name() , deferred_logger);
                break;
            }
            } // end of switch

            if (ws.bhp == 0.) {
                ws.bhp = controls.bhp_limit;
            }
        }
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    wellUnderZeroRateTarget(const Simulator& simulator,
                            const WellState<Scalar>& well_state,
                            DeferredLogger& deferred_logger) const
    {
        OPM_TIMEFUNCTION();
        // Check if well is under zero rate control, either directly or from group
        const bool isGroupControlled = this->wellUnderGroupControl(well_state.well(this->index_of_well_));
        if (!isGroupControlled) {
            // well is not under group control, check "individual" version
            const auto& summaryState = simulator.vanguard().summaryState();
            return this->wellUnderZeroRateTargetIndividual(summaryState, well_state);
        } else {
            return this->wellUnderZeroGroupRateTarget(simulator, well_state, deferred_logger, isGroupControlled);
        }
    }

    template <typename TypeTag>
    bool
    WellInterface<TypeTag>::wellUnderZeroGroupRateTarget(const Simulator& simulator,
                                                         const WellState<Scalar>& well_state,
                                                         DeferredLogger& deferred_logger,
                                                         const std::optional<bool> group_control) const
    {
        // Check if well is under zero rate target from group
        const bool isGroupControlled = group_control.value_or(this->wellUnderGroupControl(well_state.well(this->index_of_well_)));
        if (isGroupControlled) {
            const auto& summaryState = simulator.vanguard().summaryState();
            const auto& group_state = simulator.problem().wellModel().groupState();
            const auto& schedule = simulator.vanguard().schedule();
            return this->zeroGroupRateTarget(summaryState, schedule, well_state, group_state, deferred_logger);
        }
        return false;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    stoppedOrZeroRateTarget(const Simulator& simulator,
                            const WellState<Scalar>& well_state,
                            DeferredLogger& deferred_logger) const
    {
        // Check if well is stopped or under zero rate control, either
        // directly or from group.
        return this->wellIsStopped()
            || this->wellUnderZeroRateTarget(simulator, well_state, deferred_logger);
    }

    template<typename TypeTag>
    std::vector<typename WellInterface<TypeTag>::Scalar>
    WellInterface<TypeTag>::
    initialWellRateFractions(const Simulator& simulator,
                             const WellState<Scalar>& well_state) const
    {
        OPM_TIMEFUNCTION();
        const int np = this->number_of_phases_;
        std::vector<Scalar> scaling_factor(np);
        const auto& ws = well_state.well(this->index_of_well_);

        Scalar total_potentials = 0.0;
        for (int p = 0; p<np; ++p) {
            total_potentials += ws.well_potentials[p];
        }
        if (total_potentials > 0) {
            for (int p = 0; p<np; ++p) {
                scaling_factor[p] = ws.well_potentials[p] / total_potentials;
            }
            return scaling_factor;
        }
        // if we don't have any potentials we weight it using the mobilites
        // We only need approximation so we don't bother with the vapporized oil and dissolved gas
        Scalar total_tw = 0;
        const int nperf = this->number_of_perforations_;
        for (int perf = 0; perf < nperf; ++perf) {
            total_tw += this->well_index_[perf];
        }
        for (int perf = 0; perf < nperf; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            const Scalar well_tw_fraction = this->well_index_[perf] / total_tw;
            Scalar total_mobility = 0.0;
            for (int p = 0; p < np; ++p) {
                int modelPhaseIdx = this->flowPhaseToModelPhaseIdx(p);
                total_mobility += fs.invB(modelPhaseIdx).value() * intQuants.mobility(modelPhaseIdx).value();
            }
            for (int p = 0; p < np; ++p) {
                int modelPhaseIdx = this->flowPhaseToModelPhaseIdx(p);
                scaling_factor[p] += well_tw_fraction * fs.invB(modelPhaseIdx).value() * intQuants.mobility(modelPhaseIdx).value() / total_mobility;
            }
        }
        return scaling_factor;
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateRates(const Simulator& simulator,
                         WellState<Scalar>& well_state,
                         DeferredLogger& deferred_logger) const
    {
        OPM_TIMEFUNCTION();
        // Check if the rates of this well only are single-phase, do nothing
        // if more than one nonzero rate.
        auto& ws = well_state.well(this->index_of_well_);
        int nonzero_rate_index = -1;
        const Scalar floating_point_error_epsilon = 1e-14;
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (std::abs(ws.surface_rates[p]) > floating_point_error_epsilon) {
                if (nonzero_rate_index == -1) {
                    nonzero_rate_index = p;
                } else {
                    // More than one nonzero rate.
                    return;
                }
            }
        }

        // Calculate the rates that follow from the current primary variables.
        std::vector<Scalar> well_q_s = computeCurrentWellRates(simulator, deferred_logger);

        if (nonzero_rate_index == -1) {
            // No nonzero rates.
            // Use the computed rate directly
            for (int p = 0; p < this->number_of_phases_; ++p) {
               ws.surface_rates[p] = well_q_s[this->flowPhaseToModelCompIdx(p)];
            }
            return;
        }

        // Set the currently-zero phase flows to be nonzero in proportion to well_q_s.
        const Scalar initial_nonzero_rate = ws.surface_rates[nonzero_rate_index];
        const int comp_idx_nz = this->flowPhaseToModelCompIdx(nonzero_rate_index);
        if (std::abs(well_q_s[comp_idx_nz]) > floating_point_error_epsilon) {
            for (int p = 0; p < this->number_of_phases_; ++p) {
                if (p != nonzero_rate_index) {
                    const int comp_idx = this->flowPhaseToModelCompIdx(p);
                    Scalar& rate = ws.surface_rates[p];
                    rate = (initial_nonzero_rate / well_q_s[comp_idx_nz]) * (well_q_s[comp_idx]);
                }
            }
        }
    }

    template <typename TypeTag>
    std::vector<typename WellInterface<TypeTag>::Scalar>
    WellInterface<TypeTag>::
    wellIndex(const int                      perf,
              const IntensiveQuantities&     intQuants,
              const Scalar                   trans_mult,
              const SingleWellState<Scalar>& ws) const
    {
        OPM_TIMEFUNCTION_LOCAL();
        // Add a Forchheimer term to the gas phase CTF if the run uses
        // either of the WDFAC or the WDFACCOR keywords.
        if (static_cast<std::size_t>(perf) >= this->well_cells_.size()) {
            OPM_THROW(std::invalid_argument,"The perforation index exceeds the size of the local containers - possibly wellIndex was called with a global instead of a local perforation index!");
        }
        auto wi = std::vector<Scalar>
            (this->num_components_, this->well_index_[perf] * trans_mult);

        if constexpr (! Indices::gasEnabled) {
            return wi;
        }

        const auto& wdfac = this->well_ecl_.getWDFAC();

        if (! wdfac.useDFactor() || (this->well_index_[perf] == 0.0)) {
            return wi;
        }

        const Scalar d = this->computeConnectionDFactor(perf, intQuants, ws);
        if (d < 1.0e-15) {
            return wi;
        }

        // Solve quadratic equations for connection rates satisfying the ipr and the flow-dependent skin.
        // If more than one solution, pick the one corresponding to lowest absolute rate (smallest skin).
        const auto& connection = this->well_ecl_.getConnections()[ws.perf_data.ecl_index[perf]];
        const Scalar Kh = connection.Kh();
        const Scalar scaling = 3.141592653589 * Kh * connection.wpimult();
        const unsigned gas_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

        const Scalar connection_pressure = ws.perf_data.pressure[perf];
        const Scalar cell_pressure = getValue(intQuants.fluidState().pressure(FluidSystem::gasPhaseIdx));
        const Scalar drawdown = cell_pressure - connection_pressure;
        const Scalar invB = getValue(intQuants.fluidState().invB(FluidSystem::gasPhaseIdx));
        const Scalar mob_g = getValue(intQuants.mobility(FluidSystem::gasPhaseIdx)) * invB;
        const Scalar a = d;
        const Scalar b = 2*scaling/wi[gas_comp_idx];
        const Scalar c = -2*scaling*mob_g*drawdown;

        Scalar consistent_Q = -1.0e20;
        // Find and check negative solutions (a --> -a)
        const Scalar r2n = b*b + 4*a*c;
        if (r2n >= 0) {
            const Scalar rn = std::sqrt(r2n);
            const Scalar xn1 = (b-rn)*0.5/a;
            if (xn1 <= 0) {
                consistent_Q = xn1;
            }
            const Scalar xn2 = (b+rn)*0.5/a;
            if (xn2 <= 0 && xn2 > consistent_Q) {
                consistent_Q = xn2;
            }
        }
        // Find and check positive solutions
        consistent_Q *= -1;
        const Scalar r2p = b*b - 4*a*c;
        if (r2p >= 0) {
            const Scalar rp = std::sqrt(r2p);
            const Scalar xp1 = (rp-b)*0.5/a;
            if (xp1 > 0 && xp1 < consistent_Q) {
                consistent_Q = xp1;
            }
            const Scalar xp2 = -(rp+b)*0.5/a;
            if (xp2 > 0 && xp2 < consistent_Q) {
                consistent_Q = xp2;
            }
        }
        wi[gas_comp_idx]  = 1.0/(1.0/(trans_mult * this->well_index_[perf]) + (consistent_Q/2 * d / scaling));

        return wi;
    }

    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateConnectionDFactor(const Simulator& simulator,
                            SingleWellState<Scalar>& ws) const
    {
        if (! this->well_ecl_.getWDFAC().useDFactor()) {
            return;
        }

        auto& d_factor = ws.perf_data.connection_d_factor;

        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);

            d_factor[perf] = this->computeConnectionDFactor(perf, intQuants, ws);
        }
    }

    template <typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::
    computeConnectionDFactor(const int                      perf,
                             const IntensiveQuantities&     intQuants,
                             const SingleWellState<Scalar>& ws) const
    {
        auto rhoGS = [regIdx = this->pvtRegionIdx()]() {
            return FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regIdx);
        };

        // Viscosity is evaluated at connection pressure.
        auto gas_visc = [connection_pressure = ws.perf_data.pressure[perf],
                         temperature = ws.temperature,
                         regIdx = this->pvtRegionIdx(), &intQuants]()
        {
            const auto rv = getValue(intQuants.fluidState().Rv());

            const auto& gasPvt = FluidSystem::gasPvt();

            // Note that rv here is from grid block with typically
            //    p_block > connection_pressure
            // so we may very well have rv > rv_sat
            const Scalar rv_sat = gasPvt.saturatedOilVaporizationFactor
                (regIdx, temperature, connection_pressure);

            if (! (rv < rv_sat)) {
                return gasPvt.saturatedViscosity(regIdx, temperature,
                                                 connection_pressure);
            }

            return gasPvt.viscosity(regIdx, temperature, connection_pressure,
                                    rv, getValue(intQuants.fluidState().Rvw()));
        };

        const auto& connection = this->well_ecl_.getConnections()
            [ws.perf_data.ecl_index[perf]];

        return this->well_ecl_.getWDFAC().getDFactor(rhoGS, gas_visc, connection);
    }


    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateConnectionTransmissibilityFactor(const Simulator& simulator,
                                           SingleWellState<Scalar>& ws) const
    {
        auto connCF = [&connIx = std::as_const(ws.perf_data.ecl_index),
                       &conns = this->well_ecl_.getConnections()]
            (const int perf)
        {
            return conns[connIx[perf]].CF();
        };

        auto& tmult = ws.perf_data.connection_compaction_tmult;
        auto& ctf   = ws.perf_data.connection_transmissibility_factor;

        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];

            const auto& intQuants = simulator.model()
                .intensiveQuantities(cell_idx, /*timeIdx=*/ 0);

            tmult[perf] = simulator.problem()
                .template wellTransMultiplier<double>(intQuants, cell_idx);

            ctf[perf] = connCF(perf) * tmult[perf];
        }
    }


    template<typename TypeTag>
    typename WellInterface<TypeTag>::Eval
    WellInterface<TypeTag>::getPerfCellPressure(const typename WellInterface<TypeTag>::FluidState& fs) const
    {
        if constexpr (Indices::oilEnabled) {
            return fs.pressure(FluidSystem::oilPhaseIdx);
        } else if constexpr (Indices::gasEnabled) {
            return fs.pressure(FluidSystem::gasPhaseIdx);
        } else {
            return fs.pressure(FluidSystem::waterPhaseIdx);
        }
    }

    template <typename TypeTag>
    template<class Value, class Callback>
    void
    WellInterface<TypeTag>::
    getMobility(const Simulator& simulator,
                const int perf,
                std::vector<Value>& mob,
                Callback& extendEval,
                [[maybe_unused]] DeferredLogger& deferred_logger) const
    {
        auto relpermArray = []()
                            {
                                if constexpr (std::is_same_v<Value, Scalar>) {
                                    return std::array<Scalar,3>{};
                                } else {
                                    return std::array<Eval,3>{};
                                }
                            };
        if (static_cast<std::size_t>(perf) >= this->well_cells_.size()) {
            OPM_THROW(std::invalid_argument,"The perforation index exceeds the size of the local containers - possibly getMobility was called with a global instead of a local perforation index!");
        }
        const int cell_idx = this->well_cells_[perf];
        assert (int(mob.size()) == this->num_components_);
        const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
        const auto& materialLawManager = simulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calculate its own
        // based on passing the saturation table index
        const int satid = this->saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if (satid == satid_elem) { // the same saturation number is used. i.e. just use the mobilty from the cell
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = extendEval(intQuants.mobility(phaseIdx));
            }
            if constexpr (has_solvent) {
                mob[Indices::contiSolventEqIdx] = extendEval(intQuants.solventMobility());
            }
        } else {
            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            auto relativePerms = relpermArray();
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = extendEval(relativePerms[phaseIdx] / intQuants.fluidState().viscosity(phaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if constexpr (has_solvent) {
                OPM_DEFLOG_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent", deferred_logger);
            }
        }

        if (this->isInjector() && !this->inj_fc_multiplier_.empty()) {
            const auto perf_ecl_index = this->perforationData()[perf].ecl_index;
            const auto& connections = this->well_ecl_.getConnections();
            const auto& connection = connections[perf_ecl_index];
            if (connection.filterCakeActive()) {
                for (auto& val : mob) {
                    val *= this->inj_fc_multiplier_[perf];
                }
            }
        }
    }


    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellStateWithTHPTargetProd(const Simulator& simulator,
                                     WellState<Scalar>& well_state,
                                     DeferredLogger& deferred_logger) const
    {
        OPM_TIMEFUNCTION();
        const auto& summary_state = simulator.vanguard().summaryState();

        auto bhp_at_thp_limit = computeBhpAtThpLimitProdWithAlq(
            simulator, summary_state, this->getALQ(well_state), deferred_logger, /*iterate_if_no_solution */ true);
        if (bhp_at_thp_limit) {
            std::vector<Scalar> rates(this->number_of_phases_, 0.0);
            if (thp_update_iterations) {
                computeWellRatesWithBhpIterations(simulator, *bhp_at_thp_limit,
                                                  rates, deferred_logger);
            } else {
                computeWellRatesWithBhp(simulator, *bhp_at_thp_limit,
                                        rates, deferred_logger);
            }
            auto& ws = well_state.well(this->name());
            ws.surface_rates = rates;
            ws.bhp = *bhp_at_thp_limit;
            ws.thp = this->getTHPConstraint(summary_state);
            return true;
        } else {
            return false;
        }
    }

    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    computeConnLevelProdInd(const FluidState& fs,
                            const std::function<Scalar(const Scalar)>& connPICalc,
                            const std::vector<Scalar>& mobility,
                            Scalar* connPI) const
    {
        const auto& pu = this->phaseUsage();
        const int   np = this->number_of_phases_;
        for (int p = 0; p < np; ++p) {
            // Note: E100's notion of PI value phase mobility includes
            // the reciprocal FVF.
            const auto connMob =
                mobility[this->flowPhaseToModelCompIdx(p)]
                    * fs.invB(this->flowPhaseToModelPhaseIdx(p)).value();

            connPI[p] = connPICalc(connMob);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
        {
            const auto io = pu.phase_pos[Oil];
            const auto ig = pu.phase_pos[Gas];

            const auto vapoil = connPI[ig] * fs.Rv().value();
            const auto disgas = connPI[io] * fs.Rs().value();

            connPI[io] += vapoil;
            connPI[ig] += disgas;
        }
    }


    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    computeConnLevelInjInd(const FluidState& fs,
                           const Phase preferred_phase,
                           const std::function<Scalar(const Scalar)>& connIICalc,
                           const std::vector<Scalar>& mobility,
                           Scalar* connII,
                           DeferredLogger& deferred_logger) const
    {
        // Assumes single phase injection
        const auto& pu = this->phaseUsage();

        auto phase_pos = 0;
        if (preferred_phase == Phase::GAS) {
            phase_pos = pu.phase_pos[Gas];
        }
        else if (preferred_phase == Phase::OIL) {
            phase_pos = pu.phase_pos[Oil];
        }
        else if (preferred_phase == Phase::WATER) {
            phase_pos = pu.phase_pos[Water];
        }
        else {
            OPM_DEFLOG_THROW(NotImplemented,
                             fmt::format("Unsupported Injector Type ({}) "
                                         "for well {} during connection I.I. calculation",
                                         static_cast<int>(preferred_phase), this->name()),
                             deferred_logger);
        }

        const auto mt     = std::accumulate(mobility.begin(), mobility.end(), 0.0);
        connII[phase_pos] = connIICalc(mt * fs.invB(this->flowPhaseToModelPhaseIdx(phase_pos)).value());
    }

    template<typename TypeTag>
    template<class GasLiftSingleWell>
    std::unique_ptr<GasLiftSingleWell> 
    WellInterface<TypeTag>::
    initializeGliftWellTest_(const Simulator& simulator,
                             WellState<Scalar>& well_state,
                             const GroupState<Scalar>& group_state,
                             const PhaseUsage& phase_usage,
                             GLiftEclWells& ecl_well_map,
                             DeferredLogger& deferred_logger)
    {
        // Instantiate group info object (without initialization) since it is needed in GasLiftSingleWell
        auto& comm = simulator.vanguard().grid().comm();
        ecl_well_map.try_emplace(this->name(),  &(this->wellEcl()), this->indexOfWell());
        GasLiftGroupInfo<Scalar> group_info {
                ecl_well_map,
                simulator.vanguard().schedule(),
                simulator.vanguard().summaryState(),
                simulator.episodeIndex(),
                simulator.model().newtonMethod().numIterations(),
                phase_usage,
                deferred_logger,
                well_state,
                group_state,
                comm,
                false
        };

        // Return GasLiftSingleWell object to use the wellTestALQ() function
        std::set<int> sync_groups;
        const auto& summary_state = simulator.vanguard().summaryState();
        return std::make_unique<GasLiftSingleWell>(*this, 
                                                    simulator, 
                                                    summary_state,
                                                    deferred_logger, 
                                                    well_state, 
                                                    group_state,
                                                    group_info, 
                                                    sync_groups, 
                                                    comm, 
                                                    false);
        
    }

} // namespace Opm

#endif
