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

#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <dune/common/version.hh>

namespace Opm
{


    template<typename TypeTag>
    WellInterface<TypeTag>::
    WellInterface(const Well& well,
                  const ParallelWellInfo& pw_info,
                  const int time_step,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components,
                  const int num_phases,
                  const int index_of_well,
                  const std::vector<PerforationData>& perf_data)
      : WellInterfaceIndices<FluidSystem,Indices,Scalar>(well,
                                                         pw_info,
                                                         time_step,
                                                         rate_converter,
                                                         pvtRegionIdx,
                                                         num_components,
                                                         num_phases,
                                                         index_of_well,
                                                         perf_data)
      , param_(param)
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
         const std::vector<double>& /* depth_arg */,
         const double gravity_arg,
         const int /* num_cells */,
         const std::vector< Scalar >& B_avg)
    {
        this->phase_usage_ = phase_usage_arg;
        this->gravity_ = gravity_arg;
        B_avg_ = B_avg;
    }




    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wpolymer() const
    {
        if constexpr (has_polymer) {
            auto injectorType = this->well_ecl_.injectorType();

            if (injectorType == InjectorType::WATER) {
                WellPolymerProperties polymer = this->well_ecl_.getPolymerProperties();
                const double polymer_injection_concentration = polymer.m_polymerConcentration;
                return polymer_injection_concentration;
            } else {
                // Not a water injection well => no polymer.
                return 0.0;
            }
        }

        return 0.0;
    }





    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wfoam() const
    {
        if constexpr (has_foam) {
            auto injectorType = this->well_ecl_.injectorType();

            if (injectorType == InjectorType::GAS) {
                WellFoamProperties fprop = this->well_ecl_.getFoamProperties();
                return fprop.m_foamConcentration;
            } else {
                // Not a gas injection well => no foam.
                return 0.0;
            }
        }

        return 0.0;
    }



    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wsalt() const
    {
        if constexpr (has_brine) {
            auto injectorType = this->well_ecl_.injectorType();

            if (injectorType == InjectorType::WATER) {
                WellBrineProperties fprop = this->well_ecl_.getBrineProperties();
                return fprop.m_saltConcentration;
            } else {
                // Not a water injection well => no salt (?).
                return 0.0;
            }
        }

        return 0.0;
    }

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wmicrobes() const
    {
      if constexpr (has_micp) {
          auto injectorType = this->well_ecl_.injectorType();

          if (injectorType == InjectorType::WATER) {
              WellMICPProperties microbes = this->well_ecl_.getMICPProperties();
              const double microbial_injection_concentration = microbes.m_microbialConcentration;
              return microbial_injection_concentration;
          } else {
              // Not a water injection well => no microbes.
              return 0.0;
          }
      }

      return 0.0;
    }

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    woxygen() const
    {
      if constexpr (has_micp) {
          auto injectorType = this->well_ecl_.injectorType();

          if (injectorType == InjectorType::WATER) {
              WellMICPProperties oxygen = this->well_ecl_.getMICPProperties();
              const double oxygen_injection_concentration = oxygen.m_oxygenConcentration;
              return oxygen_injection_concentration;
          } else {
              // Not a water injection well => no oxygen.
              return 0.0;
          }
      }

      return 0.0;
    }

    // The urea injection concentration is scaled down by a factor of 10, since its value
    // can be much bigger than 1 (not doing this slows the simulations). The
    // corresponding values are scaled accordingly in blackoilmicpmodules.hh when computing
    // the reactions and also when writing the output files (vtk and eclipse format, i.e.,
    // vtkblackoilmicpmodule.hh and ecloutputblackoilmodel.hh respectively).

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wurea() const
    {
      if constexpr (has_micp) {
          auto injectorType = this->well_ecl_.injectorType();

          if (injectorType == InjectorType::WATER) {
              WellMICPProperties urea = this->well_ecl_.getMICPProperties();
              const double urea_injection_concentration = urea.m_ureaConcentration / 10.; //Dividing by scaling factor 10
              return urea_injection_concentration;
          } else {
              // Not a water injection well => no urea.
              return 0.0;
          }
      }

      return 0.0;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellControl(const Simulator& ebos_simulator,
                      const IndividualOrGroup iog,
                      WellState& well_state,
                      const GroupState& group_state,
                      DeferredLogger& deferred_logger) /* const */
    {
        if (this->wellIsStopped()) {
            return false;
        }

        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto& schedule = ebos_simulator.vanguard().schedule();
        const auto& well = this->well_ecl_;
        auto& ws = well_state.well(this->index_of_well_);
        std::string from;
        if (well.isInjector()) {
            from = Well::InjectorCMode2String(ws.injection_cmode);
        } else {
            from = Well::ProducerCMode2String(ws.production_cmode);
        }
        bool oscillating = std::count(this->well_control_log_.begin(), this->well_control_log_.end(), from) >= param_.max_number_of_well_switches_;

        if (oscillating) {
            // only output frist time
            bool output = std::count(this->well_control_log_.begin(), this->well_control_log_.end(), from) == param_.max_number_of_well_switches_;
            if (output) {
                std::ostringstream ss;
                ss << "    The control model for well " << this->name()
                   << " is oscillating\n"
                   << "    We don't allow for more than "
                   << param_.max_number_of_well_switches_
                   << " switches. The control is kept at " << from;
                deferred_logger.info(ss.str());
                // add one more to avoid outputting the same info again
                this->well_control_log_.push_back(from);
            }
            return false;
        }
        bool changed = false;
        if (iog == IndividualOrGroup::Individual) {
            changed = this->checkIndividualConstraints(ws, summaryState);
        } else if (iog == IndividualOrGroup::Group) {
            changed = this->checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        } else {
            assert(iog == IndividualOrGroup::Both);
            changed = this->checkConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        }
        Parallel::Communication cc = ebos_simulator.vanguard().grid().comm();
        // checking whether control changed
        if (changed) {
            std::string to;
            if (well.isInjector()) {
                to = Well::InjectorCMode2String(ws.injection_cmode);
            } else {
                to = Well::ProducerCMode2String(ws.production_cmode);
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
            updateWellStateWithTarget(ebos_simulator, group_state, well_state, deferred_logger);
            updatePrimaryVariables(well_state, deferred_logger);
        }

        return changed;
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTesting(const Simulator& simulator,
                const double simulation_time,
                /* const */ WellState& well_state,
                const GroupState& group_state,
                WellTestState& well_test_state,
                DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + this->name() + " is being tested");

        WellState well_state_copy = well_state;
        auto& ws = well_state_copy.well(this->indexOfWell());

        updateWellStateWithTarget(simulator, group_state, well_state_copy, deferred_logger);
        calculateExplicitQuantities(simulator, well_state_copy, deferred_logger);
        updatePrimaryVariables(well_state_copy, deferred_logger);
        initPrimaryVariablesEvaluation();

        WellTestState welltest_state_temp;

        bool testWell = true;
        // if a well is closed because all completions are closed, we need to check each completion
        // individually. We first open all completions, then we close one by one by calling updateWellTestState
        // untill the number of closed completions do not increase anymore.
        while (testWell) {
            const size_t original_number_closed_completions = welltest_state_temp.num_closed_completions();
            bool converged = solveWellForTesting(simulator, well_state_copy, group_state, deferred_logger);
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

            std::vector<double> potentials;
            try {
                computeWellPotentials(simulator, well_state_copy, potentials, deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg = std::string("well ") + this->name() + std::string(": computeWellPotentials() failed during testing for re-opening: ") + e.what();
                deferred_logger.info(msg);
                return;
            }
            const int np = well_state_copy.numPhases();
            for (int p = 0; p < np; ++p) {
                ws.well_potentials[p] = std::abs(potentials[p]);
            }
            this->updateWellTestState(well_state_copy.well(this->indexOfWell()), simulation_time, /*writeMessageToOPMLog=*/ false, welltest_state_temp, deferred_logger);
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
            for (auto& completion : this->well_ecl_.getCompletions()) {
                if (!welltest_state_temp.completion_is_closed(this->name(), completion.first))
                    well_test_state.open_completion(this->name(), completion.first);
            }
            // set the status of the well_state to open
            ws.open();
            well_state = well_state_copy;
        }
    }




    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    iterateWellEquations(const Simulator& ebosSimulator,
                         const double dt,
                         WellState& well_state,
                         const GroupState& group_state,
                         DeferredLogger& deferred_logger)
    {
        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        const auto inj_controls = this->well_ecl_.isInjector() ? this->well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = this->well_ecl_.isProducer() ? this->well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        bool converged = false;
        try {
            converged = this->iterateWellEqWithControl(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
        } catch (NumericalIssue& e ) {
            const std::string msg = "Inner well iterations failed for well " + this->name() + " Treat the well as unconverged. ";
            deferred_logger.warning("INNER_ITERATION_FAILED", msg);
            converged = false;
        }
        return converged;
    }


    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    solveWellForTesting(const Simulator& ebosSimulator, WellState& well_state, const GroupState& group_state,
                        DeferredLogger& deferred_logger)
    {
        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const double dt = ebosSimulator.timeStepSize();
        const bool converged = iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("WellTest: Well equation for well " + this->name() +  " converged");
            return true;
        }
        const int max_iter = param_.max_welleq_iter_;
        deferred_logger.debug("WellTest: Well equation for well " + this->name() + " failed converging in "
                              + std::to_string(max_iter) + " iterations");
        well_state = well_state0;
        return false;
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellEquation(const Simulator& ebosSimulator,
                      WellState& well_state,
                      const GroupState& group_state,
                      DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable())
            return;

        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const double dt = ebosSimulator.timeStepSize();
        const bool converged = iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        if (!converged) {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("Compute initial well solution for well " + this->name() + ". Failed to converge in "
                                  + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    assembleWellEq(const Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   const GroupState& group_state,
                   DeferredLogger& deferred_logger)
    {
        const bool old_well_operable = this->operability_status_.isOperableAndSolvable();

        if (param_.check_well_operability_iter_)
            checkWellOperability(ebosSimulator, well_state, deferred_logger);

        // only use inner well iterations for the first newton iterations.
        const int iteration_idx = ebosSimulator.model().newtonMethod().numIterations();
        if (iteration_idx < param_.max_niter_inner_well_iter_) {
            this->operability_status_.solvable = true;
            bool converged = this->iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);

            // unsolvable wells are treated as not operable and will not be solved for in this iteration.
            if (!converged) {
                if (param_.shut_unsolvable_wells_)
                    this->operability_status_.solvable = false;
            }
        }

        const bool well_operable = this->operability_status_.isOperableAndSolvable();
        if (!well_operable && old_well_operable) {
            if (this->well_ecl_.getAutomaticShutIn()) {
                deferred_logger.info(" well " + this->name() + " gets SHUT during iteration ");
            } else {
                if (!this->wellIsStopped()) {
                    deferred_logger.info(" well " + this->name() + " gets STOPPED during iteration ");
                    this->stopWell();
                    changed_to_stopped_this_step_ = true;
                }
            }
        } else if (well_operable && !old_well_operable) {
            deferred_logger.info(" well " + this->name() + " gets REVIVED during iteration ");
            this->openWell();
            changed_to_stopped_this_step_ = false;
        }

        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        const auto inj_controls = this->well_ecl_.isInjector() ? this->well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = this->well_ecl_.isProducer() ? this->well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::addCellRates(RateVector& rates, int cellIdx) const
    {
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
    WellInterface<TypeTag>::volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const {
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
    checkWellOperability(const Simulator& ebos_simulator,
                         const WellState& well_state,
                         DeferredLogger& deferred_logger)
    {

        if (!param_.check_well_operability_) {
            return;
        }

        if (this->wellIsStopped() && !changed_to_stopped_this_step_) {
            return;
        }

        updateWellOperability(ebos_simulator, well_state, deferred_logger);
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& well_state,
                          DeferredLogger& deferred_logger)
    {
        this->operability_status_.resetOperability();

        bool thp_controled = this->isInjector() ? well_state.well(this->index_of_well_).injection_cmode == Well::InjectorCMode::THP:
                                              well_state.well(this->index_of_well_).production_cmode == Well::ProducerCMode::THP;
        bool bhp_controled = this->isInjector() ? well_state.well(this->index_of_well_).injection_cmode == Well::InjectorCMode::BHP:
                                              well_state.well(this->index_of_well_).production_cmode == Well::ProducerCMode::BHP;

        // Operability checking is not free
        // Only check wells under BHP and THP control
        if(bhp_controled || thp_controled) {
            updateIPR(ebos_simulator, deferred_logger);
            checkOperabilityUnderBHPLimit(well_state, ebos_simulator, deferred_logger);
        }
        // we do some extra checking for wells under THP control.
        if (thp_controled) {
            checkOperabilityUnderTHPLimit(ebos_simulator, well_state, deferred_logger);
        }
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              const GroupState& group_state,
                              WellState& well_state,
                              DeferredLogger& deferred_logger) const
    {

        // only bhp and wellRates are used to initilize the primaryvariables for standard wells
        const auto& well = this->well_ecl_;
        const int well_index = this->index_of_well_;
        auto& ws = well_state.well(well_index);
        const auto& pu = this->phaseUsage();
        const int np = well_state.numPhases();
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto& schedule = ebos_simulator.vanguard().schedule();

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

            switch(current) {
            case Well::InjectorCMode::RATE:
            {
                ws.surface_rates[phasePos] = controls.surface_rate;
                break;
            }

            case Well::InjectorCMode::RESV:
            {
                std::vector<double> convert_coeff(this->number_of_phases_, 1.0);
                this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);
                const double coeff = convert_coeff[phasePos];
                ws.surface_rates[phasePos] = controls.reservoir_rate/coeff;
                break;
            }

            case Well::InjectorCMode::THP:
            {
                auto rates = ws.surface_rates;
                double bhp = this->calculateBhpFromThp(well_state, rates, well, summaryState, this->getRefDensity(), deferred_logger);
                ws.bhp = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0)
                    ws.surface_rates = ws.well_potentials;

                break;
            }
            case Well::InjectorCMode::BHP:
            {
                ws.bhp = controls.bhp_limit;
                double total_rate = 0.0;
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
                const double efficiencyFactor = well.getEfficiencyFactor();
                std::optional<double> target =
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
        }
        //Producer
        else
        {
            const auto current = ws.production_cmode;
            const auto& controls = well.productionControls(summaryState);
            switch (current) {
            case Well::ProducerCMode::ORAT:
            {
                double current_rate = -ws.surface_rates[ pu.phase_pos[Oil] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.oil_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                double current_rate = -ws.surface_rates[ pu.phase_pos[Water] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.water_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
                    double control_fraction = fractions[pu.phase_pos[Water]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.water_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRAT:
            {
                double current_rate = -ws.surface_rates[pu.phase_pos[Gas] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.gas_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
                    double control_fraction = fractions[pu.phase_pos[Gas]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.gas_rate/control_fraction;
                        }
                    }
                }

                break;

            }
            case Well::ProducerCMode::LRAT:
            {
                double current_rate = -ws.surface_rates[ pu.phase_pos[Water] ]
                        - ws.surface_rates[ pu.phase_pos[Oil] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= controls.liquid_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
                    double control_fraction = fractions[pu.phase_pos[Water]] + fractions[pu.phase_pos[Oil]];
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
                OPM_DEFLOG_THROW(std::runtime_error, "CRAT control not supported " << this->name(), deferred_logger);
            }
            case Well::ProducerCMode::RESV:
            {
                std::vector<double> convert_coeff(this->number_of_phases_, 1.0);
                this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);
                double total_res_rate = 0.0;
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
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] = - fractions[p] * controls.resv_rate / convert_coeff[p];
                        }
                    }
                } else {
                    std::vector<double> hrates(this->number_of_phases_,0.);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        hrates[pu.phase_pos[Water]] = controls.water_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        hrates[pu.phase_pos[Oil]] = controls.oil_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        hrates[pu.phase_pos[Gas]] = controls.gas_rate;
                    }
                    std::vector<double> hrates_resv(this->number_of_phases_,0.);
                    this->rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, this->pvtRegionIdx_, hrates, hrates_resv);
                    double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            ws.surface_rates[p] *= target/total_res_rate;
                        }
                    } else {
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                double total_rate = 0.0;
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
                auto rates = ws.surface_rates;
                double bhp = this->calculateBhpFromThp(well_state, rates, well, summaryState, this->getRefDensity(), deferred_logger);
                ws.bhp = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = -std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] = -ws.well_potentials[p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRUP:
            {
                assert(well.isAvailableForGroupControl());
                const auto& group = schedule.getGroup(well.groupName(), this->currentStep());
                const double efficiencyFactor = well.getEfficiencyFactor();
                double scale = this->getGroupProductionTargetRate(group,
                                                          well_state,
                                                          group_state,
                                                          schedule,
                                                          summaryState,
                                                          efficiencyFactor);

                // we don't want to scale with zero and get zero rates.
                if (scale > 0) {
                    for (int p = 0; p<np; ++p) {
                        ws.surface_rates[p] *= scale;
                    }
                }
                break;
            }
            case Well::ProducerCMode::CMODE_UNDEFINED:
            case Well::ProducerCMode::NONE:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + this->name() , deferred_logger);
            }

                break;
            } // end of switch
        }
    }

    template<typename TypeTag>
    std::vector<double>
    WellInterface<TypeTag>::
    initialWellRateFractions(const Simulator& ebosSimulator, const WellState& well_state) const
    {
        const int np = this->number_of_phases_;
        std::vector<double> scaling_factor(np);
        const auto& ws = well_state.well(this->index_of_well_);

        double total_potentials = 0.0;
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
        double total_tw = 0;
        const int nperf = this->number_of_perforations_;
        for (int perf = 0; perf < nperf; ++perf) {
            total_tw += this->well_index_[perf];
        }
        for (int perf = 0; perf < nperf; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            const double well_tw_fraction = this->well_index_[perf] / total_tw;
            double total_mobility = 0.0;
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = this->flowPhaseToEbosPhaseIdx(p);
                total_mobility += fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value();
            }
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = this->flowPhaseToEbosPhaseIdx(p);
                scaling_factor[p] += well_tw_fraction * fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value() / total_mobility;
            }
        }
        return scaling_factor;
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateRates(const Simulator& ebosSimulator,
                         WellState& well_state,
                         DeferredLogger& deferred_logger) const
    {
        // Check if the rates of this well only are single-phase, do nothing
        // if more than one nonzero rate.
        auto& ws = well_state.well(this->index_of_well_);
        int nonzero_rate_index = -1;
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (ws.surface_rates[p] != 0.0) {
                if (nonzero_rate_index == -1) {
                    nonzero_rate_index = p;
                } else {
                    // More than one nonzero rate.
                    return;
                }
            }
        }
        if (nonzero_rate_index == -1) {
            // No nonzero rates.
            return;
        }

        // Calculate the rates that follow from the current primary variables.
        std::vector<double> well_q_s = computeCurrentWellRates(ebosSimulator, deferred_logger);

        // Set the currently-zero phase flows to be nonzero in proportion to well_q_s.
        const double initial_nonzero_rate = ws.surface_rates[nonzero_rate_index];
        const int comp_idx_nz = this->flowPhaseToEbosCompIdx(nonzero_rate_index);
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (p != nonzero_rate_index) {
                const int comp_idx = this->flowPhaseToEbosCompIdx(p);
                double& rate = ws.surface_rates[p];
                rate = (initial_nonzero_rate/well_q_s[comp_idx_nz]) * (well_q_s[comp_idx]);
            }
        }
    }


} // namespace Opm
