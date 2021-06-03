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
#include <opm/simulators/wells/TargetCalculator.hpp>

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
                  const int first_perf_index,
                  const std::vector<PerforationData>& perf_data)
      : WellInterfaceFluidSystem<FluidSystem>(well, pw_info, time_step, rate_converter,
                                              pvtRegionIdx, num_components, num_phases,
                                              index_of_well, first_perf_index, perf_data)
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
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const auto& pu = this->phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

        // for other phases return the index
        return phaseIdx;
    }

    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const auto& pu = this->phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return FluidSystem::waterPhaseIdx;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return FluidSystem::oilPhaseIdx;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return FluidSystem::gasPhaseIdx;

        // for other phases return the index
        return phaseIdx;
    }

    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const
    {
        const auto& pu = this->phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == compIdx)
            return pu.phase_pos[Water];
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == compIdx)
            return pu.phase_pos[Oil];
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == compIdx)
            return pu.phase_pos[Gas];

        // for other phases return the index
        return compIdx;
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
        std::string from;
        if (well.isInjector()) {
            from = Well::InjectorCMode2String(well_state.currentInjectionControl(this->index_of_well_));
        } else {
            from = Well::ProducerCMode2String(well_state.currentProductionControl(this->index_of_well_));
        }

        bool changed = false;
        if (iog == IndividualOrGroup::Individual) {
            changed = this->checkIndividualConstraints(well_state, summaryState);
        } else if (iog == IndividualOrGroup::Group) {
            changed = this->checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        } else {
            assert(iog == IndividualOrGroup::Both);
            changed = this->checkConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();

        // checking whether control changed
        if (changed) {
            std::string to;
            if (well.isInjector()) {
                to = Well::InjectorCMode2String(well_state.currentInjectionControl(this->index_of_well_));
            } else {
                to = Well::ProducerCMode2String(well_state.currentProductionControl(this->index_of_well_));
            }
            std::ostringstream ss;
            ss << "    Switching control mode for well " << this->name()
               << " from " << from
               << " to " <<  to;
            if (cc.size() > 1) {
               ss << " on rank " << cc.rank();
            }
            deferred_logger.info(ss.str());
            updateWellStateWithTarget(ebos_simulator, well_state, deferred_logger);
            updatePrimaryVariables(well_state, deferred_logger);
        }

        return changed;
    }



    template<typename TypeTag>
    template<class ValueType>
    ValueType
    WellInterface<TypeTag>::
    calculateBhpFromThp(const WellState &well_state,
                        const std::vector<ValueType>& rates,
                        const Well& well,
                        const SummaryState& summaryState,
                        DeferredLogger& deferred_logger) const
    {
        // TODO: when well is under THP control, the BHP is dependent on the rates,
        // the well rates is also dependent on the BHP, so it might need to do some iteration.
        // However, when group control is involved, change of the rates might impacts other wells
        // so iterations on a higher level will be required. Some investigation might be needed when
        // we face problems under THP control.

        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        const ValueType aqua = rates[Water];
        const ValueType liquid = rates[Oil];
        const ValueType vapour = rates[Gas];

        // pick the reference density
        // typically the reference in the top layer
        const double rho = getRefDensity();

        if (this->isInjector() )
        {
            const auto& controls = well.injectionControls(summaryState);
            const double vfp_ref_depth = this->vfp_properties_->getInj()->getTable(controls.vfp_table_number).getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(this->ref_depth_, vfp_ref_depth, rho, this->gravity_);
            return this->vfp_properties_->getInj()->bhp(controls.vfp_table_number, aqua, liquid, vapour, this->getTHPConstraint(summaryState)) - dp;
         }
         else if (this->isProducer()) {
             const auto& controls = well.productionControls(summaryState);
             const double vfp_ref_depth = this->vfp_properties_->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
             const double dp = wellhelpers::computeHydrostaticCorrection(this->ref_depth_, vfp_ref_depth, rho, this->gravity_);
             return this->vfp_properties_->getProd()->bhp(controls.vfp_table_number, aqua, liquid, vapour, this->getTHPConstraint(summaryState), this->getALQ(well_state)) - dp;
         }
         else {
             OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER for well " + this->name(), deferred_logger);
         }



    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTesting(const Simulator& simulator,
                const double simulation_time, const int report_step,
                const WellTestConfig::Reason testing_reason,
                /* const */ WellState& well_state,
                const GroupState& group_state,
                WellTestState& well_test_state,
                DeferredLogger& deferred_logger)
    {
        if (testing_reason == WellTestConfig::Reason::PHYSICAL) {
            wellTestingPhysical(simulator, simulation_time, report_step,
                                well_state, group_state, well_test_state, deferred_logger);
        }

        if (testing_reason == WellTestConfig::Reason::ECONOMIC) {
            wellTestingEconomic(simulator, simulation_time,
                                well_state, group_state, well_test_state, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTestingEconomic(const Simulator& simulator,
                        const double simulation_time, const WellState& well_state, const GroupState& group_state,
                        WellTestState& welltest_state, DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + this->name() + " is being tested for economic limits");

        WellState well_state_copy = well_state;

        updateWellStateWithTarget(simulator, well_state_copy, deferred_logger);
        calculateExplicitQuantities(simulator, well_state_copy, deferred_logger);
        updatePrimaryVariables(well_state_copy, deferred_logger);
        initPrimaryVariablesEvaluation();

        WellTestState welltest_state_temp;

        bool testWell = true;
        // if a well is closed because all completions are closed, we need to check each completion
        // individually. We first open all completions, then we close one by one by calling updateWellTestState
        // untill the number of closed completions do not increase anymore.
        while (testWell) {
            const size_t original_number_closed_completions = welltest_state_temp.sizeCompletions();
            solveWellForTesting(simulator, well_state_copy, group_state, deferred_logger);
            this->updateWellTestState(well_state_copy, simulation_time, /*writeMessageToOPMLog=*/ false, welltest_state_temp, deferred_logger);
            this->closeCompletions(welltest_state_temp);

            // Stop testing if the well is closed or shut due to all completions shut
            // Also check if number of completions has increased. If the number of closed completions do not increased
            // we stop the testing.
            // TODO: it can be tricky here, if the well is shut/closed due to other reasons
            if ( welltest_state_temp.sizeWells() > 0 ||
                (original_number_closed_completions == welltest_state_temp.sizeCompletions()) ) {
                    testWell = false; // this terminates the while loop
            }
        }

        // update wellTestState if the well test succeeds
        if (!welltest_state_temp.hasWellClosed(this->name(), WellTestConfig::Reason::ECONOMIC)) {
            welltest_state.openWell(this->name(), WellTestConfig::Reason::ECONOMIC);
            const std::string msg = std::string("well ") + this->name() + std::string(" is re-opened through ECONOMIC testing");
            deferred_logger.info(msg);

            // also reopen completions
            for (auto& completion : this->well_ecl_.getCompletions()) {
                if (!welltest_state_temp.hasCompletion(this->name(), completion.first)) {
                    welltest_state.dropCompletion(this->name(), completion.first);
                }
            }
        }
    }



    template<typename TypeTag>
    double
    WellInterface<TypeTag>::scalingFactor(const int phaseIdx) const
    {
        const auto& pu = this->phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return 0.01;
        if (has_solvent && phaseIdx == contiSolventEqIdx )
            return 0.01;

        // we should not come this far
        assert(false);
        return 1.0;
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

        return this->iterateWellEqWithControl(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
    }


    template<typename TypeTag>
    void
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
        } else {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("WellTest: Well equation for well " + this->name() + " failed converging in "
                          + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellEquation(const Simulator& ebosSimulator,
                      WellState& well_state,
                      const GroupState& group_state,
                      DeferredLogger& deferred_logger)
    {
        if (!this->isOperable())
            return;

        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const double dt = ebosSimulator.timeStepSize();
        const bool converged = iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("Compute initial well solution for well " + this->name() +  ". Converged");
        } else {
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

        checkWellOperability(ebosSimulator, well_state, deferred_logger);
        const int iterationIdx = ebosSimulator.model().newtonMethod().numIterations();

        if (this->useInnerIterations() && iterationIdx < 3) {
            this->iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
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
    wellTestingPhysical(const Simulator& ebos_simulator,
                        const double /* simulation_time */, const int /* report_step */,
                        WellState& well_state,
                        const GroupState& group_state,
                        WellTestState& welltest_state,
                        DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + this->name() + " is being tested for physical limits");

        // some most difficult things are the explicit quantities, since there is no information
        // in the WellState to do a decent initialization

        // TODO: Let us assume that the simulator is updated

        // Let us try to do a normal simualtion running, to keep checking the operability status
        // If the well is not operable during any of the time. It means it does not pass the physical
        // limit test.

        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState well_state_copy = well_state;

        // TODO: well state for this well is kind of all zero status
        // we should be able to provide a better initialization
        calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);

        updateWellOperability(ebos_simulator, well_state_copy, deferred_logger);

        if ( !this->isOperable() ) {
            const std::string msg = " well " + this->name() + " is not operable during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        updateWellStateWithTarget(ebos_simulator, well_state_copy, deferred_logger);

        calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);

        const double dt = ebos_simulator.timeStepSize();
        const bool converged = this->iterateWellEquations(ebos_simulator, dt, well_state_copy, group_state, deferred_logger);

        if (!converged) {
            const std::string msg = " well " + this->name() + " did not get converged during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        if (this->isOperable() ) {
            welltest_state.openWell(this->name(), WellTestConfig::PHYSICAL );
            const std::string msg = " well " + this->name() + " is re-opened through well testing for physical reason";
            deferred_logger.info(msg);
            well_state = well_state_copy;
        } else {
            const std::string msg = " well " + this->name() + " is not operable during well testing for physical reason";
            deferred_logger.debug(msg);
        }
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkWellOperability(const Simulator& ebos_simulator,
                         const WellState& well_state,
                         DeferredLogger& deferred_logger)
    {

        const bool checkOperability = EWOMS_GET_PARAM(TypeTag, bool, EnableWellOperabilityCheck);
        if (!checkOperability) {
            return;
        }

        // focusing on PRODUCER for now
        if (this->isInjector()) {
            return;
        }

        if (!this->underPredictionMode() ) {
            return;
        }

        if (this->wellIsStopped() && !changed_to_stopped_this_step_) {
            return;
        }

        const bool old_well_operable = this->operability_status_.isOperable();

        updateWellOperability(ebos_simulator, well_state, deferred_logger);

        const bool well_operable = this->operability_status_.isOperable();

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
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& well_state,
                          DeferredLogger& deferred_logger)
    {
        this->operability_status_.reset();

        auto current_control = well_state.currentProductionControl(this->index_of_well_);
        // Operability checking is not free
        // Only check wells under BHP and THP control
        if(current_control == Well::ProducerCMode::BHP || current_control == Well::ProducerCMode::THP) {
            updateIPR(ebos_simulator, deferred_logger);
            checkOperabilityUnderBHPLimitProducer(well_state, ebos_simulator, deferred_logger);
        }
        // we do some extra checking for wells under THP control.
        if (current_control == Well::ProducerCMode::THP) {
            checkOperabilityUnderTHPLimitProducer(ebos_simulator, well_state, deferred_logger);
        }
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              WellState& well_state,
                              DeferredLogger& deferred_logger) const
    {

        // only bhp and wellRates are used to initilize the primaryvariables for standard wells
        const auto& well = this->well_ecl_;
        const int well_index = this->index_of_well_;
        const auto& pu = this->phaseUsage();
        const int np = well_state.numPhases();
        const auto& summaryState = ebos_simulator.vanguard().summaryState();

        if (this->wellIsStopped()) {
            for (int p = 0; p<np; ++p) {
                well_state.wellRates(well_index)[p] = 0.0;
            }
            well_state.update_thp(well_index, 0.0);
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

            auto current = well_state.currentInjectionControl(well_index);

            switch(current) {
            case Well::InjectorCMode::RATE:
            {
                well_state.wellRates(well_index)[phasePos] = controls.surface_rate;
                break;
            }

            case Well::InjectorCMode::RESV:
            {
                std::vector<double> convert_coeff(this->number_of_phases_, 1.0);
                this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);
                const double coeff = convert_coeff[phasePos];
                well_state.wellRates(well_index)[phasePos] = controls.reservoir_rate/coeff;
                break;
            }

            case Well::InjectorCMode::THP:
            {
                std::vector<double> rates(3, 0.0);
                for (int p = 0; p<np; ++p) {
                    rates[p] = well_state.wellRates(well_index)[p];
                }
                double bhp = calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
                well_state.update_bhp(well_index, bhp);

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::InjectorCMode::BHP:
            {
                well_state.update_bhp(well_index, controls.bhp_limit);
                double total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate += well_state.wellRates(well_index)[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::InjectorCMode::GRUP:
            {
                //do nothing at the moment
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
            auto current = well_state.currentProductionControl(well_index);
            const auto& controls = well.productionControls(summaryState);
            switch (current) {
            case Well::ProducerCMode::ORAT:
            {
                double current_rate = -well_state.wellRates(well_index)[ pu.phase_pos[Oil] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] *= controls.oil_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * controls.oil_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::WRAT:
            {
                double current_rate = -well_state.wellRates(well_index)[ pu.phase_pos[Water] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] *= controls.water_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Water]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * controls.water_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRAT:
            {
                double current_rate = -well_state.wellRates(well_index)[pu.phase_pos[Gas] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] *= controls.gas_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Gas]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * controls.gas_rate/control_fraction;
                        }
                    }
                }

                break;

            }
            case Well::ProducerCMode::LRAT:
            {
                double current_rate = -well_state.wellRates(well_index)[ pu.phase_pos[Water] ]
                        - well_state.wellRates(well_index)[ pu.phase_pos[Oil] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] *= controls.liquid_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Water]] + fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * controls.liquid_rate / control_fraction;
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
                    total_res_rate -= well_state.wellRates(well_index)[p] * convert_coeff[p];
                }
                if (controls.prediction_mode) {
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] *= controls.resv_rate/total_res_rate;
                        }
                    } else {
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * controls.resv_rate / convert_coeff[p];
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
                            well_state.wellRates(well_index)[p] *= target/total_res_rate;
                        }
                    } else {
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state.wellPotentials());
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * target / convert_coeff[p];
                        }
                    }

                }
                break;
            }
            case Well::ProducerCMode::BHP:
            {
                well_state.update_bhp(well_index, controls.bhp_limit);
                double total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate -= well_state.wellRates(well_index)[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = -well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::THP:
            {
                std::vector<double> rates(3, 0.0);
                for (int p = 0; p<np; ++p) {
                    rates[p] = well_state.wellRates(well_index)[p];
                }
                double bhp = calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
                well_state.update_bhp(well_index, bhp);

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = -std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = -well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRUP:
            {
                //do nothing at the moment
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
    initialWellRateFractions(const Simulator& ebosSimulator, const std::vector<double>& potentials) const
    {
        const int np = this->number_of_phases_;
        std::vector<double> scaling_factor(np);

        double total_potentials = 0.0;
        for (int p = 0; p<np; ++p) {
            total_potentials += potentials[this->index_of_well_*np + p];
        }
        if (total_potentials > 0) {
            for (int p = 0; p<np; ++p) {
                scaling_factor[p] = potentials[this->index_of_well_*np + p] / total_potentials;
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
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(p);
                total_mobility += fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value();
            }
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(p);
                scaling_factor[p] += well_tw_fraction * fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value() / total_mobility;
            }
        }
        return scaling_factor;
    }



    template <typename TypeTag>
    template <class EvalWell, class BhpFromThpFunc>
    void
    WellInterface<TypeTag>::assembleControlEqInj(const WellState& well_state,
                                                 const GroupState& group_state,
                                                 const Schedule& schedule,
                                                 const SummaryState& summaryState,
                                                 const Well::InjectionControls& controls,
                                                 const EvalWell& bhp,
                                                 const EvalWell& injection_rate,
                                                 BhpFromThpFunc bhp_from_thp,
                                                 EvalWell& control_eq,
                                                 DeferredLogger& deferred_logger)
    {
        auto current = well_state.currentInjectionControl(this->index_of_well_);
        const InjectorType injectorType = controls.injector_type;
        const auto& pu = this->phaseUsage();
        const double efficiencyFactor = this->well_ecl_.getEfficiencyFactor();

        switch (current) {
        case Well::InjectorCMode::RATE: {
            control_eq = injection_rate - controls.surface_rate;
            break;
        }
        case Well::InjectorCMode::RESV: {
            std::vector<double> convert_coeff(this->number_of_phases_, 1.0);
            this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);

            double coeff = 1.0;

            switch (injectorType) {
            case InjectorType::WATER: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Aqua]];
                break;
            }
            case InjectorType::OIL: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Liquid]];
                break;
            }
            case InjectorType::GAS: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Vapour]];
                break;
            }
            default:
                throw("Expected WATER, OIL or GAS as type for injectors " + this->well_ecl_.name());
            }

            control_eq = coeff * injection_rate - controls.reservoir_rate;
            break;
        }
        case Well::InjectorCMode::THP: {
            control_eq = bhp - bhp_from_thp();
            break;
        }
        case Well::InjectorCMode::BHP: {
            control_eq = bhp - controls.bhp_limit;
            break;
        }
        case Well::InjectorCMode::GRUP: {
            assert(this->well_ecl_.isAvailableForGroupControl());
            const auto& group = schedule.getGroup(this->well_ecl_.groupName(), this->current_step_);
            getGroupInjectionControl(group,
                                     well_state,
                                     group_state,
                                     schedule,
                                     summaryState,
                                     injectorType,
                                     bhp,
                                     injection_rate,
                                     control_eq,
                                     efficiencyFactor,
                                     deferred_logger);
            break;
        }
        case Well::InjectorCMode::CMODE_UNDEFINED: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + this->name(), deferred_logger);
        }
        }
    }




    template <typename TypeTag>
    template <class EvalWell, class BhpFromThpFunc>
    void
    WellInterface<TypeTag>::assembleControlEqProd(const WellState& well_state,
                                                  const GroupState& group_state,
                                                  const Schedule& schedule,
                                                  const SummaryState& summaryState,
                                                  const Well::ProductionControls& controls,
                                                  const EvalWell& bhp,
                                                  const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                                  BhpFromThpFunc bhp_from_thp,
                                                  EvalWell& control_eq,
                                                  DeferredLogger& deferred_logger)
    {
        auto current = well_state.currentProductionControl(this->index_of_well_);
        const auto& pu = this->phaseUsage();
        const double efficiencyFactor = this->well_ecl_.getEfficiencyFactor();

        switch (current) {
        case Well::ProducerCMode::ORAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Liquid];
            control_eq = rate - controls.oil_rate;
            break;
        }
        case Well::ProducerCMode::WRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Aqua];
            control_eq = rate - controls.water_rate;
            break;
        }
        case Well::ProducerCMode::GRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Vapour];
            control_eq = rate - controls.gas_rate;
            break;
        }
        case Well::ProducerCMode::LRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            EvalWell rate = -rates[BlackoilPhases::Aqua] - rates[BlackoilPhases::Liquid];
            control_eq = rate - controls.liquid_rate;
            break;
        }
        case Well::ProducerCMode::CRAT: {
            OPM_DEFLOG_THROW(std::runtime_error, "CRAT control not supported " << this->name(), deferred_logger);
        }
        case Well::ProducerCMode::RESV: {
            auto total_rate = rates[0]; // To get the correct type only.
            total_rate = 0.0;
            std::vector<double> convert_coeff(this->number_of_phases_, 1.0);
            this->rateConverter_.calcCoeff(/*fipreg*/ 0, this->pvtRegionIdx_, convert_coeff);
            for (int phase = 0; phase < 3; ++phase) {
                if (pu.phase_used[phase]) {
                    const int pos = pu.phase_pos[phase];
                    total_rate -= rates[phase] * convert_coeff[pos]; // Note different indices.
                }
            }
            if (controls.prediction_mode) {
                control_eq = total_rate - controls.resv_rate;
            } else {
                std::vector<double> hrates(this->number_of_phases_, 0.);
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    hrates[pu.phase_pos[Water]] = controls.water_rate;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    hrates[pu.phase_pos[Oil]] = controls.oil_rate;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    hrates[pu.phase_pos[Gas]] = controls.gas_rate;
                }
                std::vector<double> hrates_resv(this->number_of_phases_, 0.);
                this->rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, this->pvtRegionIdx_, hrates, hrates_resv);
                double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                control_eq = total_rate - target;
            }
            break;
        }
        case Well::ProducerCMode::BHP: {
            control_eq = bhp - controls.bhp_limit;
            break;
        }
        case Well::ProducerCMode::THP: {
            control_eq = bhp - bhp_from_thp();
            break;
        }
        case Well::ProducerCMode::GRUP: {
            assert(this->well_ecl_.isAvailableForGroupControl());
            const auto& group = schedule.getGroup(this->well_ecl_.groupName(), this->current_step_);
            // Annoying thing: the rates passed to this function are
            // always of size 3 and in canonical (for PhaseUsage)
            // order. This is what is needed for VFP calculations if
            // they are required (THP controlled well). But for the
            // group production control things we must pass only the
            // active phases' rates.
            std::vector<EvalWell> active_rates(pu.num_phases);
            for (int canonical_phase = 0; canonical_phase < 3; ++canonical_phase) {
                if (pu.phase_used[canonical_phase]) {
                    active_rates[pu.phase_pos[canonical_phase]] = rates[canonical_phase];
                }
            }
            getGroupProductionControl(group, well_state, group_state, schedule, summaryState, bhp, active_rates, control_eq, efficiencyFactor);
            break;
        }
        case Well::ProducerCMode::CMODE_UNDEFINED: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + this->name(), deferred_logger);
        }
        case Well::ProducerCMode::NONE: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + this->name(), deferred_logger);
        }
        }
    }



    template <typename TypeTag>
    template <class EvalWell>
    void
    WellInterface<TypeTag>::getGroupInjectionControl(const Group& group,
                                                      const WellState& well_state,
                                                     const GroupState& group_state,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const InjectorType& injectorType,
                                                      const EvalWell& bhp,
                                                      const EvalWell& injection_rate,
                                                      EvalWell& control_eq,
                                                      double efficiencyFactor,
                                                      DeferredLogger& deferred_logger)
    {
        // Setting some defaults to silence warnings below.
        // Will be overwritten in the switch statement.
        Phase injectionPhase = Phase::WATER;
        switch (injectorType) {
        case InjectorType::WATER:
        {
            injectionPhase = Phase::WATER;
            break;
        }
        case InjectorType::OIL:
        {
            injectionPhase = Phase::OIL;
            break;
        }
        case InjectorType::GAS:
        {
            injectionPhase = Phase::GAS;
            break;
        }
        default:
            // Should not be here.
            assert(false);
        }

        auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
        if (currentGroupControl == Group::InjectionCMode::FLD ||
            currentGroupControl == Group::InjectionCMode::NONE) {
            if (!group.injectionGroupControlAvailable(injectionPhase)) {
                // We cannot go any further up the hierarchy. This could
                // be the FIELD group, or any group for which this has
                // been set in GCONINJE or GCONPROD. If we are here
                // anyway, it is likely that the deck set inconsistent
                // requirements, such as GRUP control mode on a well with
                // no appropriate controls defined on any of its
                // containing groups. We will therefore use the wells' bhp
                // limit equation as a fallback.
                const auto& controls = this->well_ecl_.injectionControls(summaryState);
                control_eq = bhp - controls.bhp_limit;
                return;
            } else {
                // Inject share of parents control
                const auto& parent = schedule.getGroup( group.parent(), this->current_step_ );
                efficiencyFactor *= group.getGroupEfficiencyFactor();
                getGroupInjectionControl(parent, well_state, group_state, schedule, summaryState, injectorType, bhp, injection_rate, control_eq, efficiencyFactor, deferred_logger);
                return;
            }
        }

        efficiencyFactor *= group.getGroupEfficiencyFactor();
        const auto& well = this->well_ecl_;
        const auto pu = this->phaseUsage();

        if (!group.isInjectionGroup()) {
            // use bhp as control eq and let the updateControl code find a valid control
            const auto& controls = well.injectionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        }

        // If we are here, we are at the topmost group to be visited in the recursion.
        // This is the group containing the control we will check against.

        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(this->phaseUsage().num_phases, 1.0);
        this->rateConverter_.calcCoeff(0, this->pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        double sales_target = 0;
        if (schedule[this->current_step_].gconsale().has(group.name())) {
            const auto& gconsale = schedule[this->current_step_].gconsale().get(group.name(), summaryState);
            sales_target = gconsale.sales_target;
        }
        WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu, resv_coeff, group.name(), sales_target, group_state, injectionPhase, deferred_logger);
        WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, this->current_step_, this->guide_rate_, tcalc.guideTargetMode(), pu, false, injectionPhase);

        auto localFraction = [&](const std::string& child) {
            return fcalc.localFraction(child, "");
        };

        auto localReduction = [&](const std::string& group_name) {
            const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
            return tcalc.calcModeRateFromRates(groupTargetReductions);
        };

        const double orig_target = tcalc.groupTarget(group.injectionControls(injectionPhase, summaryState), deferred_logger);
        const auto chain = WellGroupHelpers::groupChainTopBot(this->name(), group.name(), schedule, this->current_step_);
        // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
        const size_t num_ancestors = chain.size() - 1;
        double target = orig_target;
        for (size_t ii = 0; ii < num_ancestors; ++ii) {
            if ((ii == 0) || this->guide_rate_->has(chain[ii], injectionPhase)) {
                // Apply local reductions only at the control level
                // (top) and for levels where we have a specified
                // group guide rate.
                target -= localReduction(chain[ii]);
            }
            target *= localFraction(chain[ii+1]);
        }
        // Avoid negative target rates coming from too large local reductions.
        const double target_rate = std::max(0.0, target / efficiencyFactor);
        const auto current_rate = injection_rate; // Switch sign since 'rates' are negative for producers.
        control_eq = current_rate - target_rate;
    }



    template <typename TypeTag>
    template <class EvalWell>
    void
    WellInterface<TypeTag>::getGroupProductionControl(const Group& group,
                                                      const WellState& well_state,
                                                      const GroupState& group_state,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const EvalWell& bhp,
                                                      const std::vector<EvalWell>& rates,
                                                      EvalWell& control_eq,
                                                      double efficiencyFactor)
    {
        const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
        if (currentGroupControl == Group::ProductionCMode::FLD ||
            currentGroupControl == Group::ProductionCMode::NONE) {
            if (!group.productionGroupControlAvailable()) {
                // We cannot go any further up the hierarchy. This could
                // be the FIELD group, or any group for which this has
                // been set in GCONINJE or GCONPROD. If we are here
                // anyway, it is likely that the deck set inconsistent
                // requirements, such as GRUP control mode on a well with
                // no appropriate controls defined on any of its
                // containing groups. We will therefore use the wells' bhp
                // limit equation as a fallback.
                const auto& controls = this->well_ecl_.productionControls(summaryState);
                control_eq = bhp - controls.bhp_limit;
                return;
            } else {
                // Produce share of parents control
                const auto& parent = schedule.getGroup( group.parent(), this->current_step_ );
                efficiencyFactor *= group.getGroupEfficiencyFactor();
                getGroupProductionControl(parent, well_state, group_state, schedule, summaryState, bhp, rates, control_eq, efficiencyFactor);
                return;
            }
        }

        efficiencyFactor *= group.getGroupEfficiencyFactor();
        const auto& well = this->well_ecl_;
        const auto pu = this->phaseUsage();

        if (!group.isProductionGroup()) {
            // use bhp as control eq and let the updateControl code find a valid control
            const auto& controls = well.productionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        }

        // If we are here, we are at the topmost group to be visited in the recursion.
        // This is the group containing the control we will check against.

        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(this->phaseUsage().num_phases, 1.0);
        this->rateConverter_.calcCoeff(0, this->pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        // gconsale may adjust the grat target.
        // the adjusted rates is send to the targetCalculator
        double gratTargetFromSales = 0.0;
        if (group_state.has_grat_sales_target(group.name()))
            gratTargetFromSales = group_state.grat_sales_target(group.name());

        WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff, gratTargetFromSales);
        WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, this->current_step_, this->guide_rate_, tcalc.guideTargetMode(), pu, true, Phase::OIL);

        auto localFraction = [&](const std::string& child) {
            return fcalc.localFraction(child, "");
        };

        auto localReduction = [&](const std::string& group_name) {
            const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
            return tcalc.calcModeRateFromRates(groupTargetReductions);
        };

        const double orig_target = tcalc.groupTarget(group.productionControls(summaryState));
        const auto chain = WellGroupHelpers::groupChainTopBot(this->name(), group.name(), schedule, this->current_step_);
        // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
        const size_t num_ancestors = chain.size() - 1;
        double target = orig_target;
        for (size_t ii = 0; ii < num_ancestors; ++ii) {
            if ((ii == 0) || this->guide_rate_->has(chain[ii])) {
                // Apply local reductions only at the control level
                // (top) and for levels where we have a specified
                // group guide rate.
                target -= localReduction(chain[ii]);
            }
            target *= localFraction(chain[ii+1]);
        }
        // Avoid negative target rates coming from too large local reductions.
        const double target_rate = std::max(0.0, target / efficiencyFactor);
        const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
        control_eq = current_rate - target_rate;
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
        int nonzero_rate_index = -1;
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (well_state.wellRates(this->index_of_well_)[p] != 0.0) {
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
        const double initial_nonzero_rate = well_state.wellRates(this->index_of_well_)[nonzero_rate_index];
        const int comp_idx_nz = flowPhaseToEbosCompIdx(nonzero_rate_index);
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (p != nonzero_rate_index) {
                const int comp_idx = flowPhaseToEbosCompIdx(p);
                double& rate = well_state.wellRates(this->index_of_well_)[p];
                rate = (initial_nonzero_rate/well_q_s[comp_idx_nz]) * (well_q_s[comp_idx]);
            }
        }
    }


} // namespace Opm
