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
            deferred_logger.info(ss.str());
            updateWellStateWithTarget(ebos_simulator, group_state, well_state, deferred_logger);
            updatePrimaryVariables(well_state, deferred_logger);
        }

        return changed;
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
                        const double simulation_time, WellState& well_state, const GroupState& group_state,
                        WellTestState& welltest_state, DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + this->name() + " is being tested for economic limits");

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
            const size_t original_number_closed_completions = welltest_state_temp.sizeCompletions();
            solveWellForTesting(simulator, well_state_copy, group_state, deferred_logger);
            std::vector<double> potentials;
            try {
                computeWellPotentials(simulator, well_state_copy, potentials, deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg = std::string("well ") + this->name() + std::string(": computeWellPotentials() failed during testing for re-opening: ") + e.what();
                deferred_logger.info(msg);
            }
            const int np = well_state_copy.numPhases();
            for (int p = 0; p < np; ++p) {
                ws.well_potentials[p] = std::abs(potentials[p]);
            }
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
        if (!converged) {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("Compute initial well solution for well " + this->name() + ". Failed to converge in "
                                  + std::to_string(max_iter) + " iterations");
            // the well operability system currently works only for producers in prediction mode
            if (this->shutUnsolvableWells())
                this->operability_status_.solvable = false;

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
        const bool old_well_operable = this->operability_status_.isOperable();
        checkWellOperability(ebosSimulator, well_state, deferred_logger);

        // only use inner well iterations for the first newton iterations.
        const int iteration_idx = ebosSimulator.model().newtonMethod().numIterations();
        bool converged = true;
        if (iteration_idx < param_.max_niter_inner_well_iter_)
            converged = this->iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);

        // unsolvable wells are treated as not operable and will not be solved for in this iteration.
        if (!converged) {
            if (this->shutUnsolvableWells())
                this->operability_status_.solvable = false;
        }
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
        auto& ws = well_state_copy.well(this->indexOfWell());

        // TODO: well state for this well is kind of all zero status
        // we should be able to provide a better initialization
        calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);

        updateWellOperability(ebos_simulator, well_state_copy, deferred_logger);

        if ( !this->isOperable() ) {
            const std::string msg = " well " + this->name() + " is not operable during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        updateWellStateWithTarget(ebos_simulator, group_state, well_state_copy, deferred_logger);

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
            // we need to populate the new well with potentials
            std::vector<double> potentials;
            try {
                computeWellPotentials(ebos_simulator, well_state_copy, potentials, deferred_logger);
            } catch (const std::exception& e) {
                const std::string msg2 = std::string("well ") + this->name() + std::string(": computeWellPotentials() failed during testing for re-opening: ") + e.what();
                deferred_logger.info(msg2);
            }
            const int np = well_state_copy.numPhases();
            for (int p = 0; p < np; ++p) {
                ws.well_potentials[p] = std::abs(potentials[p]);
            }
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

        updateWellOperability(ebos_simulator, well_state, deferred_logger);
    }


    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    shutUnsolvableWells() const
    {
        bool shut_unsolvable_wells = param_.shut_unsolvable_wells_;
        // the well operability system currently works only for producers in prediction mode
        return shut_unsolvable_wells && !this->isInjector() && this->underPredictionMode();
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& well_state,
                          DeferredLogger& deferred_logger)
    {
        this->operability_status_.reset();

        auto current_control = well_state.well(this->index_of_well_).production_cmode;
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
                well_state.wellRates(well_index)[p] = 0.0;
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
                double bhp = this->calculateBhpFromThp(well_state, rates, well, summaryState, this->getRefDensity(), deferred_logger);
                ws.bhp = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = ws.well_potentials[p];
                    }
                }
                break;
            }
            case Well::InjectorCMode::BHP:
            {
                ws.bhp = controls.bhp_limit;
                double total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate += well_state.wellRates(well_index)[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = ws.well_potentials[p];
                    }
                }
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
                    well_state.wellRates(well_index)[phasePos] = *target;
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
                double current_rate = -well_state.wellRates(well_index)[ pu.phase_pos[Oil] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] *= controls.oil_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                    const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
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
                        const std::vector<double> fractions = initialWellRateFractions(ebos_simulator, well_state);
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates(well_index)[p] = - fractions[p] * target / convert_coeff[p];
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
                    total_rate -= well_state.wellRates(well_index)[p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = -ws.well_potentials[p];
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
                double bhp = this->calculateBhpFromThp(well_state, rates, well, summaryState, this->getRefDensity(), deferred_logger);
                ws.bhp = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = -std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates(well_index)[p] = -ws.well_potentials[p];
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
                        well_state.wellRates(well_index)[p] *= scale;
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
        const int comp_idx_nz = this->flowPhaseToEbosCompIdx(nonzero_rate_index);
        for (int p = 0; p < this->number_of_phases_; ++p) {
            if (p != nonzero_rate_index) {
                const int comp_idx = this->flowPhaseToEbosCompIdx(p);
                double& rate = well_state.wellRates(this->index_of_well_)[p];
                rate = (initial_nonzero_rate/well_q_s[comp_idx_nz]) * (well_q_s[comp_idx]);
            }
        }
    }


} // namespace Opm
