/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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

#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/MSW/Segment.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/material/densead/EvaluationFormat.hpp>

#include <opm/simulators/wells/MultisegmentWellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <algorithm>
#include <cstddef>
#include <string>

#if HAVE_CUDA || HAVE_OPENCL
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#endif

namespace Opm
{


    template <typename TypeTag>
    MultisegmentWell<TypeTag>::
    MultisegmentWell(const Well& well,
                     const ParallelWellInfo& pw_info,
                     const int time_step,
                     const ModelParameters& param,
                     const RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components,
                     const int num_phases,
                     const int index_of_well,
                     const std::vector<PerforationData>& perf_data)
    : Base(well, pw_info, time_step, param, rate_converter, pvtRegionIdx, num_components, num_phases, index_of_well, perf_data)
    , MSWEval(static_cast<WellInterfaceIndices<FluidSystem,Indices,Scalar>&>(*this))
    , regularize_(false)
    , segment_fluid_initial_(this->numberOfSegments(), std::vector<double>(this->num_components_, 0.0))
    {
        // not handling solvent or polymer for now with multisegment well
        if constexpr (has_solvent) {
            OPM_THROW(std::runtime_error, "solvent is not supported by multisegment well yet");
        }

        if constexpr (has_polymer) {
            OPM_THROW(std::runtime_error, "polymer is not supported by multisegment well yet");
        }

        if constexpr (Base::has_energy) {
            OPM_THROW(std::runtime_error, "energy is not supported by multisegment well yet");
        }

        if constexpr (Base::has_foam) {
            OPM_THROW(std::runtime_error, "foam is not supported by multisegment well yet");
        }

        if constexpr (Base::has_brine) {
            OPM_THROW(std::runtime_error, "brine is not supported by multisegment well yet");
        }

        if constexpr (Base::has_watVapor) {
            OPM_THROW(std::runtime_error, "water evaporation is not supported by multisegment well yet");
        }

        if(this->rsRvInj() > 0) {
            OPM_THROW(std::runtime_error,
                      "dissolved gas/ vapporized oil in injected oil/gas not supported by multisegment well yet."
                      " \n See  (WCONINJE item 10 / WCONHIST item 8)");
        }
        if constexpr (!Indices::oilEnabled && Indices::numPhases > 1) {
            OPM_THROW(std::runtime_error, "water + gas case not supported by multisegment well yet");
        }

        this->thp_update_iterations = true;
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells,
         const std::vector< Scalar >& B_avg,
         const bool changed_to_open_this_step)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, num_cells, B_avg, changed_to_open_this_step);

        // TODO: for StandardWell, we need to update the perf depth here using depth_arg.
        // for MultisegmentWell, it is much more complicated.
        // It can be specified directly, it can be calculated from the segment depth,
        // it can also use the cell center, which is the same for StandardWell.
        // For the last case, should we update the depth with the depth_arg? For the
        // future, it can be a source of wrong result with Multisegment well.
        // An indicator from the opm-parser should indicate what kind of depth we should use here.

        // \Note: we do not update the depth here. And it looks like for now, we only have the option to use
        // specified perforation depth
        this->initMatrixAndVectors(num_cells);

        // calculate the depth difference between the perforations and the perforated grid block
        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            this->cell_perforation_depth_diffs_[perf] = depth_arg[cell_idx] - this->perf_depth_[perf];
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initPrimaryVariablesEvaluation()
    {
        this->primary_variables_.init();
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updatePrimaryVariables(const SummaryState& summary_state,
                           const WellState& well_state,
                           DeferredLogger& /* deferred_logger */)
    {
        const bool stop_or_zero_rate_target = this->stopppedOrZeroRateTarget(summary_state, well_state);
        this->primary_variables_.update(well_state, stop_or_zero_rate_target);
    }






    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              const GroupState& group_state,
                              WellState& well_state,
                              DeferredLogger&  deferred_logger) const
    {
        Base::updateWellStateWithTarget(ebos_simulator, group_state, well_state, deferred_logger);
        // scale segment rates based on the wellRates
        // and segment pressure based on bhp
        this->scaleSegmentRatesWithWellRates(this->segments_.inlets(),
                                             this->segments_.perforations(),
                                             well_state);
        this->scaleSegmentPressuresWithBhp(well_state);
    }





    template <typename TypeTag>
    ConvergenceReport
    MultisegmentWell<TypeTag>::
    getWellConvergence(const SummaryState& /* summary_state */,
                       const WellState& well_state,
                       const std::vector<double>& B_avg,
                       DeferredLogger& deferred_logger,
                       const bool relax_tolerance) const
    {
        return this->MSWEval::getWellConvergence(well_state,
                                                 B_avg,
                                                 deferred_logger,
                                                 this->param_.max_residual_allowed_,
                                                 this->param_.tolerance_wells_,
                                                 this->param_.relaxed_tolerance_flow_well_,
                                                 this->param_.tolerance_pressure_ms_wells_,
                                                 this->param_.relaxed_tolerance_pressure_ms_well_,
                                                 relax_tolerance, 
                                                 this->wellIsStopped());

    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) {
            return;
        }

        if (this->param_.matrix_add_well_contributions_) {
            // Contributions are already in the matrix itself
            return;
        }

        this->linSys_.apply(x, Ax);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(BVector& r) const
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) {
            return;
        }

        this->linSys_.apply(r);
    }



    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const SummaryState& summary_state,
                                          const BVector& x,
                                          WellState& well_state,
                                          DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) {
            return;
        }

        BVectorWell xw(1);
        this->linSys_.recoverSolutionWell(x, xw);
        updateWellState(summary_state, xw, well_state, deferred_logger);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials,
                          DeferredLogger& deferred_logger)
    {
        const auto [compute_potential, bhp_controlled_well] =
            this->WellInterfaceGeneric::computeWellPotentials(well_potentials, well_state);

        if (!compute_potential) {
            return;
        }

        debug_cost_counter_ = 0;
        bool converged_implicit = false;
        if (this->param_.local_well_solver_control_switching_) {
            converged_implicit = computeWellPotentialsImplicit(ebosSimulator, well_potentials, deferred_logger);
            if (!converged_implicit) {
                deferred_logger.debug("Implicit potential calculations failed for well "
                                       + this->name() + ",  reverting to original aproach.");
            }
        }
        if (!converged_implicit) {
            // does the well have a THP related constraint?
            const auto& summaryState = ebosSimulator.vanguard().summaryState();
            if (!Base::wellHasTHPConstraints(summaryState) || bhp_controlled_well) {
                computeWellRatesAtBhpLimit(ebosSimulator, well_potentials, deferred_logger);
            } else {
                well_potentials = computeWellPotentialWithTHP(
                    well_state, ebosSimulator, deferred_logger);
            }
        }
        deferred_logger.debug("Cost in iterations of finding well potential for well "
                              + this->name() + ": " + std::to_string(debug_cost_counter_));

        this->checkNegativeWellPotentials(well_potentials,
                                          this->param_.check_well_operability_,
                                          deferred_logger);
    }




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesAtBhpLimit(const Simulator& ebosSimulator,
                               std::vector<double>& well_flux,
                               DeferredLogger& deferred_logger) const
    {
        if (this->well_ecl_.isInjector()) {
            const auto controls = this->well_ecl_.injectionControls(ebosSimulator.vanguard().summaryState());
            computeWellRatesWithBhpIterations(ebosSimulator, controls.bhp_limit, well_flux, deferred_logger);
        } else {
            const auto controls = this->well_ecl_.productionControls(ebosSimulator.vanguard().summaryState());
            computeWellRatesWithBhpIterations(ebosSimulator, controls.bhp_limit, well_flux, deferred_logger);
        }
    }

    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const double& bhp,
                            std::vector<double>& well_flux,
                            DeferredLogger& deferred_logger) const
    {

        const int np = this->number_of_phases_;

        well_flux.resize(np, 0.0);
        const bool allow_cf = this->getAllowCrossFlow();
        const int nseg = this->numberOfSegments();
        const WellState& well_state = ebosSimulator.problem().wellModel().wellState();
        const auto& ws = well_state.well(this->indexOfWell());
        auto segments_copy = ws.segments;
        segments_copy.scale_pressure(bhp);
        const auto& segment_pressure = segments_copy.pressure;
        for (int seg = 0; seg < nseg; ++seg) {
            for (const int perf : this->segments_.perforations()[seg]) {
                const int cell_idx = this->well_cells_[perf];
                const auto& intQuants = ebosSimulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                // flux for each perforation
                std::vector<Scalar> mob(this->num_components_, 0.);
                getMobility(ebosSimulator, perf, mob, deferred_logger);
                const double trans_mult = ebosSimulator.problem().template wellTransMultiplier<double>(intQuants, cell_idx);
                const auto& wellstate_nupcol = ebosSimulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
                const std::vector<Scalar> Tw = this->wellIndex(perf, intQuants, trans_mult, wellstate_nupcol);
                const Scalar seg_pressure = segment_pressure[seg];
                std::vector<Scalar> cq_s(this->num_components_, 0.);
                Scalar perf_press = 0.0;
                PerforationRates perf_rates;
                computePerfRate(intQuants, mob, Tw, seg, perf, seg_pressure,
                                allow_cf, cq_s, perf_press, perf_rates, deferred_logger);

                for(int p = 0; p < np; ++p) {
                    well_flux[this->ebosCompIdxToFlowCompIdx(p)] += cq_s[p];
                }
            }
        }
        this->parallel_well_info_.communication().sum(well_flux.data(), well_flux.size());
    }


    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                      const Scalar& bhp,
                                      std::vector<double>& well_flux,
                                      DeferredLogger& deferred_logger) const
    {
        // creating a copy of the well itself, to avoid messing up the explicit information
        // during this copy, the only information not copied properly is the well controls
        MultisegmentWell<TypeTag> well_copy(*this);
        well_copy.debug_cost_counter_ = 0;

        // store a copy of the well state, we don't want to update the real well state
        WellState well_state_copy = ebosSimulator.problem().wellModel().wellState();
        const auto& group_state = ebosSimulator.problem().wellModel().groupState();
        auto& ws = well_state_copy.well(this->index_of_well_);

        // Get the current controls.
        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        auto inj_controls = well_copy.well_ecl_.isInjector()
            ? well_copy.well_ecl_.injectionControls(summary_state)
            : Well::InjectionControls(0);
        auto prod_controls = well_copy.well_ecl_.isProducer()
            ? well_copy.well_ecl_.productionControls(summary_state) :
            Well::ProductionControls(0);

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        if (well_copy.well_ecl_.isInjector()) {
            inj_controls.bhp_limit = bhp;
            ws.injection_cmode = Well::InjectorCMode::BHP;
        } else {
            prod_controls.bhp_limit = bhp;
            ws.production_cmode = Well::ProducerCMode::BHP;
        }
        ws.bhp = bhp;
        well_copy.scaleSegmentPressuresWithBhp(well_state_copy);

        // initialized the well rates with the potentials i.e. the well rates based on bhp
        const int np = this->number_of_phases_;
        bool trivial = true;
        for (int phase = 0; phase < np; ++phase){
            trivial = trivial && (ws.well_potentials[phase] == 0.0) ;
        }
        if (!trivial) {
            const double sign = well_copy.well_ecl_.isInjector() ? 1.0 : -1.0;
            for (int phase = 0; phase < np; ++phase) {
                ws.surface_rates[phase] = sign * ws.well_potentials[phase];
            }
        }
        well_copy.scaleSegmentRatesWithWellRates(this->segments_.inlets(),
                                                 this->segments_.perforations(),
                                                 well_state_copy);

        well_copy.calculateExplicitQuantities(ebosSimulator, well_state_copy, deferred_logger);
        const double dt = ebosSimulator.timeStepSize();
        // iterate to get a solution at the given bhp.
        well_copy.iterateWellEqWithControl(ebosSimulator, dt, inj_controls, prod_controls, well_state_copy, group_state,
                                           deferred_logger);

        // compute the potential and store in the flux vector.
        well_flux.clear();
        well_flux.resize(np, 0.0);
        for (int compIdx = 0; compIdx < this->num_components_; ++compIdx) {
            const EvalWell rate = well_copy.primary_variables_.getQs(compIdx);
            well_flux[this->ebosCompIdxToFlowCompIdx(compIdx)] = rate.value();
        }
        debug_cost_counter_ += well_copy.debug_cost_counter_;
    }



    template<typename TypeTag>
    std::vector<double>
    MultisegmentWell<TypeTag>::
    computeWellPotentialWithTHP(
          const WellState& well_state,
          const Simulator& ebos_simulator,
          DeferredLogger& deferred_logger) const
    {
        std::vector<double> potentials(this->number_of_phases_, 0.0);
        const auto& summary_state = ebos_simulator.vanguard().summaryState();

        const auto& well = this->well_ecl_;
        if (well.isInjector()){
            auto bhp_at_thp_limit = computeBhpAtThpLimitInj(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well.injectionControls(summary_state);
                const double bhp = std::min(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhpIterations(ebos_simulator, bhp, potentials, deferred_logger);
                deferred_logger.debug("Converged thp based potential calculation for well "
                                      + this->name() + ", at bhp = " + std::to_string(bhp));
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + this->name() + ". Instead the bhp based value is used");
                const auto& controls = well.injectionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhpIterations(ebos_simulator, bhp, potentials, deferred_logger);
            }
        } else {
            auto bhp_at_thp_limit = computeBhpAtThpLimitProd(
                  well_state, ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well.productionControls(summary_state);
                const double bhp = std::max(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhpIterations(ebos_simulator, bhp, potentials, deferred_logger);
                deferred_logger.debug("Converged thp based potential calculation for well "
                                      + this->name() + ", at bhp = " + std::to_string(bhp));
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + this->name() + ". Instead the bhp based value is used");
                const auto& controls = well.productionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhpIterations(ebos_simulator, bhp, potentials, deferred_logger);
            }
        }

        return potentials;
    }

    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    computeWellPotentialsImplicit(const Simulator& ebos_simulator,
                                  std::vector<double>& well_potentials,
                                  DeferredLogger& deferred_logger) const
    {
        // Create a copy of the well.
        // TODO: check if we can avoid taking multiple copies. Call from updateWellPotentials 
        // is allready a copy, but not from other calls. 
        MultisegmentWell<TypeTag> well_copy(*this);
        well_copy.debug_cost_counter_ = 0;

        // store a copy of the well state, we don't want to update the real well state
        WellState well_state_copy = ebos_simulator.problem().wellModel().wellState();
        const auto& group_state = ebos_simulator.problem().wellModel().groupState();
        auto& ws = well_state_copy.well(this->index_of_well_);
        
        // get current controls        
        const auto& summary_state = ebos_simulator.vanguard().summaryState();
        auto inj_controls = well_copy.well_ecl_.isInjector()
            ? well_copy.well_ecl_.injectionControls(summary_state)
            : Well::InjectionControls(0);
        auto prod_controls = well_copy.well_ecl_.isProducer()
            ? well_copy.well_ecl_.productionControls(summary_state)
            : Well::ProductionControls(0);
        
        // prepare/modify well state and control
        well_copy.prepareForPotentialCalculations(summary_state, well_state_copy, inj_controls, prod_controls);
        
        well_copy.scaleSegmentPressuresWithBhp(well_state_copy);

        // initialize rates from previous potentials
        const int np = this->number_of_phases_;
        bool trivial = true;
        for (int phase = 0; phase < np; ++phase){
            trivial = trivial && (ws.well_potentials[phase] == 0.0) ;
        }
        if (!trivial) {
            const double sign = well_copy.well_ecl_.isInjector() ? 1.0 : -1.0;
            for (int phase = 0; phase < np; ++phase) {
                ws.surface_rates[phase] = sign * ws.well_potentials[phase];
            }
        }
        well_copy.scaleSegmentRatesWithWellRates(this->segments_.inlets(),
                                                 this->segments_.perforations(),
                                                 well_state_copy);

        well_copy.calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);
        const double dt = ebos_simulator.timeStepSize();
        // solve equations
        bool converged = false;
        if (this->well_ecl_.isProducer() && this->wellHasTHPConstraints(summary_state)) {
            converged = well_copy.solveWellWithTHPConstraint(ebos_simulator, dt, inj_controls, prod_controls, well_state_copy, group_state, deferred_logger);
        } else {
            converged = well_copy.iterateWellEqWithSwitching(ebos_simulator, dt, inj_controls, prod_controls, well_state_copy, group_state, deferred_logger);
        }

        // fetch potentials (sign is updated on the outside).
        well_potentials.clear();
        well_potentials.resize(np, 0.0);
        for (int compIdx = 0; compIdx < this->num_components_; ++compIdx) {
            const EvalWell rate = well_copy.primary_variables_.getQs(compIdx);
            well_potentials[this->ebosCompIdxToFlowCompIdx(compIdx)] = rate.value();
        }
        debug_cost_counter_ += well_copy.debug_cost_counter_;
        return converged;
    }

    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    solveEqAndUpdateWellState(const SummaryState& summary_state,
                              WellState& well_state,
                              DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        try{
            const BVectorWell dx_well = this->linSys_.solve();

            updateWellState(summary_state, dx_well, well_state, deferred_logger);
        }
        catch(const NumericalProblem& exp) {
            // Add information about the well and log to deferred logger
            // (Logging done inside of solve() method will only be seen if
            // this is the process with rank zero)
            deferred_logger.problem("In MultisegmentWell::solveEqAndUpdateWellState for well "
                                    + this->name() +": "+exp.what());
            throw;
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computePerfCellPressDiffs(const Simulator& ebosSimulator)
    {
        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {

            std::vector<double> kr(this->number_of_phases_, 0.0);
            std::vector<double> density(this->number_of_phases_, 0.0);

            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = ebosSimulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
            const auto& fs = intQuants.fluidState();

            double sum_kr = 0.;

            const PhaseUsage& pu = this->phaseUsage();
            if (pu.phase_used[Water]) {
                const int water_pos = pu.phase_pos[Water];
                kr[water_pos] = intQuants.relativePermeability(FluidSystem::waterPhaseIdx).value();
                sum_kr += kr[water_pos];
                density[water_pos] = fs.density(FluidSystem::waterPhaseIdx).value();
            }

            if (pu.phase_used[Oil]) {
                const int oil_pos = pu.phase_pos[Oil];
                kr[oil_pos] = intQuants.relativePermeability(FluidSystem::oilPhaseIdx).value();
                sum_kr += kr[oil_pos];
                density[oil_pos] = fs.density(FluidSystem::oilPhaseIdx).value();
            }

            if (pu.phase_used[Gas]) {
                const int gas_pos = pu.phase_pos[Gas];
                kr[gas_pos] = intQuants.relativePermeability(FluidSystem::gasPhaseIdx).value();
                sum_kr += kr[gas_pos];
                density[gas_pos] = fs.density(FluidSystem::gasPhaseIdx).value();
            }

            assert(sum_kr != 0.);

            // calculate the average density
            double average_density = 0.;
            for (int p = 0; p < this->number_of_phases_; ++p) {
                average_density += kr[p] * density[p];
            }
            average_density /= sum_kr;

            this->cell_perforation_pressure_diffs_[perf] = this->gravity_ * average_density * this->cell_perforation_depth_diffs_[perf];
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeInitialSegmentFluids(const Simulator& ebos_simulator)
    {
        for (int seg = 0; seg < this->numberOfSegments(); ++seg) {
            // TODO: trying to reduce the times for the surfaceVolumeFraction calculation
            const double surface_volume = getSegmentSurfaceVolume(ebos_simulator, seg).value();
            for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                segment_fluid_initial_[seg][comp_idx] = surface_volume * this->primary_variables_.surfaceVolumeFraction(seg, comp_idx).value();
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellState(const SummaryState& summary_state,
                    const BVectorWell& dwells,
                    WellState& well_state,
                    DeferredLogger& deferred_logger,
                    const double relaxation_factor)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        const double dFLimit = this->param_.dwell_fraction_max_;
        const double max_pressure_change = this->param_.max_pressure_change_ms_wells_;
        const bool stop_or_zero_rate_target = this->stopppedOrZeroRateTarget(summary_state, well_state);
        this->primary_variables_.updateNewton(dwells,
                                              relaxation_factor,
                                              dFLimit,
                                              stop_or_zero_rate_target,
                                              max_pressure_change);

        this->primary_variables_.copyToWellState(*this, getRefDensity(), stop_or_zero_rate_target,
                                                 well_state, summary_state, deferred_logger);

        {
            auto& ws = well_state.well(this->index_of_well_);
            this->segments_.copyPhaseDensities(ws.pu, ws.segments);
        }

        Base::calculateReservoirRates(well_state.well(this->index_of_well_));
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state,
                                DeferredLogger& deferred_logger)
    {
        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        updatePrimaryVariables(summary_state, well_state, deferred_logger);
        initPrimaryVariablesEvaluation();
        computePerfCellPressDiffs(ebosSimulator);
        computeInitialSegmentFluids(ebosSimulator);
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateProductivityIndex(const Simulator& ebosSimulator,
                            const WellProdIndexCalculator& wellPICalc,
                            WellState& well_state,
                            DeferredLogger& deferred_logger) const
    {
        auto fluidState = [&ebosSimulator, this](const int perf)
        {
            const auto cell_idx = this->well_cells_[perf];
            return ebosSimulator.model()
               .intensiveQuantities(cell_idx, /*timeIdx=*/ 0).fluidState();
        };

        const int np = this->number_of_phases_;
        auto setToZero = [np](double* x) -> void
        {
            std::fill_n(x, np, 0.0);
        };

        auto addVector = [np](const double* src, double* dest) -> void
        {
            std::transform(src, src + np, dest, dest, std::plus<>{});
        };

        auto& ws = well_state.well(this->index_of_well_);
        auto& perf_data = ws.perf_data;
        auto* connPI = perf_data.prod_index.data();
        auto* wellPI = ws.productivity_index.data();

        setToZero(wellPI);

        const auto preferred_phase = this->well_ecl_.getPreferredPhase();
        auto subsetPerfID   = 0;

        for ( const auto& perf : *this->perf_data_){
            auto allPerfID = perf.ecl_index;

            auto connPICalc = [&wellPICalc, allPerfID](const double mobility) -> double
            {
                return wellPICalc.connectionProdIndStandard(allPerfID, mobility);
            };

            std::vector<Scalar> mob(this->num_components_, 0.0);
            getMobility(ebosSimulator, static_cast<int>(subsetPerfID), mob, deferred_logger);

            const auto& fs = fluidState(subsetPerfID);
            setToZero(connPI);

            if (this->isInjector()) {
                this->computeConnLevelInjInd(fs, preferred_phase, connPICalc,
                                             mob, connPI, deferred_logger);
            }
            else {  // Production or zero flow rate
                this->computeConnLevelProdInd(fs, connPICalc, mob, connPI);
            }

            addVector(connPI, wellPI);

            ++subsetPerfID;
            connPI += np;
        }

        assert (static_cast<int>(subsetPerfID) == this->number_of_perforations_ &&
                "Internal logic error in processing connections for PI/II");
    }





    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    connectionDensity(const int globalConnIdx,
                      [[maybe_unused]] const int openConnIdx) const
    {
        // Simple approximation: Mixture density at reservoir connection is
        // mixture density at connection's segment.

        const auto segNum = this->wellEcl()
            .getConnections()[globalConnIdx].segment();

        const auto segIdx = this->wellEcl()
            .getSegments().segmentNumberToIndex(segNum);

        return this->segments_.density(segIdx).value();
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    addWellContributions(SparseMatrixAdapter& jacobian) const
    {
        this->linSys_.extract(jacobian);
    }


   template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    addWellPressureEquations(PressureMatrix& jacobian,
                             const BVector& weights,
                             const int pressureVarIndex,
                             const bool use_well_weights,
                             const WellState& well_state) const
    {
        // Add the pressure contribution to the cpr system for the well
        this->linSys_.extractCPRPressureMatrix(jacobian,
                                               weights,
                                               pressureVarIndex,
                                               use_well_weights,
                                               *this,
                                               this->SPres,
                                               well_state);
    }


    template<typename TypeTag>
    template<class Value>
    void
    MultisegmentWell<TypeTag>::
    computePerfRate(const Value& pressure_cell,
                    const Value& rs,
                    const Value& rv,
                    const std::vector<Value>& b_perfcells,
                    const std::vector<Value>& mob_perfcells,
                    const std::vector<Scalar>& Tw,
                    const int perf,
                    const Value& segment_pressure,
                    const Value& segment_density,
                    const bool& allow_cf,
                    const std::vector<Value>& cmix_s,
                    std::vector<Value>& cq_s,
                    Value& perf_press,
                    PerforationRates& perf_rates,
                    DeferredLogger& deferred_logger) const
    {
        // pressure difference between the segment and the perforation
        const Value perf_seg_press_diff = this->gravity() * segment_density *
                                          this->segments_.perforation_depth_diff(perf);
        // pressure difference between the perforation and the grid cell
        const double cell_perf_press_diff = this->cell_perforation_pressure_diffs_[perf];

        // perforation pressure is the wellbore pressure corrected to perforation depth
        // (positive sign due to convention in segments_.perforation_depth_diff() )
        perf_press = segment_pressure + perf_seg_press_diff;

        // cell pressure corrected to perforation depth
        const Value cell_press_at_perf = pressure_cell - cell_perf_press_diff;

        // Pressure drawdown (also used to determine direction of flow)
        const Value drawdown = cell_press_at_perf - perf_press;

        // producing perforations
        if (drawdown > 0.0) {
            // Do nothing if crossflow is not allowed
            if (!allow_cf && this->isInjector()) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int comp_idx = 0; comp_idx < this->numComponents(); ++comp_idx) {
                const Value cq_p = - Tw[comp_idx] * (mob_perfcells[comp_idx] * drawdown);
                cq_s[comp_idx] = b_perfcells[comp_idx] * cq_p;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const Value cq_s_oil = cq_s[oilCompIdx];
                const Value cq_s_gas = cq_s[gasCompIdx];
                cq_s[gasCompIdx] += rs * cq_s_oil;
                cq_s[oilCompIdx] += rv * cq_s_gas;
            }
        } else { // injecting perforations
            // Do nothing if crossflow is not allowed
            if (!allow_cf && this->isProducer()) {
                return;
            }

            // for injecting perforations, we use total mobility
            Value total_mob = mob_perfcells[0];
            for (int comp_idx = 1; comp_idx < this->numComponents(); ++comp_idx) {
                total_mob += mob_perfcells[comp_idx];
            }

            // compute volume ratio between connection and at standard conditions
            Value volume_ratio = 0.0;
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                volume_ratio += cmix_s[waterCompIdx] / b_perfcells[waterCompIdx];
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

                // Incorporate RS/RV factors if both oil and gas active
                // TODO: not sure we use rs rv from the perforation cells when handling injecting perforations
                // basically, for injecting perforations, the wellbore is the upstreaming side.
                const Value d = 1.0 - rv * rs;

                if (getValue(d) == 0.0) {
                    OPM_DEFLOG_PROBLEM(NumericalProblem,
                                       fmt::format("Zero d value obtained for well {} "
                                                   "during flux calculation with rs {} and rv {}",
                                                   this->name(), rs, rv),
                                       deferred_logger);
                }

                const Value tmp_oil = (cmix_s[oilCompIdx] - rv * cmix_s[gasCompIdx]) / d;
                volume_ratio += tmp_oil / b_perfcells[oilCompIdx];

                const Value tmp_gas = (cmix_s[gasCompIdx] - rs * cmix_s[oilCompIdx]) / d;
                volume_ratio += tmp_gas / b_perfcells[gasCompIdx];
            } else { // not having gas and oil at the same time
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    volume_ratio += cmix_s[oilCompIdx] / b_perfcells[oilCompIdx];
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    volume_ratio += cmix_s[gasCompIdx] / b_perfcells[gasCompIdx];
                }
            }
            // injecting connections total volumerates at standard conditions
            for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
                const Value cqt_i = - Tw[componentIdx] * (total_mob * drawdown);
                Value cqt_is = cqt_i / volume_ratio;
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is;
            }
        } // end for injection perforations

        // calculating the perforation solution gas rate and solution oil rates
        if (this->isProducer()) {
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                // TODO: the formulations here remain to be tested with cases with strong crossflow through production wells
                // s means standard condition, r means reservoir condition
                // q_os = q_or * b_o + rv * q_gr * b_g
                // q_gs = q_gr * g_g + rs * q_or * b_o
                // d = 1.0 - rs * rv
                // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)

                const double d = 1.0 - getValue(rv) * getValue(rs);
                // vaporized oil into gas
                // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
                perf_rates.vap_oil = getValue(rv) * (getValue(cq_s[gasCompIdx]) - getValue(rs) * getValue(cq_s[oilCompIdx])) / d;
                // dissolved of gas in oil
                // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
                perf_rates.dis_gas = getValue(rs) * (getValue(cq_s[oilCompIdx]) - getValue(rv) * getValue(cq_s[gasCompIdx])) / d;
            }
        }
    }

    template <typename TypeTag>
    template<class Value>
    void
    MultisegmentWell<TypeTag>::
    computePerfRate(const IntensiveQuantities& int_quants,
                    const std::vector<Value>& mob_perfcells,
                    const std::vector<Scalar>& Tw,
                    const int seg,
                    const int perf,
                    const Value& segment_pressure,
                    const bool& allow_cf,
                    std::vector<Value>& cq_s,
                    Value& perf_press,
                    PerforationRates& perf_rates,
                    DeferredLogger& deferred_logger) const

    {
        auto obtain = [this](const Eval& value)
                      {
                          if constexpr (std::is_same_v<Value, Scalar>) {
                              static_cast<void>(this); // suppress clang warning
                              return getValue(value);
                          } else {
                              return this->extendEval(value);
                          }
                      };
        auto obtainN = [](const auto& value)
        {
            if constexpr (std::is_same_v<Value, Scalar>) {
                return getValue(value);
            } else {
                return value;
            }
        };
        const auto& fs = int_quants.fluidState();

        const Value pressure_cell = obtain(this->getPerfCellPressure(fs));
        const Value rs = obtain(fs.Rs());
        const Value rv = obtain(fs.Rv());

        // not using number_of_phases_ because of solvent
        std::vector<Value> b_perfcells(this->num_components_, 0.0);

        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells[compIdx] = obtain(fs.invB(phaseIdx));
        }

        std::vector<Value> cmix_s(this->numComponents(), 0.0);
        for (int comp_idx = 0; comp_idx < this->numComponents(); ++comp_idx) {
            cmix_s[comp_idx] = obtainN(this->primary_variables_.surfaceVolumeFraction(seg, comp_idx));
        }

        this->computePerfRate(pressure_cell,
                              rs,
                              rv,
                              b_perfcells,
                              mob_perfcells,
                              Tw,
                              perf,
                              segment_pressure,
                              obtainN(this->segments_.density(seg)),
                              allow_cf,
                              cmix_s,
                              cq_s,
                              perf_press,
                              perf_rates,
                              deferred_logger);
    }

    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeSegmentFluidProperties(const Simulator& ebosSimulator, DeferredLogger& deferred_logger)
    {
        // TODO: the concept of phases and components are rather confusing in this function.
        // needs to be addressed sooner or later.

        // get the temperature for later use. It is only useful when we are not handling
        // thermal related simulation
        // basically, it is a single value for all the segments

        EvalWell temperature;
        EvalWell saltConcentration;
        // not sure how to handle the pvt region related to segment
        // for the current approach, we use the pvt region of the first perforated cell
        // although there are some text indicating using the pvt region of the lowest
        // perforated cell
        // TODO: later to investigate how to handle the pvt region
        int pvt_region_index;
        {
            // using the first perforated cell
            const int cell_idx = this->well_cells_[0];
            const auto& intQuants = ebosSimulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            temperature.setValue(fs.temperature(FluidSystem::oilPhaseIdx).value());
            saltConcentration = this->extendEval(fs.saltConcentration());
            pvt_region_index = fs.pvtRegionIndex();
        }

        this->segments_.computeFluidProperties(temperature,
                                               saltConcentration,
                                               this->primary_variables_,
                                               pvt_region_index,
                                               deferred_logger);
    }

    template <typename TypeTag>
    template<class Value>
    void
    MultisegmentWell<TypeTag>::
    getMobility(const Simulator& ebosSimulator,
                const int perf,
                std::vector<Value>& mob,
                DeferredLogger& deferred_logger) const
    {
        auto obtain = [this](const Eval& value)
                      {
                          if constexpr (std::is_same_v<Value, Scalar>) {
                              static_cast<void>(this); // suppress clang warning
                              return getValue(value);
                          } else {
                              return this->extendEval(value);
                          }
                      };

        WellInterface<TypeTag>::getMobility(ebosSimulator, perf, mob, obtain, deferred_logger);

        if (this->isInjector() && this->well_ecl_.getInjMultMode() != Well::InjMultMode::NONE) {
            const auto perf_ecl_index = this->perforationData()[perf].ecl_index;
            const Connection& con = this->well_ecl_.getConnections()[perf_ecl_index];
            const int seg = this->segmentNumberToIndex(con.segment());
            // from the reference results, it looks like MSW uses segment pressure instead of BHP here
            // Note: this is against the documented definition.
            // we can change this depending on what we want
            const double segment_pres = this->primary_variables_.getSegmentPressure(seg).value();
            const double perf_seg_press_diff = this->gravity() * this->segments_.density(seg).value()
                                                               * this->segments_.perforation_depth_diff(perf);
            const double perf_press = segment_pres + perf_seg_press_diff;
            const double multiplier = this->getInjMult(perf, segment_pres, perf_press);
            for (std::size_t i = 0; i < mob.size(); ++i) {
                mob[i] *= multiplier;
            }
        }
    }



    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    getRefDensity() const
    {
        return this->segments_.getRefDensity();
    }

    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkOperabilityUnderBHPLimit(const WellState& /*well_state*/, const Simulator& ebos_simulator, DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const double bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints(summaryState) ) {
            // if the BHP limit is not defaulted or the well does not have a THP limit
            // we need to check the BHP limit
            double total_ipr_mass_rate = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                const double ipr_rate = this->ipr_a_[compIdx] - this->ipr_b_[compIdx] * bhp_limit;

                const double rho = FluidSystem::referenceDensity( phaseIdx, Base::pvtRegionIdx() );
                total_ipr_mass_rate += ipr_rate * rho;
            }
            if ( (this->isProducer() && total_ipr_mass_rate < 0.) || (this->isInjector() && total_ipr_mass_rate > 0.) ) {
                this->operability_status_.operable_under_only_bhp_limit = false;
            }

            // checking whether running under BHP limit will violate THP limit
            if (this->operability_status_.operable_under_only_bhp_limit && this->wellHasTHPConstraints(summaryState)) {
                // option 1: calculate well rates based on the BHP limit.
                // option 2: stick with the above IPR curve
                // we use IPR here
                std::vector<double> well_rates_bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                const double thp_limit = this->getTHPConstraint(summaryState);
                const double thp = WellBhpThpCalculator(*this).calculateThpFromBhp(well_rates_bhp_limit,
                                                                                   bhp_limit,
                                                                                   this->getRefDensity(),
                                                                                   this->wellEcl().alq_value(summaryState),
                                                                                   thp_limit,
                                                                                   deferred_logger);
                if ( (this->isProducer() && thp < thp_limit) || (this->isInjector() && thp > thp_limit) ) {
                    this->operability_status_.obey_thp_limit_under_bhp_limit = false;
                }
            }
        } else {
            // defaulted BHP and there is a THP constraint
            // default BHP limit is about 1 atm.
            // when applied the hydrostatic pressure correction dp,
            // most likely we get a negative value (bhp + dp)to search in the VFP table,
            // which is not desirable.
            // we assume we can operate under defaulted BHP limit and will violate the THP limit
            // when operating under defaulted BHP limit.
            this->operability_status_.operable_under_only_bhp_limit = true;
            this->operability_status_.obey_thp_limit_under_bhp_limit = false;
        }
    }



    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const
    {
        // TODO: not handling solvent related here for now

        // initialize all the values to be zero to begin with
        std::fill(this->ipr_a_.begin(), this->ipr_a_.end(), 0.);
        std::fill(this->ipr_b_.begin(), this->ipr_b_.end(), 0.);

        const int nseg = this->numberOfSegments();
        std::vector<double> seg_dp(nseg, 0.0);
        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the perforation rate for each perforation that belongs to this segment
            const double dp = this->getSegmentDp(seg,
                                                 this->segments_.density(seg).value(),
                                                 seg_dp);
            seg_dp[seg] = dp;
            for (const int perf : this->segments_.perforations()[seg]) {
                std::vector<Scalar> mob(this->num_components_, 0.0);

                // TODO: maybe we should store the mobility somewhere, so that we only need to calculate it one per iteration
                getMobility(ebos_simulator, perf, mob, deferred_logger);

                const int cell_idx = this->well_cells_[perf];
                const auto& int_quantities = ebos_simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                const auto& fs = int_quantities.fluidState();
                // pressure difference between the segment and the perforation
                const double perf_seg_press_diff = this->segments_.getPressureDiffSegPerf(seg, perf);
                // pressure difference between the perforation and the grid cell
                const double cell_perf_press_diff = this->cell_perforation_pressure_diffs_[perf];
                const double pressure_cell = this->getPerfCellPressure(fs).value();

                // calculating the b for the connection
                std::vector<double> b_perf(this->num_components_);
                for (std::size_t phase = 0; phase < FluidSystem::numPhases; ++phase) {
                    if (!FluidSystem::phaseIsActive(phase)) {
                        continue;
                    }
                    const unsigned comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phase));
                    b_perf[comp_idx] = fs.invB(phase).value();
                }

                // the pressure difference between the connection and BHP
                const double h_perf = cell_perf_press_diff + perf_seg_press_diff + dp;
                const double pressure_diff = pressure_cell - h_perf;

                // do not take into consideration the crossflow here.
                if ( (this->isProducer() && pressure_diff < 0.) || (this->isInjector() && pressure_diff > 0.) ) {
                    deferred_logger.debug("CROSSFLOW_IPR",
                                    "cross flow found when updateIPR for well " + this->name());
                }

                // the well index associated with the connection
                const double trans_mult = ebos_simulator.problem().template wellTransMultiplier<double>(int_quantities, cell_idx);
                const auto& wellstate_nupcol = ebos_simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
                const std::vector<Scalar> tw_perf = this->wellIndex(perf, int_quantities, trans_mult, wellstate_nupcol);  
                std::vector<double> ipr_a_perf(this->ipr_a_.size());
                std::vector<double> ipr_b_perf(this->ipr_b_.size());
                for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                    const double tw_mob = tw_perf[comp_idx] * mob[comp_idx] * b_perf[comp_idx];
                    ipr_a_perf[comp_idx] += tw_mob * pressure_diff;
                    ipr_b_perf[comp_idx] += tw_mob;
                }

                // we need to handle the rs and rv when both oil and gas are present
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned oil_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    const unsigned gas_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    const double rs = (fs.Rs()).value();
                    const double rv = (fs.Rv()).value();

                    const double dis_gas_a = rs * ipr_a_perf[oil_comp_idx];
                    const double vap_oil_a = rv * ipr_a_perf[gas_comp_idx];

                    ipr_a_perf[gas_comp_idx] += dis_gas_a;
                    ipr_a_perf[oil_comp_idx] += vap_oil_a;

                    const double dis_gas_b = rs * ipr_b_perf[oil_comp_idx];
                    const double vap_oil_b = rv * ipr_b_perf[gas_comp_idx];

                    ipr_b_perf[gas_comp_idx] += dis_gas_b;
                    ipr_b_perf[oil_comp_idx] += vap_oil_b;
                }

                for (std::size_t comp_idx = 0; comp_idx < ipr_a_perf.size(); ++comp_idx) {
                    this->ipr_a_[comp_idx] += ipr_a_perf[comp_idx];
                    this->ipr_b_[comp_idx] += ipr_b_perf[comp_idx];
                }
            }
        }
    }

    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateIPRImplicit(const Simulator& ebos_simulator, WellState& well_state, DeferredLogger& deferred_logger)
    {
        // Compute IPR based on *converged* well-equation:
        // For a component rate r the derivative dr/dbhp is obtained by 
        // dr/dbhp = - (partial r/partial x) * inv(partial Eq/partial x) * (partial Eq/partial bhp_target)
        // where Eq(x)=0 is the well equation setup with bhp control and primary variables x 

        // We shouldn't have zero rates at this stage, but check
        bool zero_rates;
        auto rates = well_state.well(this->index_of_well_).surface_rates;
        zero_rates = true;
        for (std::size_t p = 0; p < rates.size(); ++p) {
            zero_rates &= rates[p] == 0.0;
        }
        auto& ws = well_state.well(this->index_of_well_);
        if (zero_rates) {
            const auto msg = fmt::format("updateIPRImplicit: Well {} has zero rate, IPRs might be problematic", this->name());
            deferred_logger.debug(msg);
            /*
            // could revert to standard approach here:    
            updateIPR(ebos_simulator, deferred_logger);
            for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx){
                const int idx = this->ebosCompIdxToFlowCompIdx(comp_idx);
                ws.implicit_ipr_a[idx] = this->ipr_a_[comp_idx];
                ws.implicit_ipr_b[idx] = this->ipr_b_[comp_idx];
            }
            return;
            */
        }
        const auto& group_state  = ebos_simulator.problem().wellModel().groupState();

        std::fill(ws.implicit_ipr_a.begin(), ws.implicit_ipr_a.end(), 0.);
        std::fill(ws.implicit_ipr_b.begin(), ws.implicit_ipr_b.end(), 0.);
        //WellState well_state_copy = well_state;    
        auto inj_controls = Well::InjectionControls(0);
        auto prod_controls = Well::ProductionControls(0);
        prod_controls.addControl(Well::ProducerCMode::BHP);
        prod_controls.bhp_limit = well_state.well(this->index_of_well_).bhp;

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        const auto cmode = ws.production_cmode;
        ws.production_cmode = Well::ProducerCMode::BHP;
        const double dt = ebos_simulator.timeStepSize();
        assembleWellEqWithoutIteration(ebos_simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

        BVectorWell rhs(this->numberOfSegments());
        rhs = 0.0;
        rhs[0][SPres] = -1.0;

        const BVectorWell x_well = this->linSys_.solve(rhs);
        constexpr int num_eq = MSWEval::numWellEq;
        for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx){
            const EvalWell comp_rate = this->primary_variables_.getQs(comp_idx);
            const int idx = this->ebosCompIdxToFlowCompIdx(comp_idx);
            for (size_t pvIdx = 0; pvIdx < num_eq; ++pvIdx) {
                // well primary variable derivatives in EvalWell start at position Indices::numEq 
                ws.implicit_ipr_b[idx] -= x_well[0][pvIdx]*comp_rate.derivative(pvIdx+Indices::numEq);
            }
            ws.implicit_ipr_a[idx] = ws.implicit_ipr_b[idx]*ws.bhp - comp_rate.value();
        }
        // reset cmode
        ws.production_cmode = cmode;
    }    

    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkOperabilityUnderTHPLimit(
             const Simulator& ebos_simulator,
             const WellState& well_state,
             DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto obtain_bhp = this->isProducer()
            ? computeBhpAtThpLimitProd(
                        well_state, ebos_simulator, summaryState, deferred_logger)
            : computeBhpAtThpLimitInj(ebos_simulator, summaryState, deferred_logger);

        if (obtain_bhp) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const double  bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
            this->operability_status_.obey_bhp_limit_with_thp_limit = (*obtain_bhp >= bhp_limit);

            const double thp_limit = this->getTHPConstraint(summaryState);
            if (this->isProducer() && *obtain_bhp < thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is SMALLER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a producer for well " + this->name();
                deferred_logger.debug(msg);
            }
            else if (this->isInjector() && *obtain_bhp > thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is LARGER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a injector for well " + this->name();
                deferred_logger.debug(msg);
            }
        } else {
            // Shutting wells that can not find bhp value from thp
            // when under THP control
            this->operability_status_.can_obtain_bhp_with_thp_limit = false;
            this->operability_status_.obey_bhp_limit_with_thp_limit = false;
            if (!this->wellIsStopped()) {
                const double thp_limit = this->getTHPConstraint(summaryState);
                deferred_logger.debug(" could not find bhp value at thp limit "
                                      + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                      + " bar for well " + this->name() + ", the well might need to be closed ");
            }
        }
    }





    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    iterateWellEqWithControl(const Simulator& ebosSimulator,
                             const double dt,
                             const Well::InjectionControls& inj_controls,
                             const Well::ProductionControls& prod_controls,
                             WellState& well_state,
                             const GroupState& group_state,
                             DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return true;

        const int max_iter_number = this->param_.max_inner_iter_ms_wells_;

        {
            // getWellFiniteResiduals returns false for nan/inf residuals
            const auto& [isFinite, residuals] = this->getFiniteWellResiduals(Base::B_avg_, deferred_logger);
            if(!isFinite)
                return false;
        }

        std::vector<std::vector<Scalar> > residual_history;
        std::vector<double> measure_history;
        int it = 0;
        // relaxation factor
        double relaxation_factor = 1.;
        const double min_relaxation_factor = 0.6;
        bool converged = false;
        int stagnate_count = 0;
        bool relax_convergence = false;
        this->regularize_ = false;
        for (; it < max_iter_number; ++it, ++debug_cost_counter_) {

            assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            BVectorWell dx_well;
            try{
                dx_well = this->linSys_.solve();
            }
            catch(const NumericalProblem& exp) {
                // Add information about the well and log to deferred logger
                // (Logging done inside of solve() method will only be seen if
                // this is the process with rank zero)
                deferred_logger.problem("In MultisegmentWell::iterateWellEqWithControl for well "
                                        + this->name() +": "+exp.what());
                throw;
            }

            if (it > this->param_.strict_inner_iter_wells_) {
                relax_convergence = true;
                this->regularize_ = true;
            }

            const auto& summary_state = ebosSimulator.vanguard().summaryState();
            const auto report = getWellConvergence(summary_state, well_state, Base::B_avg_, deferred_logger, relax_convergence);
            if (report.converged()) {
                converged = true;
                break;
            }

            {
                // getFinteWellResiduals returns false for nan/inf residuals
                const auto& [isFinite, residuals] = this->getFiniteWellResiduals(Base::B_avg_, deferred_logger);
                if (!isFinite)
                    return false;

                residual_history.push_back(residuals);
                measure_history.push_back(this->getResidualMeasureValue(well_state,
                                                                    residual_history[it],
                                                                    this->param_.tolerance_wells_,
                                                                    this->param_.tolerance_pressure_ms_wells_,
                                                                    deferred_logger) );
            }


            bool is_oscillate = false;
            bool is_stagnate = false;

            this->detectOscillations(measure_history, is_oscillate, is_stagnate);
            // TODO: maybe we should have more sophisticated strategy to recover the relaxation factor,
            // for example, to recover it to be bigger

            if (is_oscillate || is_stagnate) {
                // HACK!
                std::ostringstream sstr;
                if (relaxation_factor == min_relaxation_factor) {
                    // Still stagnating, terminate iterations if 5 iterations pass.
                    ++stagnate_count;
                    if (stagnate_count == 6) {
                        sstr << " well " << this->name() << " observes severe stagnation and/or oscillation. We relax the tolerance and check for convergence. \n";
                        const auto reportStag = getWellConvergence(summary_state, well_state, Base::B_avg_, deferred_logger, true);
                        if (reportStag.converged()) {
                            converged = true;
                            sstr << " well " << this->name() << " manages to get converged with relaxed tolerances in " << it << " inner iterations";
                            deferred_logger.debug(sstr.str());
                            return converged;
                        }
                    }
                }

                // a factor value to reduce the relaxation_factor
                const double reduction_mutliplier = 0.9;
                relaxation_factor = std::max(relaxation_factor * reduction_mutliplier, min_relaxation_factor);

                // debug output
                if (is_stagnate) {
                    sstr << " well " << this->name() << " observes stagnation in inner iteration " << it << "\n";

                }
                if (is_oscillate) {
                    sstr << " well " << this->name() << " observes oscillation in inner iteration " << it << "\n";
                }
                sstr << " relaxation_factor is " << relaxation_factor << " now\n";

                this->regularize_ = true;
                deferred_logger.debug(sstr.str());
            }
            updateWellState(summary_state, dx_well, well_state, deferred_logger, relaxation_factor);
            initPrimaryVariablesEvaluation();
        }

        // TODO: we should decide whether to keep the updated well_state, or recover to use the old well_state
        if (converged) {
            std::ostringstream sstr;
            sstr << "     Well " << this->name() << " converged in " << it << " inner iterations.";
            if (relax_convergence)
                sstr << "      (A relaxed tolerance was used after "<< this->param_.strict_inner_iter_wells_ << " iterations)";
            deferred_logger.debug(sstr.str());
        } else {
            std::ostringstream sstr;
            sstr << "     Well " << this->name() << " did not converge in " << it << " inner iterations.";
#define EXTRA_DEBUG_MSW 0
#if EXTRA_DEBUG_MSW
            sstr << "***** Outputting the residual history for well " << this->name() << " during inner iterations:";
            for (int i = 0; i < it; ++i) {
                const auto& residual = residual_history[i];
                sstr << " residual at " << i << "th iteration ";
                for (const auto& res : residual) {
                    sstr << " " << res;
                }
                sstr << " " << measure_history[i] << " \n";
            }
#endif
#undef EXTRA_DEBUG_MSW
            deferred_logger.debug(sstr.str());
        }

        return converged;
    }


    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    iterateWellEqWithSwitching(const Simulator& ebosSimulator,
                             const double dt,
                             const Well::InjectionControls& inj_controls,
                             const Well::ProductionControls& prod_controls,
                             WellState& well_state,
                             const GroupState& group_state,
                             DeferredLogger& deferred_logger, 
                             const bool fixed_control /*false*/, 
                             const bool fixed_status /*false*/)
    {
        const int max_iter_number = this->param_.max_inner_iter_ms_wells_;

        {
            // getWellFiniteResiduals returns false for nan/inf residuals
            const auto& [isFinite, residuals] = this->getFiniteWellResiduals(Base::B_avg_, deferred_logger);
            if(!isFinite)
                return false;
        }

        std::vector<std::vector<Scalar> > residual_history;
        std::vector<double> measure_history;
        int it = 0;
        // relaxation factor
        double relaxation_factor = 1.;
        const double min_relaxation_factor = 0.6;
        bool converged = false;
        [[maybe_unused]] int stagnate_count = 0;
        bool relax_convergence = false;
        this->regularize_ = false;
        const auto& summary_state = ebosSimulator.vanguard().summaryState();

        // Always take a few (more than one) iterations after a switch before allowing a new switch
        // The optimal number here is subject to further investigation, but it has been observerved 
        // that unless this number is >1, we may get stuck in a cycle 
        const int min_its_after_switch = 3;
        int its_since_last_switch = min_its_after_switch;
        int switch_count= 0;
        // if we fail to solve eqs, we reset status/operability before leaving
        const auto well_status_orig = this->wellStatus_;
        const auto operability_orig = this->operability_status_;
        auto well_status_cur = well_status_orig;
        // don't allow opening wells that are stopped from schedule or has a stopped well state
        const bool allow_open =  this->well_ecl_.getStatus() == WellStatus::OPEN &&
                                 well_state.well(this->index_of_well_).status == WellStatus::OPEN;
        // don't allow switcing for wells under zero rate target or requested fixed status and control
        const bool allow_switching = !this->wellUnderZeroRateTarget(summary_state, well_state) &&
                                     (!fixed_control || !fixed_status) && allow_open;
        bool changed = false;
        bool final_check = false;
        // well needs to be set operable or else solving/updating of re-opened wells is skipped
        this->operability_status_.resetOperability();
        this->operability_status_.solvable = true;

        for (; it < max_iter_number; ++it, ++debug_cost_counter_) {
            its_since_last_switch++;
            if (allow_switching && its_since_last_switch >= min_its_after_switch){
                const double wqTotal = this->primary_variables_.getWQTotal().value();
                changed = this->updateWellControlAndStatusLocalIteration(ebosSimulator, well_state, group_state, inj_controls, prod_controls, wqTotal, deferred_logger, fixed_control, fixed_status);
                if (changed){
                    its_since_last_switch = 0;
                    switch_count++;
                    if (well_status_cur != this->wellStatus_) {
                        well_status_cur = this->wellStatus_;
                    }
                }
                if (!changed && final_check) {
                    break;
                } else {
                    final_check = false;
                }
            }

            assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            const BVectorWell dx_well = this->linSys_.solve();

            if (it > this->param_.strict_inner_iter_wells_) {
                relax_convergence = true;
                this->regularize_ = true;
            }

            const auto report = getWellConvergence(summary_state, well_state, Base::B_avg_, deferred_logger, relax_convergence);
            converged = report.converged();
            if (converged) {
                // if equations are sufficiently linear they might converge in less than min_its_after_switch
                // in this case, make sure all constraints are satisfied before returning
                if (switch_count > 0 && its_since_last_switch < min_its_after_switch) {
                    final_check = true;
                    its_since_last_switch = min_its_after_switch;
                } else {
                    break;
                }
            }

            // getFinteWellResiduals returns false for nan/inf residuals
            {
                const auto& [isFinite, residuals] = this->getFiniteWellResiduals(Base::B_avg_, deferred_logger);
                if (!isFinite)
                    return false;

                residual_history.push_back(residuals);
            }

            if (!converged) {
                measure_history.push_back(this->getResidualMeasureValue(well_state,
                                                                        residual_history[it],
                                                                        this->param_.tolerance_wells_,
                                                                        this->param_.tolerance_pressure_ms_wells_,
                                                                        deferred_logger));

                bool is_oscillate = false;
                bool is_stagnate = false;

                this->detectOscillations(measure_history, is_oscillate, is_stagnate);
                // TODO: maybe we should have more sophisticated strategy to recover the relaxation factor,
                // for example, to recover it to be bigger

                if (is_oscillate || is_stagnate) {
                    // HACK!
                    std::string message;
                    if (relaxation_factor == min_relaxation_factor) {
                        ++stagnate_count;
                        if (false) { // this disables the usage of the relaxed tolerance
                            fmt::format_to(std::back_inserter(message), " Well {} observes severe stagnation and/or oscillation."
                                                                        " We relax the tolerance and check for convergence. \n", this->name());
                            const auto reportStag = getWellConvergence(summary_state, well_state, Base::B_avg_,
                                                                       deferred_logger, true);
                            if (reportStag.converged()) {
                                converged = true;
                                fmt::format_to(std::back_inserter(message), " Well {}  manages to get converged with relaxed tolerances in {} inner iterations", this->name(), it);
                                deferred_logger.debug(message);
                                return converged;
                            }
                        }
                    }

                    // a factor value to reduce the relaxation_factor
                    constexpr double reduction_mutliplier = 0.9;
                    relaxation_factor = std::max(relaxation_factor * reduction_mutliplier, min_relaxation_factor);

                    // debug output
                    if (is_stagnate) {
                        fmt::format_to(std::back_inserter(message), " well {} observes stagnation in inner iteration {}\n", this->name(), it);
                    }
                    if (is_oscillate) {
                        fmt::format_to(std::back_inserter(message), " well {} observes oscillation in inner iteration {}\n", this->name(), it);
                    }
                    fmt::format_to(std::back_inserter(message), " relaxation_factor is {} now\n", relaxation_factor);

                    this->regularize_ = true;
                    deferred_logger.debug(message);
                }
            }
            updateWellState(summary_state, dx_well, well_state, deferred_logger, relaxation_factor);
            initPrimaryVariablesEvaluation();
        }

        if (converged) {
            if (allow_switching){
                // update operability if status change
                const bool is_stopped = this->wellIsStopped();
                if (this->wellHasTHPConstraints(summary_state)){
                    this->operability_status_.can_obtain_bhp_with_thp_limit = !is_stopped;
                    this->operability_status_.obey_thp_limit_under_bhp_limit = !is_stopped;
                } else {
                    this->operability_status_.operable_under_only_bhp_limit = !is_stopped;
                }
            }
            std::string message = fmt::format("   Well {} converged in {} inner iterations ("
                                                    "{} control/status switches).", this->name(), it, switch_count);
            if (relax_convergence) {
                message.append(fmt::format("   (A relaxed tolerance was used after {} iterations)",
                                           this->param_.strict_inner_iter_wells_));
            }
            deferred_logger.debug(message);
        } else {
            this->wellStatus_ = well_status_orig;
            this->operability_status_ = operability_orig;            
            const std::string message = fmt::format("   Well {} did not converge in {} inner iterations ("
                                                    "{} control/status switches).", this->name(), it, switch_count);
            deferred_logger.debug(message);
        }

        return converged;
    }


    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                   const double dt,
                                   const Well::InjectionControls& inj_controls,
                                   const Well::ProductionControls& prod_controls,
                                   WellState& well_state,
                                   const GroupState& group_state,
                                   DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // update the upwinding segments
        this->segments_.updateUpwindingSegments(this->primary_variables_);

        // calculate the fluid properties needed.
        computeSegmentFluidProperties(ebosSimulator, deferred_logger);

        // clear all entries
        this->linSys_.clear();

        auto& ws = well_state.well(this->index_of_well_);
        ws.phase_mixing_rates.fill(0.0);

        // for the black oil cases, there will be four equations,
        // the first three of them are the mass balance equations, the last one is the pressure equations.
        //
        // but for the top segment, the pressure equation will be the well control equation, and the other three will be the same.

        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);

        const int nseg = this->numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the accumulation term
            // TODO: without considering the efficiency factor for now
            {
                const EvalWell segment_surface_volume = getSegmentSurfaceVolume(ebosSimulator, seg);

                // Add a regularization_factor to increase the accumulation term
                // This will make the system less stiff and help convergence for
                // difficult cases
                const Scalar regularization_factor =  this->regularize_? this->param_.regularization_factor_wells_ : 1.0;
                // for each component
                for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                    const EvalWell accumulation_term = regularization_factor * (segment_surface_volume * this->primary_variables_.surfaceVolumeFraction(seg, comp_idx)
                                                     - segment_fluid_initial_[seg][comp_idx]) / dt;
                    MultisegmentWellAssemble(*this).
                        assembleAccumulationTerm(seg, comp_idx, accumulation_term, this->linSys_);
                }
            }
            // considering the contributions due to flowing out from the segment
            {
                const int seg_upwind = this->segments_.upwinding_segment(seg);
                for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                    const EvalWell segment_rate =
                        this->primary_variables_.getSegmentRateUpwinding(seg,
                                                                         seg_upwind,
                                                                         comp_idx) *
                        this->well_efficiency_factor_;
                    MultisegmentWellAssemble(*this).
                        assembleOutflowTerm(seg, seg_upwind, comp_idx, segment_rate, this->linSys_);
                }
            }

            // considering the contributions from the inlet segments
            {
                for (const int inlet : this->segments_.inlets()[seg]) {
                    const int inlet_upwind = this->segments_.upwinding_segment(inlet);
                    for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                        const EvalWell inlet_rate =
                            this->primary_variables_.getSegmentRateUpwinding(inlet,
                                                                             inlet_upwind,
                                                                             comp_idx) *
                            this->well_efficiency_factor_;
                        MultisegmentWellAssemble(*this).
                            assembleInflowTerm(seg, inlet, inlet_upwind, comp_idx, inlet_rate, this->linSys_);
                    }
                }
            }

            // calculating the perforation rate for each perforation that belongs to this segment
            const EvalWell seg_pressure = this->primary_variables_.getSegmentPressure(seg);
            auto& perf_data = ws.perf_data;
            auto& perf_rates = perf_data.phase_rates;
            auto& perf_press_state = perf_data.pressure;
            for (const int perf : this->segments_.perforations()[seg]) {
                const int cell_idx = this->well_cells_[perf];
                const auto& int_quants = ebosSimulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                std::vector<EvalWell> mob(this->num_components_, 0.0);
                getMobility(ebosSimulator, perf, mob, deferred_logger);
                const double trans_mult = ebosSimulator.problem().template wellTransMultiplier<double>(int_quants, cell_idx);
                const auto& wellstate_nupcol = ebosSimulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
                const std::vector<Scalar> Tw = this->wellIndex(perf, int_quants, trans_mult, wellstate_nupcol);
                std::vector<EvalWell> cq_s(this->num_components_, 0.0);
                EvalWell perf_press;
                PerforationRates perfRates;
                computePerfRate(int_quants, mob, Tw, seg, perf, seg_pressure,
                                allow_cf, cq_s, perf_press, perfRates, deferred_logger);

                // updating the solution gas rate and solution oil rate
                if (this->isProducer()) {
                    ws.phase_mixing_rates[ws.dissolved_gas] += perfRates.dis_gas;
                    ws.phase_mixing_rates[ws.vaporized_oil] += perfRates.vap_oil;
                    perf_data.phase_mixing_rates[perf][ws.dissolved_gas] = perfRates.dis_gas;
                    perf_data.phase_mixing_rates[perf][ws.vaporized_oil] = perfRates.vap_oil;
                }

                // store the perf pressure and rates
                for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                    perf_rates[perf*this->number_of_phases_ + this->ebosCompIdxToFlowCompIdx(comp_idx)] = cq_s[comp_idx].value();
                }
                perf_press_state[perf] = perf_press.value();

                for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                    // the cq_s entering mass balance equations need to consider the efficiency factors.
                    const EvalWell cq_s_effective = cq_s[comp_idx] * this->well_efficiency_factor_;

                    this->connectionRates_[perf][comp_idx] = Base::restrictEval(cq_s_effective);

                    MultisegmentWellAssemble(*this).
                        assemblePerforationEq(seg, cell_idx, comp_idx, cq_s_effective, this->linSys_);
                }
            }

            // the fourth dequation, the pressure drop equation
            if (seg == 0) { // top segment, pressure equation is the control equation
                const auto& summaryState = ebosSimulator.vanguard().summaryState();
                const Schedule& schedule = ebosSimulator.vanguard().schedule();
                MultisegmentWellAssemble(*this).
                        assembleControlEq(well_state,
                                        group_state,
                                        schedule,
                                        summaryState,
                                        inj_controls,
                                        prod_controls,
                                        getRefDensity(),
                                        this->primary_variables_,
                                        this->linSys_,
                                        deferred_logger);
            } else {
                const UnitSystem& unit_system = ebosSimulator.vanguard().eclState().getDeckUnitSystem();
                const auto& summary_state = ebosSimulator.vanguard().summaryState();
                this->assemblePressureEq(seg, unit_system, well_state, summary_state, this->param_.use_average_density_ms_wells_, deferred_logger);
            }
        }

        this->linSys_.createSolver();
    }




    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const
    {
        return !this->getAllowCrossFlow() && allDrawDownWrongDirection(ebos_simulator);
    }


    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    allDrawDownWrongDirection(const Simulator& ebos_simulator) const
    {
        bool all_drawdown_wrong_direction = true;
        const int nseg = this->numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            const EvalWell segment_pressure = this->primary_variables_.getSegmentPressure(seg);
            for (const int perf : this->segments_.perforations()[seg]) {

                const int cell_idx = this->well_cells_[perf];
                const auto& intQuants = ebos_simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                const auto& fs = intQuants.fluidState();

                // pressure difference between the segment and the perforation
                const EvalWell perf_seg_press_diff = this->segments_.getPressureDiffSegPerf(seg, perf);
                // pressure difference between the perforation and the grid cell
                const double cell_perf_press_diff = this->cell_perforation_pressure_diffs_[perf];

                const double pressure_cell = this->getPerfCellPressure(fs).value();
                const double perf_press = pressure_cell - cell_perf_press_diff;
                // Pressure drawdown (also used to determine direction of flow)
                // TODO: not 100% sure about the sign of the seg_perf_press_diff
                const EvalWell drawdown = perf_press - (segment_pressure + perf_seg_press_diff);

                // for now, if there is one perforation can produce/inject in the correct
                // direction, we consider this well can still produce/inject.
                // TODO: it can be more complicated than this to cause wrong-signed rates
                if ( (drawdown < 0. && this->isInjector()) ||
                     (drawdown > 0. && this->isProducer()) )  {
                    all_drawdown_wrong_direction = false;
                    break;
                }
            }
        }

        return all_drawdown_wrong_direction;
    }




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWaterThroughput(const double /*dt*/, WellState& /*well_state*/) const
    {
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentSurfaceVolume(const Simulator& ebos_simulator, const int seg_idx) const
    {
        EvalWell temperature;
        EvalWell saltConcentration;
        int pvt_region_index;
        {
            // using the pvt region of first perforated cell
            // TODO: it should be a member of the WellInterface, initialized properly
            const int cell_idx = this->well_cells_[0];
            const auto& intQuants = ebos_simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            temperature.setValue(fs.temperature(FluidSystem::oilPhaseIdx).value());
            saltConcentration = this->extendEval(fs.saltConcentration());
            pvt_region_index = fs.pvtRegionIndex();
        }

        return this->segments_.getSurfaceVolume(temperature,
                                                saltConcentration,
                                                this->primary_variables_,
                                                pvt_region_index,
                                                seg_idx);
    }


    template<typename TypeTag>
    std::optional<double>
    MultisegmentWell<TypeTag>::
    computeBhpAtThpLimitProd(const WellState& well_state,
                             const Simulator& ebos_simulator,
                             const SummaryState& summary_state,
                             DeferredLogger& deferred_logger) const
    {
        return this->MultisegmentWell<TypeTag>::computeBhpAtThpLimitProdWithAlq(
                                               ebos_simulator,
                                               summary_state,
                                               this->getALQ(well_state),
                                               deferred_logger);
    }



    template<typename TypeTag>
    std::optional<double>
    MultisegmentWell<TypeTag>::
    computeBhpAtThpLimitProdWithAlq(const Simulator& ebos_simulator,
                                    const SummaryState& summary_state,
                                    const double alq_value,
                                    DeferredLogger& deferred_logger) const
    {
        // Make the frates() function.
        auto frates = [this, &ebos_simulator, &deferred_logger](const double bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<double> rates(3);
            computeWellRatesWithBhp(ebos_simulator, bhp, rates, deferred_logger);
            return rates;
        };

        auto bhpAtLimit = WellBhpThpCalculator(*this).
               computeBhpAtThpLimitProd(frates,
                                        summary_state,
                                        this->maxPerfPress(ebos_simulator),
                                        this->getRefDensity(),
                                        alq_value,
                                        this->getTHPConstraint(summary_state),
                                        deferred_logger);

       if (bhpAtLimit)
           return bhpAtLimit;

       auto fratesIter = [this, &ebos_simulator, &deferred_logger](const double bhp) {
           // Solver the well iterations to see if we are
           // able to get a solution with an update
           // solution
           std::vector<double> rates(3);
           computeWellRatesWithBhpIterations(ebos_simulator, bhp, rates, deferred_logger);
           return rates;
       };

       return WellBhpThpCalculator(*this).
              computeBhpAtThpLimitProd(fratesIter,
                                       summary_state,
                                       this->maxPerfPress(ebos_simulator),
                                       this->getRefDensity(),
                                       alq_value,
                                       this->getTHPConstraint(summary_state),
                                       deferred_logger);
    }

    template<typename TypeTag>
    std::optional<double>
    MultisegmentWell<TypeTag>::
    computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                            const SummaryState& summary_state,
                            DeferredLogger& deferred_logger) const
    {
        // Make the frates() function.
        auto frates = [this, &ebos_simulator, &deferred_logger](const double bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<double> rates(3);
            computeWellRatesWithBhp(ebos_simulator, bhp, rates, deferred_logger);
            return rates;
        };

        auto bhpAtLimit = WellBhpThpCalculator(*this).
                computeBhpAtThpLimitInj(frates,
                                        summary_state,
                                        this->getRefDensity(),
                                        0.05,
                                        100,
                                        false,
                                        deferred_logger);

        if (bhpAtLimit)
            return bhpAtLimit;

       auto fratesIter = [this, &ebos_simulator, &deferred_logger](const double bhp) {
           // Solver the well iterations to see if we are
           // able to get a solution with an update
           // solution
           std::vector<double> rates(3);
           computeWellRatesWithBhpIterations(ebos_simulator, bhp, rates, deferred_logger);
           return rates;
       };

        return WellBhpThpCalculator(*this).
               computeBhpAtThpLimitInj(fratesIter,
                                       summary_state,
                                       this->getRefDensity(),
                                       0.05,
                                       100,
                                       false,
                                       deferred_logger);
    }





    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    maxPerfPress(const Simulator& ebos_simulator) const
    {
        double max_pressure = 0.0;
        const int nseg = this->numberOfSegments();
        for (int seg = 0; seg < nseg; ++seg) {
            for (const int perf : this->segments_.perforations()[seg]) {
                const int cell_idx = this->well_cells_[perf];
                const auto& int_quants = ebos_simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                const auto& fs = int_quants.fluidState();
                double pressure_cell = this->getPerfCellPressure(fs).value();
                max_pressure = std::max(max_pressure, pressure_cell);
            }
        }
        return max_pressure;
    }





    template<typename TypeTag>
    std::vector<double>
    MultisegmentWell<TypeTag>::
    computeCurrentWellRates(const Simulator& ebosSimulator,
                            DeferredLogger& deferred_logger) const
    {
        // Calculate the rates that follow from the current primary variables.
        std::vector<Scalar> well_q_s(this->num_components_, 0.0);
        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);
        const int nseg = this->numberOfSegments();
        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the perforation rate for each perforation that belongs to this segment
            const Scalar seg_pressure = getValue(this->primary_variables_.getSegmentPressure(seg));
            for (const int perf : this->segments_.perforations()[seg]) {
                const int cell_idx = this->well_cells_[perf];
                const auto& int_quants = ebosSimulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
                std::vector<Scalar> mob(this->num_components_, 0.0);
                getMobility(ebosSimulator, perf, mob, deferred_logger);
                const double trans_mult = ebosSimulator.problem().template wellTransMultiplier<double>(int_quants, cell_idx);
                const auto& wellstate_nupcol = ebosSimulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
                const std::vector<Scalar> Tw = this->wellIndex(perf, int_quants, trans_mult, wellstate_nupcol);
                std::vector<Scalar> cq_s(this->num_components_, 0.0);
                Scalar perf_press = 0.0;
                PerforationRates perf_rates;
                computePerfRate(int_quants, mob, Tw, seg, perf, seg_pressure,
                                allow_cf, cq_s, perf_press, perf_rates, deferred_logger);
                for (int comp = 0; comp < this->num_components_; ++comp) {
                    well_q_s[comp] += cq_s[comp];
                }
            }
        }
        return well_q_s;
    }


    template <typename TypeTag>
    std::vector<double>
    MultisegmentWell<TypeTag>::
    getPrimaryVars() const
    {
        const int num_seg = this->numberOfSegments();
        constexpr int num_eq = MSWEval::numWellEq;
        std::vector<double> retval(num_seg * num_eq);
        for (int ii = 0; ii < num_seg; ++ii) {
            const auto& pv = this->primary_variables_.value(ii);
            std::copy(pv.begin(), pv.end(), retval.begin() + ii * num_eq);
        }
        return retval;
    }




    template <typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    setPrimaryVars(std::vector<double>::const_iterator it)
    {
        const int num_seg = this->numberOfSegments();
        constexpr int num_eq = MSWEval::numWellEq;
        std::array<double, num_eq> tmp;
        for (int ii = 0; ii < num_seg; ++ii) {
            const auto start = it + num_seg * num_eq;
            std::copy(start, start + num_eq, tmp.begin());
            this->primary_variables_.setValue(ii, tmp);
        }
        return num_seg * num_eq;
    }

} // namespace Opm
