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


#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

namespace Opm
{


    template <typename TypeTag>
    MultisegmentWell<TypeTag>::
    MultisegmentWell(const Well2& well, const int time_step, const Wells* wells,
                     const ModelParameters& param,
                     const RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components)
    : Base(well, time_step, wells, param, rate_converter, pvtRegionIdx, num_components)
    , segment_perforations_(numberOfSegments())
    , segment_inlets_(numberOfSegments())
    , cell_perforation_depth_diffs_(number_of_perforations_, 0.0)
    , cell_perforation_pressure_diffs_(number_of_perforations_, 0.0)
    , perforation_segment_depth_diffs_(number_of_perforations_, 0.0)
    , segment_fluid_initial_(numberOfSegments(), std::vector<double>(num_components_, 0.0))
    , segment_densities_(numberOfSegments(), 0.0)
    , segment_viscosities_(numberOfSegments(), 0.0)
    , segment_mass_rates_(numberOfSegments(), 0.0)
    , segment_depth_diffs_(numberOfSegments(), 0.0)
    , upwinding_segments_(numberOfSegments(), 0)
    {
        // not handling solvent or polymer for now with multisegment well
        if (has_solvent) {
            OPM_THROW(std::runtime_error, "solvent is not supported by multisegment well yet");
        }

        if (has_polymer) {
            OPM_THROW(std::runtime_error, "polymer is not supported by multisegment well yet");
        }

        if (Base::has_energy) {
            OPM_THROW(std::runtime_error, "energy is not supported by multisegment well yet");
        }
        // since we decide to use the WellSegments from the well parser. we can reuse a lot from it.
        // for other facilities needed but not available from parser, we need to process them here

        // initialize the segment_perforations_ and update perforation_segment_depth_diffs_
        const WellConnections& completion_set = well_ecl_.getConnections();
        // index of the perforation within wells struct
        // there might be some perforations not active, which causes the number of the perforations in
        // well_ecl_ and wells struct different
        // the current implementation is a temporary solution for now, it should be corrected from the parser
        // side
        int i_perf_wells = 0;
        perf_depth_.resize(number_of_perforations_, 0.);
        for (size_t perf = 0; perf < completion_set.size(); ++perf) {
            const Connection& connection = completion_set.get(perf);
            if (connection.state() == WellCompletion::OPEN) {
                const int segment_index = segmentNumberToIndex(connection.segment());
                segment_perforations_[segment_index].push_back(i_perf_wells);
                perf_depth_[i_perf_wells] = connection.depth();
                const double segment_depth = segmentSet()[segment_index].depth();
                perforation_segment_depth_diffs_[i_perf_wells] = perf_depth_[i_perf_wells] - segment_depth;
                i_perf_wells++;
            }
        }

        // initialize the segment_inlets_
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            const Segment& segment = segmentSet()[seg];
            const int segment_number = segment.segmentNumber();
            const int outlet_segment_number = segment.outletSegment();
            if (outlet_segment_number > 0) {
                const int segment_index = segmentNumberToIndex(segment_number);
                const int outlet_segment_index = segmentNumberToIndex(outlet_segment_number);
                segment_inlets_[outlet_segment_index].push_back(segment_index);
            }
        }

        // calculating the depth difference between the segment and its oulet_segments
        // for the top segment, we will make its zero unless we find other purpose to use this value
        for (int seg = 1; seg < numberOfSegments(); ++seg) {
            const double segment_depth = segmentSet()[seg].depth();
            const int outlet_segment_number = segmentSet()[seg].outletSegment();
            const Segment& outlet_segment = segmentSet()[segmentNumberToIndex(outlet_segment_number)];
            const double outlet_depth = outlet_segment.depth();
            segment_depth_diffs_[seg] = segment_depth - outlet_depth;
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, num_cells);

        // TODO: for StandardWell, we need to update the perf depth here using depth_arg.
        // for MultisegmentWell, it is much more complicated.
        // It can be specified directly, it can be calculated from the segment depth,
        // it can also use the cell center, which is the same for StandardWell.
        // For the last case, should we update the depth with the depth_arg? For the
        // future, it can be a source of wrong result with Multisegment well.
        // An indicator from the opm-parser should indicate what kind of depth we should use here.

        // \Note: we do not update the depth here. And it looks like for now, we only have the option to use
        // specified perforation depth
        initMatrixAndVectors(num_cells);

        // calcuate the depth difference between the perforations and the perforated grid block
        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            cell_perforation_depth_diffs_[perf] = depth_arg[cell_idx] - perf_depth_[perf];
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initMatrixAndVectors(const int num_cells) const
    {
        duneB_.setBuildMode( OffDiagMatWell::row_wise );
        duneC_.setBuildMode( OffDiagMatWell::row_wise );
        duneD_.setBuildMode( DiagMatWell::row_wise );

        // set the size and patterns for all the matrices and vectors
        // [A C^T    [x    =  [ res
        //  B D] x_well]      res_well]

        // calculatiing the NNZ for duneD_
        // NNZ = number_of_segments + 2 * (number_of_inlets / number_of_outlets)
        {
            int nnz_d = numberOfSegments();
            for (const std::vector<int>& inlets : segment_inlets_) {
                nnz_d += 2 * inlets.size();
            }
            duneD_.setSize(numberOfSegments(), numberOfSegments(), nnz_d);
        }
        duneB_.setSize(numberOfSegments(), num_cells, number_of_perforations_);
        duneC_.setSize(numberOfSegments(), num_cells, number_of_perforations_);

        // we need to add the off diagonal ones
        for (auto row = duneD_.createbegin(), end = duneD_.createend(); row != end; ++row) {
            // the number of the row corrspnds to the segment now
            const int seg = row.index();
            // adding the item related to outlet relation
            const Segment& segment = segmentSet()[seg];
            const int outlet_segment_number = segment.outletSegment();
            if (outlet_segment_number > 0) { // if there is a outlet_segment
                const int outlet_segment_index = segmentNumberToIndex(outlet_segment_number);
                row.insert(outlet_segment_index);
            }

            // Add nonzeros for diagonal
            row.insert(seg);

            // insert the item related to its inlets
            for (const int& inlet : segment_inlets_[seg]) {
                row.insert(inlet);
            }
        }

        // make the C matrix
        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row != end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (const int& perf : segment_perforations_[row.index()]) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        // make the B^T matrix
        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row != end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (const int& perf : segment_perforations_[row.index()]) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize( numberOfSegments() );

        primary_variables_.resize(numberOfSegments());
        primary_variables_evaluation_.resize(numberOfSegments());
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
                primary_variables_evaluation_[seg][eq_idx] = 0.0;
                primary_variables_evaluation_[seg][eq_idx].setValue(primary_variables_[seg][eq_idx]);
                primary_variables_evaluation_[seg][eq_idx].setDerivative(eq_idx + numEq, 1.0);
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleWellEq(const Simulator& ebosSimulator,
                   const std::vector<Scalar>& B_avg,
                   const double dt,
                   WellState& well_state,
                   Opm::DeferredLogger& deferred_logger)
    {

        const bool use_inner_iterations = param_.use_inner_iterations_ms_wells_;
        if (use_inner_iterations) {

            iterateWellEquations(ebosSimulator, B_avg, dt, well_state, deferred_logger);
        }

        assembleWellEqWithoutIteration(ebosSimulator, dt, well_state, deferred_logger);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateWithTarget(const Simulator& /* ebos_simulator */,
                              WellState& well_state,
                              Opm::DeferredLogger& /* deferred_logger */) const
    {
        // Updating well state bas on well control
        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
        const int current = well_state.currentControls()[index_of_well_];
        const double target = well_controls_iget_target(well_controls_, current);
        const double* distr = well_controls_iget_distr(well_controls_, current);
        switch (well_controls_iget_type(well_controls_, current)) {
        case BHP: {
            well_state.bhp()[index_of_well_] = target;
            const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
            well_state.segPress()[top_segment_index] = well_state.bhp()[index_of_well_];
            // TODO: similar to the way below to handle THP
            // we should not something related to thp here when there is thp constraint
            break;
        }

        case THP: {
            well_state.thp()[index_of_well_] = target;

            /* const Opm::PhaseUsage& pu = phaseUsage();
            std::vector<double> rates(3, 0.0);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                rates[ Water ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Water ] ];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                 rates[ Oil ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Oil ] ];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                rates[ Gas ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Gas ] ];
            } */

            // const int table_id = well_controls_iget_vfp(well_controls_, current);
            // const double& thp    = well_controls_iget_target(well_controls_, current);
            // const double& alq    = well_controls_iget_alq(well_controls_, current);

            // TODO: implement calculateBhpFromThp function
            // well_state.bhp()[index_of_well_] = calculateBhpFromThp(rates, current);
            // also the top segment pressure
            /* const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
            well_state.segPress()[top_segment_index] = well_state.bhp()[index_of_well_]; */
            break;
        }

        case RESERVOIR_RATE: // intentional fall-through
        case SURFACE_RATE:
            // for the update of the rates, after we update the well rates, we can try to scale
            // the segment rates and perforation rates with the same factor
            // or the other way, we can use the same approach like the initialization of the well state,
            // we devide the well rates to update the perforation rates, then we update the segment rates based
            // on the perforation rates.
            // the second way is safer, since if the well control is changed, the first way will not guarantee the consistence
            // of between the segment rates and peforation rates

            // checking the number of the phases under control
            int numPhasesWithTargetsUnderThisControl = 0;
            for (int phase = 0; phase < number_of_phases_; ++phase) {
                if (distr[phase] > 0.0) {
                    numPhasesWithTargetsUnderThisControl += 1;
                }
            }

            assert(numPhasesWithTargetsUnderThisControl > 0);

            if (well_type_ == INJECTOR) {
                // assign target value as initial guess for injectors
                // only handles single phase control at the moment
                assert(numPhasesWithTargetsUnderThisControl == 1);

                for (int phase = 0; phase < number_of_phases_; ++phase) {
                    if (distr[phase] > 0.) {
                        well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] = target / distr[phase];
                    } else {
                        well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] = 0.;
                    }
                }

                initSegmentRatesWithWellRates(well_state);

            } else if (well_type_ == PRODUCER) {
                // update the rates of phases under control based on the target,
                // and also update rates of phases not under control to keep the rate ratio,
                // assuming the mobility ratio does not change for the production wells
                double original_rates_under_phase_control = 0.0;
                for (int phase = 0; phase < number_of_phases_; ++phase) {
                    if (distr[phase] > 0.0) {
                        original_rates_under_phase_control += well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] * distr[phase];
                    }
                }

                if (original_rates_under_phase_control != 0.0 ) {
                    double scaling_factor = target / original_rates_under_phase_control;

                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] *= scaling_factor;

                        // scaling the segment rates with the same way with well rates
                        const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
                        for (int seg = 0; seg < numberOfSegments(); ++seg) {
                             well_state.segRates()[number_of_phases_ * (seg + top_segment_index) + phase] *= scaling_factor;
                        }
                    }
                } else { // scaling factor is not well defined when original_rates_under_phase_control is zero
                    // separating targets equally between phases under control
                    const double target_rate_divided = target / numPhasesWithTargetsUnderThisControl;
                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        if (distr[phase] > 0.0) {
                            well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] = target_rate_divided / distr[phase];
                        } else {
                            // this only happens for SURFACE_RATE control
                            well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] = target_rate_divided;
                        }
                    }
                    initSegmentRatesWithWellRates(well_state);
                }
            }

            break;
        } // end of switch
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initSegmentRatesWithWellRates(WellState& well_state) const
    {
        for (int phase = 0; phase < number_of_phases_; ++phase) {
            const double perf_phaserate = well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] / number_of_perforations_;
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                well_state.perfPhaseRates()[number_of_phases_ * (first_perf_ + perf) + phase] = perf_phaserate;
            }
        }

        const std::vector<double> perforation_rates(well_state.perfPhaseRates().begin() + number_of_phases_ * first_perf_,
                                                    well_state.perfPhaseRates().begin() +
                                                    number_of_phases_ * (first_perf_ + number_of_perforations_) );
        std::vector<double> segment_rates;
        WellState::calculateSegmentRates(segment_inlets_, segment_perforations_, perforation_rates, number_of_phases_,
                                         0, segment_rates);
        const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
        std::copy(segment_rates.begin(), segment_rates.end(),
                  well_state.segRates().begin() + number_of_phases_ * top_segment_index );
        // we need to check the top segment rates should be same with the well rates
    }





    template <typename TypeTag>
    ConvergenceReport
    MultisegmentWell<TypeTag>::
    getWellConvergence(const std::vector<double>& B_avg, Opm::DeferredLogger& deferred_logger) const
    {
        assert(int(B_avg.size()) == num_components_);

        // checking if any residual is NaN or too large. The two large one is only handled for the well flux
        std::vector<std::vector<double>> abs_residual(numberOfSegments(), std::vector<double>(numWellEq, 0.0));
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
                abs_residual[seg][eq_idx] = std::abs(resWell_[seg][eq_idx]);
            }
        }

        std::vector<double> maximum_residual(numWellEq, 0.0);

        ConvergenceReport report;
        // TODO: the following is a little complicated, maybe can be simplified in some way?
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                if (eq_idx < num_components_) { // phase or component mass equations
                    const double flux_residual = B_avg[eq_idx] * abs_residual[seg][eq_idx];
                    if (flux_residual > maximum_residual[eq_idx]) {
                        maximum_residual[eq_idx] = flux_residual;
                    }
                } else { // pressure or control equation
                    // for the top segment (seg == 0), it is control equation, will be checked later separately
                    if (seg > 0) {
                        // Pressure equation
                        const double pressure_residual = abs_residual[seg][eq_idx];
                        if (pressure_residual > maximum_residual[eq_idx]) {
                            maximum_residual[eq_idx] = pressure_residual;
                        }
                    }
                }
            }
        }

        using CR = ConvergenceReport;
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            if (eq_idx < num_components_) { // phase or component mass equations
                const double flux_residual = maximum_residual[eq_idx];
                // TODO: the report can not handle the segment number yet.
                if (std::isnan(flux_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::NotANumber, eq_idx, name()});
                } else if (flux_residual > param_.max_residual_allowed_) {
                    report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::TooLarge, eq_idx, name()});
                } else if (flux_residual > param_.tolerance_wells_) {
                    report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, eq_idx, name()});
                }
            } else { // pressure equation
                const double pressure_residual = maximum_residual[eq_idx];
                const int dummy_component = -1;
                if (std::isnan(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::NotANumber, dummy_component, name()});
                } else if (std::isinf(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::TooLarge, dummy_component, name()});
                } else if (pressure_residual > param_.tolerance_pressure_ms_wells_) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, name()});
                }
            }
        }

        checkConvergenceControlEq(report, deferred_logger);

        return report;
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        BVectorWell Bx(duneB_.N());

        duneB_.mv(x, Bx);

        // invDBx = duneD^-1 * Bx_
        const BVectorWell invDBx = mswellhelpers::invDXDirect(duneD_, Bx);

        // Ax = Ax - duneC_^T * invDBx
        duneC_.mmtv(invDBx,Ax);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(BVector& r) const
    {
        // invDrw_ = duneD^-1 * resWell_
        const BVectorWell invDrw = mswellhelpers::invDXDirect(duneD_, resWell_);
        // r = r - duneC_^T * invDrw
        duneC_.mmtv(invDrw, r);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x,
                                          WellState& well_state,
                                          Opm::DeferredLogger& deferred_logger) const
    {
        BVectorWell xw(1);
        recoverSolutionWell(x, xw);
        updateWellState(xw, well_state, deferred_logger);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const std::vector<Scalar>& B_avg,
                          const WellState& well_state,
                          std::vector<double>& well_potentials,
                          Opm::DeferredLogger& deferred_logger)
    {
        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

        updatePrimaryVariables(well_state, deferred_logger);

        // initialize the primary variables in Evaluation, which is used in computePerfRate for computeWellPotentials
        // TODO: for computeWellPotentials, no derivative is required actually
        initPrimaryVariablesEvaluation();

        // get the bhp value based on the bhp constraints
        const double bhp = Base::mostStrictBhpFromBhpLimits(deferred_logger);

        // does the well have a THP related constraint?
        if ( !Base::wellHasTHPConstraints() ) {
            assert(std::abs(bhp) != std::numeric_limits<double>::max());

            computeWellRatesWithBhpPotential(ebosSimulator, B_avg, bhp, well_potentials, deferred_logger);
        } else {

            const std::string msg = std::string("Well potential calculation is not supported for thp controlled multisegment wells \n")
                    + "A well potential of zero is returned for output purposes. \n"
                    + "If you need well potential computed from thp to set the guide rate for group controled wells \n"
                    + "you will have to change the " + name() + " well to a standard well \n";

            deferred_logger.warning("WELL_POTENTIAL_FOR_THP_NOT_IMPLEMENTED_FOR_MULTISEG_WELLS", msg);
            return;
        }

    }


    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesWithBhpPotential(const Simulator& ebosSimulator,
                                     const std::vector<Scalar>& B_avg,
                                     const double& bhp,
                                     std::vector<double>& well_flux,
                                     Opm::DeferredLogger& deferred_logger)
    {

        WellControls* wc = well_controls_;
        const int bhp_index = Base::getControlIndex(BHP);
        const double orig_bhp = well_controls_iget_target(wc, bhp_index);
        const auto orig_current = well_controls_get_current(wc);

        well_controls_iset_target(wc, bhp_index, bhp);
        well_controls_set_current(wc, bhp_index);

        // store a copy of the well state, we don't want to update the real well state
        WellState copy = ebosSimulator.problem().wellModel().wellState();

        initPrimaryVariablesEvaluation();

        const double dt = ebosSimulator.timeStepSize();
        // iterate to get a solution that satisfies the bhp potential.
        iterateWellEquations(ebosSimulator, B_avg, dt, copy, deferred_logger);

        // compute the potential and store in the flux vector.
        const int np = number_of_phases_;
        well_flux.resize(np, 0.0);
        for(int compIdx = 0; compIdx < num_components_; ++compIdx) {
            const EvalWell rate = getSegmentRate(0, compIdx);
            well_flux[ebosCompIdxToFlowCompIdx(compIdx)] += rate.value();
        }

        // reset bhp limit
        well_controls_iset_target(wc, bhp_index, orig_bhp);
        well_controls_set_current(wc, orig_current);

    }



    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& /* deferred_logger */) const
    {
        // TODO: to test using rate conversion coefficients to see if it will be better than
        // this default one

        // the index of the top segment in the WellState
        const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
        const std::vector<double>& segment_rates = well_state.segRates();
        const PhaseUsage& pu = phaseUsage();

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // calculate the total rate for each segment
            double total_seg_rate = 0.0;
            const int seg_index = top_segment_index + seg;
            // the segment pressure
            primary_variables_[seg][SPres] = well_state.segPress()[seg_index];
            // TODO: under what kind of circustances, the following will be wrong?
            // the definition of g makes the gas phase is always the last phase
            for (int p = 0; p < number_of_phases_; p++) {
                total_seg_rate += scalingFactor(p) * segment_rates[number_of_phases_ * seg_index + p];
            }

            primary_variables_[seg][GTotal] = total_seg_rate;
            if (std::abs(total_seg_rate) > 0.) {
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    const int water_pos = pu.phase_pos[Water];
                    primary_variables_[seg][WFrac] = scalingFactor(water_pos) * segment_rates[number_of_phases_ * seg_index + water_pos] / total_seg_rate;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const int gas_pos = pu.phase_pos[Gas];
                    primary_variables_[seg][GFrac] = scalingFactor(gas_pos) * segment_rates[number_of_phases_ * seg_index + gas_pos] / total_seg_rate;
                }
            } else { // total_seg_rate == 0
                if (well_type_ == INJECTOR) {
                    // only single phase injection handled
                    const double* distr = well_controls_get_current_distr(well_controls_);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        if (distr[pu.phase_pos[Water]] > 0.0) {
                            primary_variables_[seg][WFrac] = 1.0;
                        } else {
                            primary_variables_[seg][WFrac] = 0.0;
                        }
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        if (distr[pu.phase_pos[Gas]] > 0.0) {
                            primary_variables_[seg][GFrac] = 1.0;
                        } else {
                            primary_variables_[seg][GFrac] = 0.0;
                        }
                    }
                } else if (well_type_ == PRODUCER) { // producers
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        primary_variables_[seg][WFrac] = 1.0 / number_of_phases_;
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        primary_variables_[seg][GFrac] = 1.0 / number_of_phases_;
                    }
                }
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    recoverSolutionWell(const BVector& x, BVectorWell& xw) const
    {
        BVectorWell resWell = resWell_;
        // resWell = resWell - B * x
        duneB_.mmv(x, resWell);
        // xw = D^-1 * resWell
        xw = mswellhelpers::invDXDirect(duneD_, resWell);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    solveEqAndUpdateWellState(WellState& well_state, Opm::DeferredLogger& deferred_logger)
    {
        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        const BVectorWell dx_well = mswellhelpers::invDXDirect(duneD_, resWell_);

        updateWellState(dx_well, well_state, deferred_logger);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computePerfCellPressDiffs(const Simulator& ebosSimulator)
    {
        for (int perf = 0; perf < number_of_perforations_; ++perf) {

            std::vector<double> kr(number_of_phases_, 0.0);
            std::vector<double> density(number_of_phases_, 0.0);

            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = intQuants.fluidState();

            double sum_kr = 0.;

            const PhaseUsage& pu = phaseUsage();
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
            for (int p = 0; p < number_of_phases_; ++p) {
                average_density += kr[p] * density[p];
            }
            average_density /= sum_kr;

            cell_perforation_pressure_diffs_[perf] = gravity_ * average_density * cell_perforation_depth_diffs_[perf];
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeInitialSegmentFluids(const Simulator& ebos_simulator)
    {
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // TODO: trying to reduce the times for the surfaceVolumeFraction calculation
            const double surface_volume = getSegmentSurfaceVolume(ebos_simulator, seg).value();
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                segment_fluid_initial_[seg][comp_idx] = surface_volume * surfaceVolumeFraction(seg, comp_idx).value();
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellState(const BVectorWell& dwells,
                    WellState& well_state,
                    Opm::DeferredLogger& deferred_logger,
                    const double relaxation_factor) const
    {
        const double dFLimit = param_.dwell_fraction_max_;
        const double max_pressure_change = param_.max_pressure_change_ms_wells_;
        const std::vector<std::array<double, numWellEq> > old_primary_variables = primary_variables_;

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
                // const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]), relaxation_factor * dFLimit);
                const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]) * relaxation_factor, dFLimit);
                primary_variables_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
                // const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]), relaxation_factor * dFLimit);
                const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]) * relaxation_factor, dFLimit);
                primary_variables_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
            }

            // handling the overshooting or undershooting of the fractions
            processFractions(seg);

            // update the segment pressure
            {
                const int sign = dwells[seg][SPres] > 0.? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]), relaxation_factor * max_pressure_change);
                // const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change);
                primary_variables_[seg][SPres] = std::max( old_primary_variables[seg][SPres] - dx_limited, 1e5);
            }

            // update the total rate // TODO: should we have a limitation of the total rate change?
            {
                primary_variables_[seg][GTotal] = old_primary_variables[seg][GTotal] - relaxation_factor * dwells[seg][GTotal];
            }

        }

        updateWellStateFromPrimaryVariables(well_state);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state,
                                Opm::DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(well_state, deferred_logger);
        initPrimaryVariablesEvaluation();
        computePerfCellPressDiffs(ebosSimulator);
        computeInitialSegmentFluids(ebosSimulator);
    }





    template <typename TypeTag>
    const WellSegments&
    MultisegmentWell<TypeTag>::
    segmentSet() const
    {
        return well_ecl_.getSegments();
    }





    template <typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    numberOfSegments() const
    {
        return segmentSet().size();
    }





    template <typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    numberOfPerforations() const
    {
        return segmentSet().number_of_perforations_;
    }





    template <typename TypeTag>
    WellSegment::CompPressureDropEnum
    MultisegmentWell<TypeTag>::
    compPressureDrop() const
    {
        return segmentSet().compPressureDrop();
    }





    template <typename TypeTag>
    WellSegment::MultiPhaseModelEnum
    MultisegmentWell<TypeTag>::
    multiphaseModel() const
    {
        return segmentSet().multiPhaseModel();
    }





    template <typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    segmentNumberToIndex(const int segment_number) const
    {
        return segmentSet().segmentNumberToIndex(segment_number);
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    volumeFraction(const int seg, const unsigned compIdx) const
    {

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
            return primary_variables_evaluation_[seg][WFrac];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
            return primary_variables_evaluation_[seg][GFrac];
        }

        // Oil fraction
        EvalWell oil_fraction = 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            oil_fraction -= primary_variables_evaluation_[seg][WFrac];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            oil_fraction -= primary_variables_evaluation_[seg][GFrac];
        }
        /* if (has_solvent) {
            oil_fraction -= primary_variables_evaluation_[seg][SFrac];
        } */
        return oil_fraction;
    }




    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    volumeFractionScaled(const int seg, const int comp_idx) const
    {
        // For reservoir rate control, the distr in well control is used for the
        // rate conversion coefficients. For the injection well, only the distr of the injection
        // phase is not zero.
        const double scale = scalingFactor(ebosCompIdxToFlowCompIdx(comp_idx));
        if (scale > 0.) {
            return volumeFraction(seg, comp_idx) / scale;
        }

        return volumeFraction(seg, comp_idx);
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    surfaceVolumeFraction(const int seg, const int comp_idx) const
    {
        EvalWell sum_volume_fraction_scaled = 0.;
        for (int idx = 0; idx < num_components_; ++idx) {
            sum_volume_fraction_scaled += volumeFractionScaled(seg, idx);
        }

        assert(sum_volume_fraction_scaled.value() != 0.);

        return volumeFractionScaled(seg, comp_idx) / sum_volume_fraction_scaled;
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computePerfRatePressure(const IntensiveQuantities& int_quants,
                            const std::vector<EvalWell>& mob_perfcells,
                            const int seg,
                            const int perf,
                            const EvalWell& segment_pressure,
                            const bool& allow_cf,
                            std::vector<EvalWell>& cq_s,
                            EvalWell& perf_press,
                            double& perf_dis_gas_rate,
                            double& perf_vap_oil_rate,
                            Opm::DeferredLogger& deferred_logger) const

    {
        std::vector<EvalWell> cmix_s(num_components_, 0.0);

        // the composition of the components inside wellbore
        for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
            cmix_s[comp_idx] = surfaceVolumeFraction(seg, comp_idx);
        }

        const auto& fs = int_quants.fluidState();

        const EvalWell pressure_cell = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        const EvalWell rs = extendEval(fs.Rs());
        const EvalWell rv = extendEval(fs.Rv());

        // not using number_of_phases_ because of solvent
        std::vector<EvalWell> b_perfcells(num_components_, 0.0);

        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells[compIdx] = extendEval(fs.invB(phaseIdx));
        }

        // pressure difference between the segment and the perforation
        const EvalWell perf_seg_press_diff = gravity_ * segment_densities_[seg] * perforation_segment_depth_diffs_[perf];
        // pressure difference between the perforation and the grid cell
        const double cell_perf_press_diff = cell_perforation_pressure_diffs_[perf];

        perf_press = pressure_cell - cell_perf_press_diff;
        // Pressure drawdown (also used to determine direction of flow)
        // TODO: not 100% sure about the sign of the seg_perf_press_diff
        const EvalWell drawdown = perf_press - (segment_pressure + perf_seg_press_diff);

        // producing perforations
        if ( drawdown > 0.0) {
            // Do nothing is crossflow is not allowed
            if (!allow_cf && well_type_ == INJECTOR) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell cq_p = - well_index_[perf] * (mob_perfcells[comp_idx] * drawdown);
                cq_s[comp_idx] = b_perfcells[comp_idx] * cq_p;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const EvalWell cq_s_oil = cq_s[oilCompIdx];
                const EvalWell cq_s_gas = cq_s[gasCompIdx];
                cq_s[gasCompIdx] += rs * cq_s_oil;
                cq_s[oilCompIdx] += rv * cq_s_gas;
            }
        } else { // injecting perforations
            // Do nothing if crossflow is not allowed
            if (!allow_cf && well_type_ == PRODUCER) {
                return;
            }

            // for injecting perforations, we use total mobility
            EvalWell total_mob = mob_perfcells[0];
            for (int comp_idx = 1; comp_idx < num_components_; ++comp_idx) {
                total_mob += mob_perfcells[comp_idx];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - well_index_[perf] * (total_mob * drawdown);

            // compute volume ratio between connection and at standard conditions
            EvalWell volume_ratio = 0.0;
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
                const EvalWell d = 1.0 - rv * rs;

                if (d.value() == 0.0) {
                    OPM_DEFLOG_THROW(Opm::NumericalIssue, "Zero d value obtained for well " << name() << " during flux calcuation"
                                                  << " with rs " << rs << " and rv " << rv, deferred_logger);
                }

                const EvalWell tmp_oil = (cmix_s[oilCompIdx] - rv * cmix_s[gasCompIdx]) / d;
                volume_ratio += tmp_oil / b_perfcells[oilCompIdx];

                const EvalWell tmp_gas = (cmix_s[gasCompIdx] - rs * cmix_s[oilCompIdx]) / d;
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
            EvalWell cqt_is = cqt_i / volume_ratio;
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                cq_s[comp_idx] = cmix_s[comp_idx] * cqt_is;
            }
        } // end for injection perforations

        // calculating the perforation solution gas rate and solution oil rates
        if (well_type_ == PRODUCER) {
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

                const double d = 1.0 - rv.value() * rs.value();
                // vaporized oil into gas
                // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
                perf_vap_oil_rate = rv.value() * (cq_s[gasCompIdx].value() - rs.value() * cq_s[oilCompIdx].value()) / d;
                // dissolved of gas in oil
                // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
                perf_dis_gas_rate = rs.value() * (cq_s[oilCompIdx].value() - rv.value() * cq_s[gasCompIdx].value()) / d;
            }
        }
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    extendEval(const Eval& in) const
    {
        EvalWell out = 0.0;
        out.setValue(in.value());
        for(int eq_idx = 0; eq_idx < numEq;++eq_idx) {
            out.setDerivative(eq_idx, in.derivative(eq_idx));
        }
        return out;
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeSegmentFluidProperties(const Simulator& ebosSimulator)
    {
        // TODO: the concept of phases and components are rather confusing in this function.
        // needs to be addressed sooner or later.

        // get the temperature for later use. It is only useful when we are not handling
        // thermal related simulation
        // basically, it is a single value for all the segments

        EvalWell temperature;
        // not sure how to handle the pvt region related to segment
        // for the current approach, we use the pvt region of the first perforated cell
        // although there are some text indicating using the pvt region of the lowest
        // perforated cell
        // TODO: later to investigate how to handle the pvt region
        int pvt_region_index;
        {
            // using the first perforated cell
            const int cell_idx = well_cells_[0];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            temperature.setValue(fs.temperature(FluidSystem::oilPhaseIdx).value());
            pvt_region_index = fs.pvtRegionIndex();
        }

        std::vector<double> surf_dens(num_components_);
        // Surface density.
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            surf_dens[compIdx] = FluidSystem::referenceDensity( phaseIdx, pvt_region_index );
        }

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // the compostion of the components inside wellbore under surface condition
            std::vector<EvalWell> mix_s(num_components_, 0.0);
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                mix_s[comp_idx] = surfaceVolumeFraction(seg, comp_idx);
            }

            std::vector<EvalWell> b(num_components_, 0.0);
            std::vector<EvalWell> visc(num_components_, 0.0);

            const EvalWell seg_pressure = getSegmentPressure(seg);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                b[waterCompIdx] =
                    FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                visc[waterCompIdx] =
                    FluidSystem::waterPvt().viscosity(pvt_region_index, temperature, seg_pressure);
            }

            EvalWell rv(0.0);
            // gas phase
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    const EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index, temperature, seg_pressure);
                    if (mix_s[oilCompIdx] > 0.0) {
                        if (mix_s[gasCompIdx] > 0.0) {
                            rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                        }

                        if (rv > rvmax) {
                            rv = rvmax;
                        }
                        b[gasCompIdx] =
                            FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rv);
                        visc[gasCompIdx] =
                            FluidSystem::gasPvt().viscosity(pvt_region_index, temperature, seg_pressure, rv);
                    } else { // no oil exists
                        b[gasCompIdx] =
                            FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                        visc[gasCompIdx] =
                            FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    }
                } else { // no Liquid phase
                    // it is the same with zero mix_s[Oil]
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                }
            }

            EvalWell rs(0.0);
            // oil phase
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    const EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index, temperature, seg_pressure);
                    if (mix_s[gasCompIdx] > 0.0) {
                        if (mix_s[oilCompIdx] > 0.0) {
                            rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                        }

                        if (rs > rsmax) {
                            rs = rsmax;
                        }
                        b[oilCompIdx] =
                            FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rs);
                        visc[oilCompIdx] =
                            FluidSystem::oilPvt().viscosity(pvt_region_index, temperature, seg_pressure, rs);
                    } else { // no oil exists
                        b[oilCompIdx] =
                            FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                        visc[oilCompIdx] =
                            FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    }
                } else { // no Liquid phase
                    // it is the same with zero mix_s[Oil]
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                }
            }

            std::vector<EvalWell> mix(mix_s);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

                if (rs != 0.0) { // rs > 0.0?
                    mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / (1. - rs * rv);
                }
                if (rv != 0.0) { // rv > 0.0?
                    mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / (1. - rs * rv);
                }
            }

            EvalWell volrat(0.0);
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                volrat += mix[comp_idx] / b[comp_idx];
            }

            segment_viscosities_[seg] = 0.;
            // calculate the average viscosity
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell comp_fraction = mix[comp_idx] / b[comp_idx] / volrat;
                segment_viscosities_[seg] += visc[comp_idx] * comp_fraction;
            }

            EvalWell density(0.0);
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                density += surf_dens[comp_idx] * mix_s[comp_idx];
            }
            segment_densities_[seg] = density / volrat;

            // calculate the mass rates
            // TODO: for now, we are not considering the upwinding for this amount
            // since how to address the fact that the derivatives is not trivial for now
            // and segment_mass_rates_ goes a long way with the frictional pressure loss
            // and accelerational pressure loss, which needs some work to handle
            segment_mass_rates_[seg] = 0.;
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell rate = getSegmentRate(seg, comp_idx);
                segment_mass_rates_[seg] += rate * surf_dens[comp_idx];
            }
        }
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentPressure(const int seg) const
    {
        return primary_variables_evaluation_[seg][SPres];
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentRate(const int seg,
                   const int comp_idx) const
    {
        return primary_variables_evaluation_[seg][GTotal] * volumeFractionScaled(seg, comp_idx);
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentRateUpwinding(const int seg,
                            const int comp_idx) const
    {
        const int seg_upwind = upwinding_segments_[seg];
        // the result will contain the derivative with resepct to GTotal in segment seg,
        // and the derivatives with respect to WFrac GFrac in segment seg_upwind.
        // the derivative with respect to SPres should be zero.
        const EvalWell segment_rate = primary_variables_evaluation_[seg][GTotal] * volumeFractionScaled(seg_upwind, comp_idx);

        assert(segment_rate.derivative(SPres + numEq) == 0.);

        return segment_rate;
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentGTotal(const int seg) const
    {
        return primary_variables_evaluation_[seg][GTotal];
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    getMobility(const Simulator& ebosSimulator,
                const int perf,
                std::vector<EvalWell>& mob) const
    {
        // TODO: most of this function, if not the whole function, can be moved to the base class
        const int cell_idx = well_cells_[perf];
        assert (int(mob.size()) == num_components_);
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = extendEval(intQuants.mobility(phaseIdx));
            }
            // if (has_solvent) {
            //     mob[contiSolventEqIdx] = extendEval(intQuants.solventMobility());
            // }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            Eval relativePerms[3] = { 0.0, 0.0, 0.0 };
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

        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleControlEq(Opm::DeferredLogger& deferred_logger) const
    {
        EvalWell control_eq(0.0);

        switch (well_controls_get_current_type(well_controls_)) {
            case THP: // not handling this one for now
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Not handling THP control for Multisegment wells for now", deferred_logger);
            }
            case BHP:
            {
                const double target_bhp = well_controls_get_current_target(well_controls_);
                control_eq = getSegmentPressure(0) - target_bhp;
                break;
            }
            case SURFACE_RATE:
            {
                // finding if it is a single phase control or combined phase control
                int number_phases_under_control = 0;
                const double* distr = well_controls_get_current_distr(well_controls_);
                for (int phase = 0; phase < number_of_phases_; ++phase) {
                    if (distr[phase] > 0.0) {
                        ++number_phases_under_control;
                    }
                }
                assert(number_phases_under_control > 0);

                const std::vector<double> g = {1.0, 1.0, 0.01};
                const double target_rate = well_controls_get_current_target(well_controls_);
                // TODO: the two situations below should be able to be merged to be handled as one situation

                if (number_phases_under_control == 1) { // single phase control
                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        if (distr[phase] > 0.) { // under the control of this phase
                            control_eq = getSegmentGTotal(0) * volumeFraction(0, flowPhaseToEbosCompIdx(phase)) - g[phase] * target_rate;
                            break;
                        }
                    }
                } else { // multiphase rate control
                    EvalWell rate_for_control(0.0);
                    const EvalWell G_total = getSegmentGTotal(0);
                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        if (distr[phase] > 0.) {
                            rate_for_control += G_total * volumeFractionScaled(0, flowPhaseToEbosCompIdx(phase));
                        }
                    }
                    // TODO: maybe the following equation can be scaled a little bit for gas phase
                    control_eq = rate_for_control - target_rate;
                }
                break;
            }
            case RESERVOIR_RATE:
            {
                EvalWell rate_for_control(0.0); // reservoir rate
                const double* distr = well_controls_get_current_distr(well_controls_);
                for (int phase = 0; phase < number_of_phases_; ++phase) {
                    if (distr[phase] > 0.) {
                        rate_for_control += getSegmentGTotal(0) * volumeFraction(0, flowPhaseToEbosCompIdx(phase));
                    }
                }
                const double target_rate = well_controls_get_current_target(well_controls_);
                control_eq = rate_for_control - target_rate;
                break;
            }
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
        }


        // using control_eq to update the matrix and residuals

        resWell_[0][SPres] = control_eq.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[0][0][SPres][pv_idx] = control_eq.derivative(pv_idx + numEq);
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assemblePressureEq(const int seg) const
    {
        assert(seg != 0); // not top segment

        // for top segment, the well control equation will be used.
        EvalWell pressure_equation = getSegmentPressure(seg);

        // we need to handle the pressure difference between the two segments
        // we only consider the hydrostatic pressure loss first
        pressure_equation -= getHydroPressureLoss(seg);

        if (frictionalPressureLossConsidered()) {
            pressure_equation -= getFrictionPressureLoss(seg);
        }

        resWell_[seg][SPres] = pressure_equation.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[seg][seg][SPres][pv_idx] = pressure_equation.derivative(pv_idx + numEq);
        }

        // contribution from the outlet segment
        const int outlet_segment_index = segmentNumberToIndex(segmentSet()[seg].outletSegment());
        const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

        resWell_[seg][SPres] -= outlet_pressure.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + numEq);
        }

        if (accelerationalPressureLossConsidered()) {
            handleAccelerationPressureLoss(seg);
        }
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getHydroPressureLoss(const int seg) const
    {
        return segment_densities_[seg] * gravity_ * segment_depth_diffs_[seg];
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getFrictionPressureLoss(const int seg) const
    {
        const EvalWell mass_rate = segment_mass_rates_[seg];
        const EvalWell density = segment_densities_[seg];
        const EvalWell visc = segment_viscosities_[seg];
        const int outlet_segment_index = segmentNumberToIndex(segmentSet()[seg].outletSegment());
        const double length = segmentSet()[seg].totalLength() - segmentSet()[outlet_segment_index].totalLength();
        assert(length > 0.);
        const double roughness = segmentSet()[seg].roughness();
        const double area = segmentSet()[seg].crossArea();
        const double diameter = segmentSet()[seg].internalDiameter();

        const double sign = mass_rate < 0. ? 1.0 : - 1.0;

        return sign * mswellhelpers::frictionPressureLoss(length, diameter, area, roughness, density, mass_rate, visc);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    handleAccelerationPressureLoss(const int seg) const
    {
        // TODO: this pressure loss is not significant enough to be well tested yet.
        // handle the out velcocity head
        const double area = segmentSet()[seg].crossArea();
        const EvalWell mass_rate = segment_mass_rates_[seg];
        const EvalWell density = segment_densities_[seg];
        const EvalWell out_velocity_head = mswellhelpers::velocityHead(area, mass_rate, density);

        resWell_[seg][SPres] -= out_velocity_head.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[seg][seg][SPres][pv_idx] -= out_velocity_head.derivative(pv_idx + numEq);
        }

        // calcuate the maximum cross-area among the segment and its inlet segments
        double max_area = area;
        for (const int inlet : segment_inlets_[seg]) {
            const double inlet_area = segmentSet()[inlet].crossArea();
            if (inlet_area > max_area) {
                max_area = inlet_area;
            }
        }

        // handling the velocity head of intlet segments
        for (const int inlet : segment_inlets_[seg]) {
            const EvalWell density = segment_densities_[inlet];
            const EvalWell mass_rate = segment_mass_rates_[inlet];
            const EvalWell inlet_velocity_head = mswellhelpers::velocityHead(area, mass_rate, density);
            resWell_[seg][SPres] += inlet_velocity_head.value();
            for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
                duneD_[seg][inlet][SPres][pv_idx] += inlet_velocity_head.derivative(pv_idx + numEq);
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    processFractions(const int seg) const
    {
        const PhaseUsage& pu = phaseUsage();

        std::vector<double> fractions(number_of_phases_, 0.0);

        assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
        const int oil_pos = pu.phase_pos[Oil];
        fractions[oil_pos] = 1.0;

        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            const int water_pos = pu.phase_pos[Water];
            fractions[water_pos] = primary_variables_[seg][WFrac];
            fractions[oil_pos] -= fractions[water_pos];
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            const int gas_pos = pu.phase_pos[Gas];
            fractions[gas_pos] = primary_variables_[seg][GFrac];
            fractions[oil_pos] -= fractions[gas_pos];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = pu.phase_pos[Water];
            if (fractions[water_pos] < 0.0) {
                if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
                    fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[water_pos]);
                }
                fractions[oil_pos] /= (1.0 - fractions[water_pos]);
                fractions[water_pos] = 0.0;
            }
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = pu.phase_pos[Gas];
            if (fractions[gas_pos] < 0.0) {
                if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
                    fractions[pu.phase_pos[Water]] /= (1.0 - fractions[gas_pos]);
                }
                fractions[oil_pos] /= (1.0 - fractions[gas_pos]);
                fractions[gas_pos] = 0.0;
            }
        }

        if (fractions[oil_pos] < 0.0) {
            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
                fractions[pu.phase_pos[Water]] /= (1.0 - fractions[oil_pos]);
            }
            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
                fractions[pu.phase_pos[Gas]] /= (1.0 - fractions[oil_pos]);
            }
            fractions[oil_pos] = 0.0;
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            primary_variables_[seg][WFrac] = fractions[pu.phase_pos[Water]];
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            primary_variables_[seg][GFrac] = fractions[pu.phase_pos[Gas]];
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkWellOperability(const Simulator& /* ebos_simulator */,
                         const WellState& /* well_state */,
                         Opm::DeferredLogger& deferred_logger)
    {
        const std::string msg = "Support of well operability checking for multisegment wells is not implemented "
                                "yet, checkWellOperability() for " + name() + " will do nothing";
        deferred_logger.warning("NO_OPERATABILITY_CHECKING_MS_WELLS", msg);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateFromPrimaryVariables(WellState& well_state) const
    {
        const PhaseUsage& pu = phaseUsage();
        assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
        const int oil_pos = pu.phase_pos[Oil];

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            std::vector<double> fractions(number_of_phases_, 0.0);
            fractions[oil_pos] = 1.0;

            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
                const int water_pos = pu.phase_pos[Water];
                fractions[water_pos] = primary_variables_[seg][WFrac];
                fractions[oil_pos] -= fractions[water_pos];
            }

            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
                const int gas_pos = pu.phase_pos[Gas];
                fractions[gas_pos] = primary_variables_[seg][GFrac];
                fractions[oil_pos] -= fractions[gas_pos];
            }

            // convert the fractions to be Q_p / G_total to calculate the phase rates
            for (int p = 0; p < number_of_phases_; ++p) {
                const double scale = scalingFactor(p);
                // for injection wells, there should only one non-zero scaling factor
                if (scale > 0.) {
                    fractions[p] /= scale;
                } else {
                    // this should only happens to injection wells
                    fractions[p] = 0.;
                }
            }

            // calculate the phase rates based on the primary variables
            const double g_total = primary_variables_[seg][GTotal];
            const int top_segment_index = well_state.topSegmentIndex(index_of_well_);
            for (int p = 0; p < number_of_phases_; ++p) {
                const double phase_rate = g_total * fractions[p];
                well_state.segRates()[(seg + top_segment_index) * number_of_phases_ + p] = phase_rate;
                if (seg == 0) { // top segment
                    well_state.wellRates()[index_of_well_ * number_of_phases_ + p] = phase_rate;
                }
            }

            // update the segment pressure
            well_state.segPress()[seg + top_segment_index] = primary_variables_[seg][SPres];
            if (seg == 0) { // top segment
                well_state.bhp()[index_of_well_] = well_state.segPress()[seg + top_segment_index];
            }
        }
    }




    template <typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    frictionalPressureLossConsidered() const
    {
        // HF- and HFA needs to consider frictional pressure loss
        return (segmentSet().compPressureDrop() != WellSegment::H__);
    }





    template <typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    accelerationalPressureLossConsidered() const
    {
        return (segmentSet().compPressureDrop() == WellSegment::HFA);
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    iterateWellEquations(const Simulator& ebosSimulator,
                         const std::vector<Scalar>& B_avg,
                         const double dt,
                         WellState& well_state,
                         Opm::DeferredLogger& deferred_logger)
    {
        const int max_iter_number = param_.max_inner_iter_ms_wells_;
        const WellState well_state0 = well_state;
        const std::vector<Scalar> residuals0 = getWellResiduals(B_avg);
        std::vector<std::vector<Scalar> > residual_history;
        std::vector<double> measure_history;
        int it = 0;
        // relaxation factor
        double relaxation_factor = 1.;
        const double min_relaxation_factor = 0.2;
        for (; it < max_iter_number; ++it) {

            assembleWellEqWithoutIteration(ebosSimulator, dt, well_state, deferred_logger);

            const BVectorWell dx_well = mswellhelpers::invDXDirect(duneD_, resWell_);


            const auto report = getWellConvergence(B_avg, deferred_logger);
            if (report.converged()) {
                break;
            }

            residual_history.push_back(getWellResiduals(B_avg));
            measure_history.push_back(getResidualMeasureValue(residual_history[it], deferred_logger) );

            bool is_oscillate = false;
            bool is_stagnate = false;

            detectOscillations(measure_history, it, is_oscillate, is_stagnate);
            // TODO: maybe we should have more sophiscated strategy to recover the relaxation factor,
            // for example, to recover it to be bigger

            if (is_oscillate || is_stagnate) {
                // a factor value to reduce the relaxation_factor
                const double reduction_mutliplier = 0.9;
                relaxation_factor = std::max(relaxation_factor * reduction_mutliplier, min_relaxation_factor);

                // debug output
                std::ostringstream sstr;
                if (is_stagnate) {
                    sstr << " well " << name() << " observes stagnation within " << it << "th inner iterations\n";
                }
                if (is_oscillate) {
                    sstr << " well " << name() << " osbserves oscillation within " << it <<"th inner iterations\n";
                }
                sstr << " relaxation_factor is " << relaxation_factor << " now\n";
                deferred_logger.debug(sstr.str());
            }

            updateWellState(dx_well, well_state, deferred_logger, relaxation_factor);

            // TODO: should we do something more if a switching of control happens
            this->updateWellControl(ebosSimulator, well_state, deferred_logger);

            initPrimaryVariablesEvaluation();
        }

        // TODO: we should decide whether to keep the updated well_state, or recover to use the old well_state
        if (it < max_iter_number) {
            std::ostringstream sstr;
            sstr << " well " << name() << " manage to get converged within " << it << " inner iterations";
            deferred_logger.debug(sstr.str());
        } else {
            std::ostringstream sstr;
            sstr << " well " << name() << " did not get converged within " << it << " inner iterations \n";
            sstr << " outputting the residual history for well " << name() << " during inner iterations \n";
            for (int i = 0; i < it; ++i) {
                const auto& residual = residual_history[i];
                sstr << " residual at " << i << "th iteration ";
                for (const auto& res : residual) {
                    sstr << " " << res;
                }
                sstr << " " << measure_history[i] << " \n";
            }
            deferred_logger.debug(sstr.str());
        }
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                   const double dt,
                                   WellState& well_state,
                                   Opm::DeferredLogger& deferred_logger)
    {
        // calculate the fluid properties needed.
        computeSegmentFluidProperties(ebosSimulator);

        // update the upwinding segments
        updateUpwindingSegments();

        // clear all entries
        duneB_ = 0.0;
        duneC_ = 0.0;

        duneD_ = 0.0;
        resWell_ = 0.0;

        well_state.wellVaporizedOilRates()[index_of_well_] = 0.;
        well_state.wellDissolvedGasRates()[index_of_well_] = 0.;

        // for the black oil cases, there will be four equations,
        // the first three of them are the mass balance equations, the last one is the pressure equations.
        //
        // but for the top segment, the pressure equation will be the well control equation, and the other three will be the same.

        const bool allow_cf = getAllowCrossFlow();

        const int nseg = numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the accumulation term
            // TODO: without considering the efficiencty factor for now
            {
                const EvalWell segment_surface_volume = getSegmentSurfaceVolume(ebosSimulator, seg);
                // for each component
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    const EvalWell accumulation_term = (segment_surface_volume * surfaceVolumeFraction(seg, comp_idx)
                                                     - segment_fluid_initial_[seg][comp_idx]) / dt;

                    resWell_[seg][comp_idx] += accumulation_term.value();
                    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
                        duneD_[seg][seg][comp_idx][pv_idx] += accumulation_term.derivative(pv_idx + numEq);
                    }
                }
            }
            // considering the contributions due to flowing out from the segment
            {
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    const EvalWell segment_rate = getSegmentRateUpwinding(seg, comp_idx);

                    const int seg_upwind = upwinding_segments_[seg];
                    // segment_rate contains the derivatives with respect to GTotal in seg,
                    // and WFrac and GFrac in seg_upwind
                    resWell_[seg][comp_idx] -= segment_rate.value();
                    duneD_[seg][seg][comp_idx][GTotal] -= segment_rate.derivative(GTotal + numEq);
                    duneD_[seg][seg_upwind][comp_idx][WFrac] -= segment_rate.derivative(WFrac + numEq);
                    duneD_[seg][seg_upwind][comp_idx][GFrac] -= segment_rate.derivative(GFrac + numEq);
                    // pressure derivative should be zero
                }
            }

            // considering the contributions from the inlet segments
            {
                for (const int inlet : segment_inlets_[seg]) {
                    for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                        const EvalWell inlet_rate = getSegmentRateUpwinding(inlet, comp_idx);

                        const int inlet_upwind = upwinding_segments_[inlet];
                        // inlet_rate contains the derivatives with respect to GTotal in inlet,
                        // and WFrac and GFrac in inlet_upwind
                        resWell_[seg][comp_idx] += inlet_rate.value();
                        duneD_[seg][inlet][comp_idx][GTotal] += inlet_rate.derivative(GTotal + numEq);
                        duneD_[seg][inlet_upwind][comp_idx][WFrac] += inlet_rate.derivative(WFrac + numEq);
                        duneD_[seg][inlet_upwind][comp_idx][GFrac] += inlet_rate.derivative(GFrac + numEq);
                        // pressure derivative should be zero
                    }
                }
            }

            // calculating the perforation rate for each perforation that belongs to this segment
            const EvalWell seg_pressure = getSegmentPressure(seg);
            for (const int perf : segment_perforations_[seg]) {
                const int cell_idx = well_cells_[perf];
                const auto& int_quants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                std::vector<EvalWell> mob(num_components_, 0.0);
                getMobility(ebosSimulator, perf, mob);
                std::vector<EvalWell> cq_s(num_components_, 0.0);
                EvalWell perf_press;
                double perf_dis_gas_rate = 0.;
                double perf_vap_oil_rate = 0.;

                computePerfRatePressure(int_quants, mob, seg, perf, seg_pressure, allow_cf, cq_s, perf_press, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);

                // updating the solution gas rate and solution oil rate
                if (well_type_ == PRODUCER) {
                    well_state.wellDissolvedGasRates()[index_of_well_] += perf_dis_gas_rate;
                    well_state.wellVaporizedOilRates()[index_of_well_] += perf_vap_oil_rate;
                }

                // store the perf pressure and rates
                const int rate_start_offset = (first_perf_ + perf) * number_of_phases_;
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    well_state.perfPhaseRates()[rate_start_offset + ebosCompIdxToFlowCompIdx(comp_idx)] = cq_s[comp_idx].value();
                }
                well_state.perfPress()[first_perf_ + perf] = perf_press.value();

                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    // the cq_s entering mass balance equations need to consider the efficiency factors.
                    const EvalWell cq_s_effective = cq_s[comp_idx] * well_efficiency_factor_;

                    connectionRates_[perf][comp_idx] = Base::restrictEval(cq_s_effective);

                    // subtract sum of phase fluxes in the well equations.
                    resWell_[seg][comp_idx] += cq_s_effective.value();

                    // assemble the jacobians
                    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {

                        // also need to consider the efficiency factor when manipulating the jacobians.
                        duneC_[seg][cell_idx][pv_idx][comp_idx] -= cq_s_effective.derivative(pv_idx + numEq); // intput in transformed matrix

                        // the index name for the D should be eq_idx / pv_idx
                        duneD_[seg][seg][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx + numEq);
                    }

                    for (int pv_idx = 0; pv_idx < numEq; ++pv_idx) {
                        // also need to consider the efficiency factor when manipulating the jacobians.
                        duneB_[seg][cell_idx][comp_idx][pv_idx] += cq_s_effective.derivative(pv_idx);
                    }
                }
            }

            // the fourth dequation, the pressure drop equation
            if (seg == 0) { // top segment, pressure equation is the control equation
                assembleControlEq(deferred_logger);
            } else {
                assemblePressureEq(seg);
            }
        }
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    wellTestingPhysical(const Simulator& /* simulator */, const std::vector<double>& /* B_avg */,
                        const double /* simulation_time */, const int /* report_step */,
                        WellState& /* well_state */, WellTestState& /* welltest_state */, Opm::DeferredLogger& deferred_logger)
    {
        const std::string msg = "Support of well testing for physical limits for multisegment wells is not "
                                "implemented yet, wellTestingPhysical() for " + name() + " will do nothing";
        deferred_logger.warning("NO_WELLTESTPHYSICAL_CHECKING_MS_WELLS", msg);
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWaterThroughput(const double dt OPM_UNUSED, WellState& well_state OPM_UNUSED) const
    {
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentSurfaceVolume(const Simulator& ebos_simulator, const int seg_idx) const
    {
        EvalWell temperature;
        int pvt_region_index;
        {
            // using the pvt region of first perforated cell
            // TODO: it should be a member of the WellInterface, initialized properly
            const int cell_idx = well_cells_[0];
            const auto& intQuants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            temperature.setValue(fs.temperature(FluidSystem::oilPhaseIdx).value());
            pvt_region_index = fs.pvtRegionIndex();
        }

        const EvalWell seg_pressure = getSegmentPressure(seg_idx);

        std::vector<EvalWell> mix_s(num_components_, 0.0);
        for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
            mix_s[comp_idx] = surfaceVolumeFraction(seg_idx, comp_idx);
        }

        std::vector<EvalWell> b(num_components_, 0.);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            b[waterCompIdx] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
        }

        EvalWell rv(0.0);
        // gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[oilCompIdx] > 0.0) {
                    if (mix_s[gasCompIdx] > 0.0) {
                        rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                    }

                    if (rv > rvmax) {
                        rv = rvmax;
                    }
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rv);
                } else { // no oil exists
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                }
            } else { // no Liquid phase
                // it is the same with zero mix_s[Oil]
                b[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
            }
        }

        EvalWell rs(0.0);
        // oil phase
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[gasCompIdx] > 0.0) {
                    if (mix_s[oilCompIdx] > 0.0) {
                        rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                    }
                    // std::cout << " rs " << rs.value() << " rsmax " << rsmax.value() << std::endl;

                    if (rs > rsmax) {
                        rs = rsmax;
                    }
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rs);
                } else { // no oil exists
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                }
            } else { // no gas phase
                // it is the same with zero mix_s[Gas]
                b[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
            }
        }

        std::vector<EvalWell> mix(mix_s);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

            const EvalWell d = 1.0 - rs * rv;
            if (d <= 0.0 || d > 1.0) {
                OPM_THROW(Opm::NumericalIssue, "Problematic d value " << d << " obtained for well " << name()
                                               << " during convertion to surface volume with rs " << rs
                                               << ", rv " << rv << " and pressure " << seg_pressure
                                               << " obtaining d " << d);
            }

            if (rs > 0.0) { // rs > 0.0?
                mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
            }
            if (rv > 0.0) { // rv > 0.0?
                mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
            }
        }

        EvalWell vol_ratio(0.0);
        for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
            vol_ratio += mix[comp_idx] / b[comp_idx];
        }

        // segment volume
        const double volume = segmentSet()[seg_idx].volume();

        return volume / vol_ratio;
    }





    template<typename TypeTag>
    std::vector<typename MultisegmentWell<TypeTag>::Scalar>
    MultisegmentWell<TypeTag>::
    getWellResiduals(const std::vector<Scalar>& B_avg) const
    {
        assert(int(B_avg.size() ) == num_components_);
        std::vector<Scalar> residuals(numWellEq + 1, 0.0);

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
                double residual = 0.;
                if (eq_idx < num_components_) {
                    residual = std::abs(resWell_[seg][eq_idx]) * B_avg[eq_idx];
                } else {
                    if (seg > 0) {
                        residual = std::abs(resWell_[seg][eq_idx]);
                    }
                }
                if (std::isnan(residual) || std::isinf(residual)) {
                    OPM_THROW(Opm::NumericalIssue, "nan or inf value for residal get for well " << name()
                                                    << " segment " << seg << " eq_idx " << eq_idx);
                }

                if (residual > residuals[eq_idx]) {
                    residuals[eq_idx] = residual;
                }
            }
        }

        // handling the control equation residual
        {
            const double control_residual = std::abs(resWell_[0][numWellEq - 1]);
            if (std::isnan(control_residual) || std::isinf(control_residual)) {
               OPM_THROW(Opm::NumericalIssue, "nan or inf value for control residal get for well " << name());
            }
            residuals[numWellEq] = control_residual;
        }

        return residuals;
    }





    /// Detect oscillation or stagnation based on the residual measure history
    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    detectOscillations(const std::vector<double>& measure_history,
                       const int it, bool& oscillate, bool& stagnate) const
    {
        if ( it < 2 ) {
            oscillate = false;
            stagnate = false;
            return;
        }

        stagnate = true;
        const double F0 = measure_history[it];
        const double F1 = measure_history[it - 1];
        const double F2 = measure_history[it - 2];
        const double d1 = std::abs((F0 - F2) / F0);
        const double d2 = std::abs((F0 - F1) / F0);

        const double oscillaton_rel_tol = 0.2;
        oscillate = (d1 < oscillaton_rel_tol) && (oscillaton_rel_tol < d2);

        const double stagnation_rel_tol = 1.e-2;
        stagnate = std::abs((F1 - F2) / F2) <= stagnation_rel_tol;
    }





    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    getResidualMeasureValue(const std::vector<double>& residuals,
                            DeferredLogger& deferred_logger) const
    {
        assert(int(residuals.size()) == numWellEq + 1);

        const double rate_tolerance = param_.tolerance_wells_;
        int count = 0;
        double sum = 0;
        for (int eq_idx = 0; eq_idx < numWellEq - 1; ++eq_idx) {
            if (residuals[eq_idx] > rate_tolerance) {
                sum += residuals[eq_idx] / rate_tolerance;
                ++count;
            }
        }

        const double pressure_tolerance = param_.tolerance_pressure_ms_wells_;
        if (residuals[SPres] > pressure_tolerance) {
            sum += residuals[SPres] / pressure_tolerance;
            ++count;
        }

        const double control_tolerance = getControlTolerance(deferred_logger);
        if (residuals[SPres + 1] > control_tolerance) {
            sum += residuals[SPres + 1] / control_tolerance;
            ++count;
        }

        // if (count == 0), it should be converged.
        assert(count != 0);

        return sum;
    }





    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    getControlTolerance(DeferredLogger& deferred_logger) const
    {
        double control_tolerance = 0.;
        switch(well_controls_get_current_type(well_controls_) ) {
            case BHP:
                control_tolerance = param_.tolerance_wells_;
                break;
            case THP:
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case RESERVOIR_RATE:
            case SURFACE_RATE:
                control_tolerance = param_.tolerance_wells_;
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
        }

        return control_tolerance;
    }




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkConvergenceControlEq(ConvergenceReport& report,
                              DeferredLogger& deferred_logger) const
    {
        using CR = ConvergenceReport;
        CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;
        switch(well_controls_get_current_type(well_controls_)) {
            case THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                break;
            case BHP:
                ctrltype = CR::WellFailure::Type::ControlBHP;
                break;
            case RESERVOIR_RATE:
            case SURFACE_RATE:
                ctrltype = CR::WellFailure::Type::ControlRate;
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
        }
        assert(ctrltype != CR::WellFailure::Type::Invalid);

        const double control_residual = std::abs(resWell_[0][SPres]);
        const double control_tolerance = getControlTolerance(deferred_logger);

        const int dummy_component = -1;
        if (std::isnan(control_residual)) {
            report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, name()});
        } else if (control_residual > param_.max_residual_allowed_) {
            report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, name()});
        } else if (control_residual > control_tolerance) {
            report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, name()});
        }
    }






    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateUpwindingSegments()
    {
        // not considering upwinding for the injectors for now
        // but we should
        // and upwinding segment for top segment is itself
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            if ( (well_type_ == INJECTOR) || (seg == 0) ) {
                upwinding_segments_[seg] = seg;
                continue;
            }

            // for other normal segments
            if (primary_variables_evaluation_[seg][GTotal] <= 0.) {
                upwinding_segments_[seg] = seg;
            } else {
                const int outlet_segment_index = segmentNumberToIndex(segmentSet()[seg].outletSegment());
                upwinding_segments_[seg] = outlet_segment_index;
            }
        }
    }

}
