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


#include <opm/autodiff/MSWellHelpers.hpp>

namespace Opm
{


    template <typename TypeTag>
    MultisegmentWell<TypeTag>::
    MultisegmentWell(const Well* well, const int time_step, const Wells* wells,
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
    , segment_comp_initial_(numberOfSegments(), std::vector<double>(num_components_, 0.0))
    , segment_densities_(numberOfSegments(), 0.0)
    , segment_viscosities_(numberOfSegments(), 0.0)
    , segment_mass_rates_(numberOfSegments(), 0.0)
    , segment_depth_diffs_(numberOfSegments(), 0.0)
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

        // initialize the segment_perforations_
        const WellConnections& completion_set = well_ecl_->getConnections(current_step_);
        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const Connection& connection = completion_set.get(perf);
            const int segment_index = segmentNumberToIndex(connection.segment());
            segment_perforations_[segment_index].push_back(perf);
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

        // callcuate the depth difference between perforations and their segments
        perf_depth_.resize(number_of_perforations_, 0.);
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            const double segment_depth = segmentSet()[seg].depth();
            for (const int perf : segment_perforations_[seg]) {
                perf_depth_[perf] = completion_set.get(perf).depth();
                perforation_segment_depth_diffs_[perf] = perf_depth_[perf] - segment_depth;
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
                   const double dt,
                   WellState& well_state,
                   Opm::DeferredLogger& deferred_logger)
    {

        const bool use_inner_iterations = param_.use_inner_iterations_ms_wells_;
        if (use_inner_iterations) {
            iterateWellEquations(ebosSimulator, dt, well_state);
        }

        assembleWellEqWithoutIteration(ebosSimulator, dt, well_state);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              WellState& well_state,
                              Opm::DeferredLogger& deferred_logger) const
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
    getWellConvergence(const std::vector<double>& B_avg) const
    {
        assert(int(B_avg.size()) == num_components_);

        // checking if any residual is NaN or too large. The two large one is only handled for the well flux
        std::vector<std::vector<double>> abs_residual(numberOfSegments(), std::vector<double>(numWellEq, 0.0));
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
                abs_residual[seg][eq_idx] = std::abs(resWell_[seg][eq_idx]);
            }
        }

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
                OPM_THROW(std::runtime_error, "Unknown well control control types for well " << name());
        }
        assert(ctrltype != CR::WellFailure::Type::Invalid);

        std::vector<double> maximum_residual(numWellEq, 0.0);

        ConvergenceReport report;
        const int dummy_component = -1;
        // TODO: the following is a little complicated, maybe can be simplified in some way?
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                if (eq_idx < num_components_) { // phase or component mass equations
                    const double flux_residual = B_avg[eq_idx] * abs_residual[seg][eq_idx];
                    if (flux_residual > maximum_residual[eq_idx]) {
                        maximum_residual[eq_idx] = flux_residual;
                    }
                } else { // pressure or control equation
                    if (seg == 0) {
                        // Control equation
                        const double control_residual = abs_residual[seg][eq_idx];
                        if (std::isnan(control_residual)) {
                            report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, name()});
                        } else if (control_residual > param_.max_residual_allowed_) {
                            report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, name()});
                        } else if (control_residual > param_.tolerance_wells_) {
                            report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, name()});
                        }
                    } else {
                        // Pressure equation
                        const double pressure_residual = abs_residual[seg][eq_idx];
                        if (pressure_residual > maximum_residual[eq_idx]) {
                            maximum_residual[eq_idx] = pressure_residual;
                        }
                    }
                }
            }
        }

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
                if (std::isnan(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::NotANumber, dummy_component, name()});
                } else if (std::isinf(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::TooLarge, dummy_component, name()});
                } else if (pressure_residual > param_.tolerance_pressure_ms_wells_) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, name()});
                }
            }
        }

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
                                          WellState& well_state) const
    {
        BVectorWell xw(1);
        recoverSolutionWell(x, xw);
        updateWellState(xw, false, well_state);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellPotentials(const Simulator& /* ebosSimulator */,
                          const WellState& /* well_state */,
                          std::vector<double>& well_potentials,
                          Opm::DeferredLogger& deferred_logger)
    {
        const std::string msg = std::string("Well potential calculation is not supported for multisegment wells \n")
                + "A well potential of zero is returned for output purposes. \n"
                + "If you need well potential to set the guide rate for group controled wells \n"
                + "you will have to change the " + name() + " well to a standard well \n";

        deferred_logger.warning("WELL_POTENTIAL_NOT_IMPLEMENTED_FOR_MULTISEG_WELLS", msg);

        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state) const
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
    solveEqAndUpdateWellState(WellState& well_state)
    {
        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        const BVectorWell dx_well = mswellhelpers::invDXDirect(duneD_, resWell_);

        updateWellState(dx_well, false, well_state);
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
    computeInitialComposition()
    {
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // TODO: probably it should be numWellEq -1 more accurately,
            // while by meaning it should be num_comp
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                segment_comp_initial_[seg][comp_idx] = surfaceVolumeFraction(seg, comp_idx).value();
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellState(const BVectorWell& dwells,
                    const bool inner_iteration,
                    WellState& well_state) const
    {
        const bool use_inner_iterations = param_.use_inner_iterations_ms_wells_;

        const double relaxation_factor = (use_inner_iterations && inner_iteration) ? 0.2 : 1.0;

        const double dFLimit = param_.dwell_fraction_max_;
        const double max_pressure_change = param_.max_pressure_change_ms_wells_;
        const std::vector<std::array<double, numWellEq> > old_primary_variables = primary_variables_;

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]), relaxation_factor * dFLimit);
                primary_variables_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]), relaxation_factor * dFLimit);
                primary_variables_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
            }

            // handling the overshooting or undershooting of the fractions
            processFractions(seg);

            // update the segment pressure
            {
                const int sign = dwells[seg][SPres] > 0.? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]), relaxation_factor * max_pressure_change);
                primary_variables_[seg][SPres] = old_primary_variables[seg][SPres] - dx_limited;
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
                                const WellState& /* well_state */)
    {
        computePerfCellPressDiffs(ebosSimulator);
        computeInitialComposition();
    }





    template <typename TypeTag>
    const WellSegments&
    MultisegmentWell<TypeTag>::
    segmentSet() const
    {
        return well_ecl_->getWellSegments(current_step_);
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
    computePerfRate(const IntensiveQuantities& int_quants,
                    const std::vector<EvalWell>& mob_perfcells,
                    const int seg,
                    const int perf,
                    const EvalWell& segment_pressure,
                    const bool& allow_cf,
                    std::vector<EvalWell>& cq_s) const
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

        // Pressure drawdown (also used to determine direction of flow)
        // TODO: not 100% sure about the sign of the seg_perf_press_diff
        const EvalWell drawdown = (pressure_cell + cell_perf_press_diff) - (segment_pressure + perf_seg_press_diff);

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
                    OPM_THROW(Opm::NumericalIssue, "Zero d value obtained for well " << name() << " during flux calcuation"
                                                  << " with rs " << rs << " and rv " << rv);
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
    assembleControlEq() const
    {
        EvalWell control_eq(0.0);

        switch (well_controls_get_current_type(well_controls_)) {
            case THP: // not handling this one for now
            {
                OPM_THROW(std::runtime_error, "Not handling THP control for Multisegment wells for now");
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
                OPM_THROW(std::runtime_error, "Unknown well control control types for well " << name());
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
                         const double dt,
                         WellState& well_state)
    {
        // basically, it only iterate through the equations.
        // we update the primary variables
        // if converged, we can update the well_state.
        // the function updateWellState() should have a flag to show
        // if we will update the well state.
        const int max_iter_number = param_.max_inner_iter_ms_wells_;
        int it = 0;
        for (; it < max_iter_number; ++it) {

            assembleWellEqWithoutIteration(ebosSimulator, dt, well_state);

            const BVectorWell dx_well = mswellhelpers::invDXDirect(duneD_, resWell_);

            // TODO: use these small values for now, not intend to reach the convergence
            // in this stage, but, should we?
            // We should try to avoid hard-code values in the code.
            // If we want to use the real one, we need to find a way to get them.
            // const std::vector<double> B {0.8, 0.8, 0.008};
            const std::vector<double> B {0.5, 0.5, 0.005};

            const auto report = getWellConvergence(B);
            if (report.converged()) {
                break;
            }

            updateWellState(dx_well, true, well_state);

            initPrimaryVariablesEvaluation();
        }
        // TODO: maybe we should not use these values if they are not converged.
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                   const double dt,
                                   WellState& well_state)
    {
        // calculate the fluid properties needed.
        computeSegmentFluidProperties(ebosSimulator);

        // clear all entries
        duneB_ = 0.0;
        duneC_ = 0.0;

        duneD_ = 0.0;
        resWell_ = 0.0;

        // for the black oil cases, there will be four equations,
        // the first three of them are the mass balance equations, the last one is the pressure equations.
        //
        // but for the top segment, the pressure equation will be the well control equation, and the other three will be the same.

        const bool allow_cf = getAllowCrossFlow();

        const int nseg = numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the accumulation term // TODO: without considering the efficiencty factor for now
            // volume of the segment
            {
                const double volume = segmentSet()[seg].volume();
                // for each component
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    EvalWell accumulation_term = volume / dt * (surfaceVolumeFraction(seg, comp_idx) - segment_comp_initial_[seg][comp_idx])
                                               + getSegmentRate(seg, comp_idx);

                    resWell_[seg][comp_idx] += accumulation_term.value();
                    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
                        duneD_[seg][seg][comp_idx][pv_idx] += accumulation_term.derivative(pv_idx + numEq);
                    }
                }
            }

            // considering the contributions from the inlet segments
            {
                for (const int inlet : segment_inlets_[seg]) {
                    for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                        const EvalWell inlet_rate = getSegmentRate(inlet, comp_idx);
                        resWell_[seg][comp_idx] -= inlet_rate.value();
                        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
                            duneD_[seg][inlet][comp_idx][pv_idx] -= inlet_rate.derivative(pv_idx + numEq);
                        }
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
                computePerfRate(int_quants, mob, seg, perf, seg_pressure, allow_cf, cq_s);

                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    // the cq_s entering mass balance equations need to consider the efficiency factors.
                    const EvalWell cq_s_effective = cq_s[comp_idx] * well_efficiency_factor_;

                    connectionRates_[perf][comp_idx] = Base::restrictEval(cq_s_effective);

                    // subtract sum of phase fluxes in the well equations.
                    resWell_[seg][comp_idx] -= cq_s_effective.value();

                    // assemble the jacobians
                    for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {

                        // also need to consider the efficiency factor when manipulating the jacobians.
                        duneC_[seg][cell_idx][pv_idx][comp_idx] -= cq_s_effective.derivative(pv_idx + numEq); // intput in transformed matrix
                        
                        // the index name for the D should be eq_idx / pv_idx
                        duneD_[seg][seg][comp_idx][pv_idx] -= cq_s_effective.derivative(pv_idx + numEq);
                    }

                    for (int pv_idx = 0; pv_idx < numEq; ++pv_idx) {
                        // also need to consider the efficiency factor when manipulating the jacobians.
                        duneB_[seg][cell_idx][comp_idx][pv_idx] -= cq_s_effective.derivative(pv_idx);   
                    }
                }
                // TODO: we should save the perforation pressure and preforation rates?
                // we do not use it in the simulation for now, while we might need them if
                // we handle the pressure in SEG mode.
            }

            // the fourth dequation, the pressure drop equation
            if (seg == 0) { // top segment, pressure equation is the control equation
                assembleControlEq();
            } else {
                assemblePressureEq(seg);
            }
        }
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    wellTestingPhysical(Simulator& simulator, const std::vector<double>& B_avg,
                        const double simulation_time, const int report_step,
                        WellState& well_state, WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger)
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

}
