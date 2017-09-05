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




namespace Opm
{

    template<typename TypeTag>
    MultisegmentWell<TypeTag>::
    MultisegmentWell(const Well* well, const int time_step, const Wells* wells)
    : Base(well, time_step, wells)
    , segment_cell_(numberOfSegments())
    , segment_perforations_(numberOfSegments())
    , segment_inlets_(numberOfSegments())
    {
        // TODO: to see what information we need to process here later.
        // const auto& completion_set = well->getCompletions(time_step);
        // const auto& segment_set = well->getSegmentSet(time_step);

        // since we decide to use the SegmentSet from the well parser. we can reuse a lot from it.
        // other facilities needed we need to process them here

    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<bool>* active_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells)
    {
        Base::init(phase_usage_arg, active_arg,
                   depth_arg, gravity_arg, num_cells);

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

        // initialize the segment_perforations_
        const CompletionSet& completion_set = well_ecl_->getCompletions(current_step_);
        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const Completion& completion = completion_set.get(perf);
            const int segment_number = completion.getSegmentNumber();
            const int segment_location = numberToLocation(segment_number);
            segment_perforations_[segment_location].push_back(perf);
        }

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            const Segment& segment = segmentSet()[seg];
            const int segment_number = segment.segmentNumber();
            const int outlet_segment_number = segment.outletSegment();
            if (outlet_segment_number > 0) {
                const int segment_location = numberToLocation(segment_number);
                const int outlet_segment_location = numberToLocation(outlet_segment_number);
                segment_inlets_[outlet_segment_location].push_back(segment_location);
            }
        }
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initMatrixAndVectors(const int num_cells) const
    {
        duneB_.setBuildMode( OffDiagMatWell::row_wise );
        duneC_.setBuildMode( OffDiagMatWell::row_wise );
        invDuneD_.setBuildMode( DiagMatWell::row_wise );

        // set the size and patterns for all the matrices and vectors
        // [A C^T    [x    =  [ res
        //  B D] x_well]      res_well]

        // the number of the nnz should be numSegment() + numberOfOutlet()
        invDuneD_.setSize(numberOfSegments(), numberOfSegments(), 100000);
        duneB_.setSize(numberOfSegments(), num_cells, number_of_perforations_);
        duneC_.setSize(numberOfSegments(), num_cells, number_of_perforations_);

        // we need to add the off diagonal ones
        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (int perf = 0 ; perf < number_of_perforations_; ++perf) { // the segments hold some perforations
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        // make the B^T matrix
        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize( numberOfSegments() );

        // resize temporary class variables
        Bx_.resize( duneC_.N() );
        invDrw_.resize( invDuneD_.N() );

        // TODO: maybe this function need a different name for better meaning
        primary_variables_.resize(numberOfSegments());
        primary_variables_evaluation_.resize(numberOfSegments());
    }





    template<typename TypeTag>
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





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleWellEq(Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells)
    {
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateWithTarget(const int current,
                              WellState& well_state) const
    {
        // TODO: it can be challenging, when updating the segment and perforation related,
        // well rates needs to be okay.

        // Updating well state bas on well control
        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
        const double target = well_controls_iget_target(well_controls_, current);
        const double* distr = well_controls_iget_distr(well_controls_, current);
        switch (well_controls_iget_type(well_controls_, current)) {
        case BHP: {
            well_state.bhp()[index_of_well_] = target;
            const int top_segment_location = well_state.topSegmentLocation(index_of_well_);
            well_state.segPress()[top_segment_location] = well_state.bhp()[index_of_well_];
            // TODO: similar to the way below to handle THP
            // we should not something related to thp here when there is thp constraint
            break;
        }

        case THP: {
            well_state.thp()[index_of_well_] = target;

            const Opm::PhaseUsage& pu = phaseUsage();
            std::vector<double> rates(3, 0.0);
            if (active()[ Water ]) {
                rates[ Water ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Water ] ];
            }
            if (active()[ Oil ]) {
                 rates[ Oil ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Oil ] ];
            }
            if (active()[ Gas ]) {
                rates[ Gas ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Gas ] ];
            }

            const int table_id = well_controls_iget_vfp(well_controls_, current);
            const double& thp    = well_controls_iget_target(well_controls_, current);
            const double& alq    = well_controls_iget_alq(well_controls_, current);

            // TODO: implement calculateBhpFromThp function
            // well_state.bhp()[index_of_well_] = calculateBhpFromThp(rates, current);
            // also the top segment pressure
            /* const int top_segment_location = well_state.topSegmentLocation(index_of_well_);
            well_state.segPress()[top_segment_location] = well_state.bhp()[index_of_well_]; */
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
                    }
                } else { // scaling factor is not well defied when original_rates_under_phase_control is zero
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
                }
            }

            // update the perforation rates and then segment rates
            {
                for (int phase = 0; phase < number_of_phases_; ++phase) {
                    const double perf_phaserate = well_state.wellRates()[number_of_phases_ * index_of_well_ + phase] / number_of_perforations_;
                    for (int perf = 0; perf < number_of_perforations_; ++perf) {
                        well_state.perfPhaseRates()[number_of_phases_ * (first_perf_ + perf) + phase] = perf_phaserate;
                    }
                }

                const std::vector<double> perforation_rates(well_state.perfPhaseRates().begin() + number_of_phases_ * first_perf_,
                                             well_state.perfPhaseRates().begin() + number_of_phases_ * (first_perf_ + number_of_perforations_) );
                std::vector<double> segment_rates;
                WellState::calculateSegmentRates(segment_inlets_, segment_perforations_, perforation_rates, number_of_phases_,
                                                 0 /* top segment */, segment_rates);
                const int top_segment_location = well_state.topSegmentLocation(index_of_well_);
                std::copy(segment_rates.begin(), segment_rates.end(),
                          well_state.segRates().begin() + number_of_phases_ * top_segment_location );
                // we need to check the top segment rates should be same with the well rates
            }

            break;
        } // end of switch

        updatePrimaryVariables(well_state);
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::ConvergenceReport
    MultisegmentWell<TypeTag>::
    getWellConvergence(Simulator& ebosSimulator,
                       const std::vector<double>& B_avg,
                       const ModelParameters& param) const
    {
        // TODO: it will be very similar
        ConvergenceReport report;
        return report;
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeAccumWell()
    {
        // it will be vector of compositions of segments
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellConnectionPressures(const Simulator& ebosSimulator,
                                   const WellState& well_state)
    {
        // TODO: the name of the function need to change.
        // it will be calculating the pressure difference between the perforation and grid cells
        // With MS well, the depth of the perforation is not necessarily the center of the grid cells.
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {


    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(BVector& r) const
    {


    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x,
                                          const ModelParameters& param,
                                          WellState& well_state) const
    {



    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials) const
    {



    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state) const
    {
        // TODO: not tested yet.
        // TODO: not handling solvent or polymer for now.

        // TODO: to test using rate conversion coefficients to see if it will be better than
        // this default one
        std::vector<double> g = {1.0, 1.0, 0.01};

        if (well_controls_get_current_type(well_controls_) == RESERVOIR_RATE) {
            const double* distr = well_controls_get_current_distr(well_controls_);
            for (int phase = 0; phase < number_of_phases_; ++phase) {
                g[phase] = distr[phase];
            }
        }

        // the location of the top segment in the WellState
        const int top_segment_location = well_state.topSegmentLocation(index_of_well_);
        const std::vector<double>& segment_rates = well_state.segRates();
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // calculate the total rate for each segment
            double total_seg_rate = 0.0;
            const int seg_location = top_segment_location + seg;
            // the segment pressure
            primary_variables_[seg][SPres] = well_state.segPress()[seg_location];
            // TODO: under what kind of circustances, the following will be wrong?
            // the definition of g makes the gas phase is always the last phase
            for (int p = 0; p < number_of_phases_; p++) {
                total_seg_rate += g[p] * segment_rates[number_of_phases_ * seg_location + p];
            }

            primary_variables_[seg][GTotal] = total_seg_rate;
            if (std::abs(total_seg_rate) > 0.) {
                if (active()[Water]) {
                    primary_variables_[seg][WFrac] = g[Water] * segment_rates[number_of_phases_ * seg_location + Water] / total_seg_rate;
                }
                if (active()[Gas]) {
                    primary_variables_[seg][GFrac] = g[Gas] * segment_rates[number_of_phases_ * seg_location + Gas] / total_seg_rate;
                }
            } else { // total_seg_rate == 0
                if (well_type_ == INJECTOR) {
                    // only single phase injection handled
                    const double* distr = well_controls_get_current_distr(well_controls_);
                    if (active()[Water]) {
                        if (distr[Water] > 0.0) {
                            primary_variables_[seg][WFrac] = 1.0;
                        } else {
                            primary_variables_[seg][WFrac] = 0.0;
                        }
                    }

                    if (active()[Gas]) {
                        if (distr[Gas] > 0.0) {
                            // TODO: not handling solvent here yet
                            primary_variables_[seg][GFrac] = 1.0;
                        } else {
                            primary_variables_[seg][GFrac] = 0.0;
                        }
                    }
                } else if (well_type_ == PRODUCER) { // producers
                    if (active()[Water]) {
                        primary_variables_[seg][WFrac] = 1.0 / number_of_phases_;
                    }

                    if (active()[Gas]) {
                        primary_variables_[seg][GFrac] = 1.0 / number_of_phases_;
                    }
                }
            }
        }
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    solveEqAndUpdateWellState(const ModelParameters& param,
                              WellState& well_state)
    {


    }





    template<typename TypeTag>
    const SegmentSet&
    MultisegmentWell<TypeTag>::
    segmentSet() const
    {
        return well_ecl_->getSegmentSet(current_step_);
    }





    template<typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    numberOfSegments() const
    {
        return segmentSet().numberSegment();
    }





    template<typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    numberOfPerforations() const
    {
        return segmentSet().number_of_perforations_;
    }





    template<typename TypeTag>
    WellSegment::CompPressureDropEnum
    MultisegmentWell<TypeTag>::
    compPressureDrop() const
    {
        return segmentSet().compPressureDrop();
    }





    template<typename TypeTag>
    WellSegment::MultiPhaseModelEnum
    MultisegmentWell<TypeTag>::
    multiphaseModel() const
    {
        return segmentSet().multiPhaseModel();
    }





    template<typename TypeTag>
    int
    MultisegmentWell<TypeTag>::
    numberToLocation(const int segment_number) const
    {
        return segmentSet().numberToLocation(segment_number);
    }

}
