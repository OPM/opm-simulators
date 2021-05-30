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
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/Valve.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <string>
#include <algorithm>

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
                     const int first_perf_index,
                     const std::vector<PerforationData>& perf_data)
    : Base(well, pw_info, time_step, param, rate_converter, pvtRegionIdx, num_components, num_phases, index_of_well, first_perf_index, perf_data)
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
    , segment_phase_fractions_(numberOfSegments(), std::vector<EvalWell>(num_components_, 0.0)) // number of phase here?
    , segment_phase_viscosities_(numberOfSegments(), std::vector<EvalWell>(num_components_, 0.0)) // number of phase here?
    , segment_phase_densities_(numberOfSegments(), std::vector<EvalWell>(num_components_, 0.0)) // number of phase here?
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
            if (connection.state() == Connection::State::OPEN) {
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
         const int num_cells,
         const std::vector< Scalar >& B_avg)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, num_cells, B_avg);

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
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              WellState& well_state,
                              DeferredLogger&  deferred_logger) const
    {
        Base::updateWellStateWithTarget(ebos_simulator, well_state, deferred_logger);
        // scale segment rates based on the wellRates
        // and segment pressure based on bhp
        scaleSegmentPressuresWithBhp(well_state);
        scaleSegmentRatesWithWellRates(well_state);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    scaleSegmentRatesWithWellRates(WellState& well_state) const
    {
        auto segment_rates = well_state.segRates(index_of_well_);
        for (int phase = 0; phase < number_of_phases_; ++phase) {
            const double unscaled_top_seg_rate = segment_rates[phase];
            const double well_phase_rate = well_state.wellRates(index_of_well_)[phase];
            if (std::abs(unscaled_top_seg_rate) > 1e-12)
            {
                for (int seg = 0; seg < numberOfSegments(); ++seg) {
                    segment_rates[this->number_of_phases_*seg + phase] *= well_phase_rate/unscaled_top_seg_rate;
                }
            } else {
                // for newly opened wells, the unscaled rate top segment rate is zero
                // and we need to initialize the segment rates differently
                double sumTw = 0;
                for (int perf = 0; perf < number_of_perforations_; ++perf) {
                    sumTw += well_index_[perf];
                }

                std::vector<double> perforation_rates(number_of_phases_ * number_of_perforations_,0.0);
                const double perf_phaserate_scaled = well_state.wellRates(index_of_well_)[phase] / sumTw;
                for (int perf = 0; perf < number_of_perforations_; ++perf) {
                    perforation_rates[number_of_phases_ * perf + phase] = well_index_[perf] * perf_phaserate_scaled;
                }

                std::vector<double> rates;
                WellState::calculateSegmentRates(segment_inlets_, segment_perforations_, perforation_rates, number_of_phases_,
                                                 0, rates);
                std::copy(rates.begin(), rates.end(), segment_rates);
            }
        }
    }

    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    scaleSegmentPressuresWithBhp(WellState& well_state) const
    {
        //scale segment pressures
        const double bhp = well_state.bhp(index_of_well_);
        auto segment_pressure = well_state.segPress(index_of_well_);
        const double unscaled_top_seg_pressure = segment_pressure[0];
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            segment_pressure[seg] *= bhp/unscaled_top_seg_pressure;
        }
    }


    template <typename TypeTag>
    ConvergenceReport
    MultisegmentWell<TypeTag>::
    getWellConvergence(const WellState& well_state, const std::vector<double>& B_avg, DeferredLogger& deferred_logger, const bool relax_tolerance) const
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
                } else if (!relax_tolerance && flux_residual > param_.tolerance_wells_) {
                    report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, eq_idx, name()});
                } else if (flux_residual > param_.relaxed_inner_tolerance_flow_ms_well_) {
                    report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, eq_idx, name()});
                }
            } else { // pressure equation
                const double pressure_residual = maximum_residual[eq_idx];
                const int dummy_component = -1;
                if (std::isnan(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::NotANumber, dummy_component, name()});
                } else if (std::isinf(pressure_residual)) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::TooLarge, dummy_component, name()});
                } else if (!relax_tolerance && pressure_residual > param_.tolerance_pressure_ms_wells_) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, name()});
                } else if (pressure_residual > param_.relaxed_inner_tolerance_pressure_ms_well_) {
                    report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy_component, name()});
                }
            }
        }

        checkConvergenceControlEq(well_state, report, deferred_logger);

        return report;
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        if ( param_.matrix_add_well_contributions_ )
        {
            // Contributions are already in the matrix itself
            return;
        }
        BVectorWell Bx(duneB_.N());

        duneB_.mv(x, Bx);

        // invDBx = duneD^-1 * Bx_
        const BVectorWell invDBx = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, Bx);

        // Ax = Ax - duneC_^T * invDBx
        duneC_.mmtv(invDBx,Ax);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    apply(BVector& r) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        // invDrw_ = duneD^-1 * resWell_
        const BVectorWell invDrw = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell_);
        // r = r - duneC_^T * invDrw
        duneC_.mmtv(invDrw, r);
    }



#if HAVE_CUDA || HAVE_OPENCL
    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    addWellContribution(WellContributions& wellContribs) const
    {
        unsigned int Mb = duneB_.N();       // number of blockrows in duneB_, duneC_ and duneD_
        unsigned int BnumBlocks = duneB_.nonzeroes();
        unsigned int DnumBlocks = duneD_.nonzeroes();

        // duneC
        std::vector<unsigned int> Ccols;
        std::vector<double> Cvals;
        Ccols.reserve(BnumBlocks);
        Cvals.reserve(BnumBlocks * numEq * numWellEq);
        for (auto rowC = duneC_.begin(); rowC != duneC_.end(); ++rowC) {
            for (auto colC = rowC->begin(), endC = rowC->end(); colC != endC; ++colC) {
                Ccols.emplace_back(colC.index());
                for (int i = 0; i < numWellEq; ++i) {
                    for (int j = 0; j < numEq; ++j) {
                        Cvals.emplace_back((*colC)[i][j]);
                    }
                }
            }
        }

        // duneD
        Dune::UMFPack<DiagMatWell> umfpackMatrix(duneD_, 0);
        double *Dvals = umfpackMatrix.getInternalMatrix().getValues();
        auto *Dcols = umfpackMatrix.getInternalMatrix().getColStart();
        auto *Drows = umfpackMatrix.getInternalMatrix().getRowIndex();

        // duneB
        std::vector<unsigned int> Bcols;
        std::vector<unsigned int> Brows;
        std::vector<double> Bvals;
        Bcols.reserve(BnumBlocks);
        Brows.reserve(Mb+1);
        Bvals.reserve(BnumBlocks * numEq * numWellEq);
        Brows.emplace_back(0);
        unsigned int sumBlocks = 0;
        for (auto rowB = duneB_.begin(); rowB != duneB_.end(); ++rowB) {
            int sizeRow = 0;
            for (auto colB = rowB->begin(), endB = rowB->end(); colB != endB; ++colB) {
                Bcols.emplace_back(colB.index());
                for (int i = 0; i < numWellEq; ++i) {
                    for (int j = 0; j < numEq; ++j) {
                        Bvals.emplace_back((*colB)[i][j]);
                    }
                }
                sizeRow++;
            }
            sumBlocks += sizeRow;
            Brows.emplace_back(sumBlocks);
        }

        wellContribs.addMultisegmentWellContribution(numEq, numWellEq, Mb, Bvals, Bcols, Brows, DnumBlocks, Dvals, Dcols, Drows, Cvals);
    }
#endif


    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x,
                                          WellState& well_state,
                                          DeferredLogger& deferred_logger) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        BVectorWell xw(1);
        recoverSolutionWell(x, xw);
        updateWellState(xw, well_state, deferred_logger);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials,
                          DeferredLogger& deferred_logger)
    {
        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

        // Stopped wells have zero potential.
        if (this->wellIsStopped()) {
            return;
        }

        // If the well is pressure controlled the potential equals the rate.
        bool pressure_controlled_well = false;
        if (this->isInjector()) {
            const Well::InjectorCMode& current = well_state.currentInjectionControl(index_of_well_);
            if (current == Well::InjectorCMode::BHP || current == Well::InjectorCMode::THP) {
                pressure_controlled_well = true;
            }
        } else {
            const Well::ProducerCMode& current = well_state.currentProductionControl(index_of_well_);
            if (current == Well::ProducerCMode::BHP || current == Well::ProducerCMode::THP) {
                pressure_controlled_well = true;
            }
        }
        if (pressure_controlled_well) {
            // initialized the well rates with the potentials i.e. the well rates based on bhp
            const double sign = this->well_ecl_.isInjector() ? 1.0 : -1.0;
            for (int phase = 0; phase < np; ++phase){
                well_potentials[phase] = sign * well_state.wellRates(index_of_well_)[phase];
            }
            return;
        }

        debug_cost_counter_ = 0;
        // does the well have a THP related constraint?
        const auto& summaryState = ebosSimulator.vanguard().summaryState();
        if (!Base::wellHasTHPConstraints(summaryState)) {
            computeWellRatesAtBhpLimit(ebosSimulator, well_potentials, deferred_logger);
        } else {
            well_potentials = computeWellPotentialWithTHP(ebosSimulator, deferred_logger);
        }
        deferred_logger.debug("Cost in iterations of finding well potential for well "
                              + name() + ": " + std::to_string(debug_cost_counter_));
    }




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesAtBhpLimit(const Simulator& ebosSimulator,
                               std::vector<double>& well_flux,
                               DeferredLogger& deferred_logger) const
    {
        if (well_ecl_.isInjector()) {
            const auto controls = well_ecl_.injectionControls(ebosSimulator.vanguard().summaryState());
            computeWellRatesWithBhp(ebosSimulator, controls.bhp_limit, well_flux, deferred_logger);
        } else {
            const auto controls = well_ecl_.productionControls(ebosSimulator.vanguard().summaryState());
            computeWellRatesWithBhp(ebosSimulator, controls.bhp_limit, well_flux, deferred_logger);
        }
    }



    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const Scalar bhp,
                            std::vector<double>& well_flux,
                            DeferredLogger& deferred_logger) const
    {
        // creating a copy of the well itself, to avoid messing up the explicit informations
        // during this copy, the only information not copied properly is the well controls
        MultisegmentWell<TypeTag> well_copy(*this);
        well_copy.debug_cost_counter_ = 0;

        // store a copy of the well state, we don't want to update the real well state
        WellState well_state_copy = ebosSimulator.problem().wellModel().wellState();
        const auto& group_state = ebosSimulator.problem().wellModel().groupState();

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
            well_state_copy.currentInjectionControl(index_of_well_, Well::InjectorCMode::BHP);
        } else {
            prod_controls.bhp_limit = bhp;
            well_state_copy.currentProductionControl(index_of_well_, Well::ProducerCMode::BHP);
        }
        well_state_copy.update_bhp(well_copy.index_of_well_, bhp);
        well_copy.scaleSegmentPressuresWithBhp(well_state_copy);

        // initialized the well rates with the potentials i.e. the well rates based on bhp
        const int np = number_of_phases_;
        const double sign = well_copy.well_ecl_.isInjector() ? 1.0 : -1.0;
        for (int phase = 0; phase < np; ++phase){
            well_state_copy.wellRates(well_copy.index_of_well_)[phase]
                    = sign * well_state_copy.wellPotentials()[well_copy.index_of_well_*np + phase];
        }
        well_copy.scaleSegmentRatesWithWellRates(well_state_copy);

        well_copy.calculateExplicitQuantities(ebosSimulator, well_state_copy, deferred_logger);
        const double dt = ebosSimulator.timeStepSize();
        // iterate to get a solution at the given bhp.
        well_copy.iterateWellEqWithControl(ebosSimulator, dt, inj_controls, prod_controls, well_state_copy, group_state,
                                           deferred_logger);

        // compute the potential and store in the flux vector.
        well_flux.clear();
        well_flux.resize(np, 0.0);
        for (int compIdx = 0; compIdx < num_components_; ++compIdx) {
            const EvalWell rate = well_copy.getQs(compIdx);
            well_flux[ebosCompIdxToFlowCompIdx(compIdx)] = rate.value();
        }
        debug_cost_counter_ += well_copy.debug_cost_counter_;
    }



    template<typename TypeTag>
    std::vector<double>
    MultisegmentWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebos_simulator,
                                DeferredLogger& deferred_logger) const
    {
        std::vector<double> potentials(number_of_phases_, 0.0);
        const auto& summary_state = ebos_simulator.vanguard().summaryState();

        const auto& well = well_ecl_;
        if (well.isInjector()){
            auto bhp_at_thp_limit = computeBhpAtThpLimitInj(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well_ecl_.injectionControls(summary_state);
                const double bhp = std::min(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
                deferred_logger.debug("Converged thp based potential calculation for well "
                                      + name() + ", at bhp = " + std::to_string(bhp));
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const auto& controls = well_ecl_.injectionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            }
        } else {
            auto bhp_at_thp_limit = computeBhpAtThpLimitProd(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well_ecl_.productionControls(summary_state);
                const double bhp = std::max(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
                deferred_logger.debug("Converged thp based potential calculation for well "
                                      + name() + ", at bhp = " + std::to_string(bhp));
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const auto& controls = well_ecl_.productionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            }
        }

        return potentials;
    }



    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state, DeferredLogger& /* deferred_logger */) const
    {
        // TODO: to test using rate conversion coefficients to see if it will be better than
        // this default one
        if (!this->isOperable() && !this->wellIsStopped()) return;

        const Well& well = Base::wellEcl();

        // the index of the top segment in the WellState
        const auto segment_rates = well_state.segRates(index_of_well_);
        const auto segment_pressure = well_state.segPress(index_of_well_);
        const PhaseUsage& pu = phaseUsage();

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // calculate the total rate for each segment
            double total_seg_rate = 0.0;
            // the segment pressure
            primary_variables_[seg][SPres] = segment_pressure[seg];
            // TODO: under what kind of circustances, the following will be wrong?
            // the definition of g makes the gas phase is always the last phase
            for (int p = 0; p < number_of_phases_; p++) {
                total_seg_rate += scalingFactor(p) * segment_rates[number_of_phases_ * seg + p];
            }

            primary_variables_[seg][GTotal] = total_seg_rate;
            if (std::abs(total_seg_rate) > 0.) {
                if (has_wfrac_variable) {
                    const int water_pos = pu.phase_pos[Water];
                    primary_variables_[seg][WFrac] = scalingFactor(water_pos) * segment_rates[number_of_phases_ * seg + water_pos] / total_seg_rate;
                }
                if (has_gfrac_variable) {
                    const int gas_pos = pu.phase_pos[Gas];
                    primary_variables_[seg][GFrac] = scalingFactor(gas_pos) * segment_rates[number_of_phases_ * seg + gas_pos] / total_seg_rate;
                }
            } else { // total_seg_rate == 0
                if (this->isInjector()) {
                    // only single phase injection handled
                    auto phase = well.getInjectionProperties().injectorType;

                    if (has_wfrac_variable) {
                        if (phase == InjectorType::WATER) {
                            primary_variables_[seg][WFrac] = 1.0;
                        } else {
                            primary_variables_[seg][WFrac] = 0.0;
                        }
                    }

                    if (has_gfrac_variable) {
                        if (phase == InjectorType::GAS) {
                            primary_variables_[seg][GFrac] = 1.0;
                        } else {
                            primary_variables_[seg][GFrac] = 0.0;
                        }
                    }

                } else if (this->isProducer()) { // producers
                    if (has_wfrac_variable) {
                        primary_variables_[seg][WFrac] = 1.0 / number_of_phases_;
                    }

                    if (has_gfrac_variable) {
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
        if (!this->isOperable() && !this->wellIsStopped()) return;

        BVectorWell resWell = resWell_;
        // resWell = resWell - B * x
        duneB_.mmv(x, resWell);
        // xw = D^-1 * resWell
        xw = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    solveEqAndUpdateWellState(WellState& well_state, DeferredLogger& deferred_logger)
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        const BVectorWell dx_well = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell_);

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
                    DeferredLogger& deferred_logger,
                    const double relaxation_factor) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        const double dFLimit = param_.dwell_fraction_max_;
        const double max_pressure_change = param_.max_pressure_change_ms_wells_;
        const std::vector<std::array<double, numWellEq> > old_primary_variables = primary_variables_;

        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            if (has_wfrac_variable) {
                const int sign = dwells[seg][WFrac] > 0. ? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][WFrac]) * relaxation_factor, dFLimit);
                primary_variables_[seg][WFrac] = old_primary_variables[seg][WFrac] - dx_limited;
            }

            if (has_gfrac_variable) {
                const int sign = dwells[seg][GFrac] > 0. ? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][GFrac]) * relaxation_factor, dFLimit);
                primary_variables_[seg][GFrac] = old_primary_variables[seg][GFrac] - dx_limited;
            }

            // handling the overshooting or undershooting of the fractions
            processFractions(seg);

            // update the segment pressure
            {
                const int sign = dwells[seg][SPres] > 0.? 1 : -1;
                const double dx_limited = sign * std::min(std::abs(dwells[seg][SPres]) * relaxation_factor, max_pressure_change);
                primary_variables_[seg][SPres] = std::max( old_primary_variables[seg][SPres] - dx_limited, 1e5);
            }

            // update the total rate // TODO: should we have a limitation of the total rate change?
            {
                primary_variables_[seg][GTotal] = old_primary_variables[seg][GTotal] - relaxation_factor * dwells[seg][GTotal];

                // make sure that no injector produce and no producer inject
                if (seg == 0) {
                    if (this->isInjector()) {
                        primary_variables_[seg][GTotal] = std::max( primary_variables_[seg][GTotal], 0.0);
                    } else {
                        primary_variables_[seg][GTotal] = std::min( primary_variables_[seg][GTotal], 0.0);
                    }
                }
            }

        }

        updateWellStateFromPrimaryVariables(well_state, deferred_logger);
        Base::calculateReservoirRates(well_state);
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state,
                                DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(well_state, deferred_logger);
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
               .cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0)->fluidState();
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

        auto* wellPI = &well_state.productivityIndex()[this->index_of_well_*np + 0];
        auto* connPI = &well_state.connectionProductivityIndex()[this->first_perf_*np + 0];

        setToZero(wellPI);

        const auto preferred_phase = this->well_ecl_.getPreferredPhase();
        auto subsetPerfID   = 0;

        for ( const auto& perf : *this->perf_data_){
            auto allPerfID = perf.ecl_index;

            auto connPICalc = [&wellPICalc, allPerfID](const double mobility) -> double
            {
                return wellPICalc.connectionProdIndStandard(allPerfID, mobility);
            };

            std::vector<EvalWell> mob(this->num_components_, 0.0);
            getMobility(ebosSimulator, static_cast<int>(subsetPerfID), mob);

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
    void
    MultisegmentWell<TypeTag>::
    addWellContributions(SparseMatrixAdapter& jacobian) const
    {
        const auto invDuneD = mswellhelpers::invertWithUMFPack<DiagMatWell, BVectorWell>(duneD_, duneDSolver_);

        // We need to change matrix A as follows
        // A -= C^T D^-1 B
        // D is a (nseg x nseg) block matrix with (4 x 4) blocks.
        // B and C are (nseg x ncells) block matrices with (4 x 4 blocks).
        // They have nonzeros at (i, j) only if this well has a
        // perforation at cell j connected to segment i.  The code
        // assumes that no cell is connected to more than one segment,
        // i.e. the columns of B/C have no more than one nonzero.
        for (size_t rowC = 0; rowC < duneC_.N(); ++rowC) {
            for (auto colC = duneC_[rowC].begin(), endC = duneC_[rowC].end(); colC != endC; ++colC) {
                const auto row_index = colC.index();
                for (size_t rowB = 0; rowB < duneB_.N(); ++rowB) {
                    for (auto colB = duneB_[rowB].begin(), endB = duneB_[rowB].end(); colB != endB; ++colB) {
                        const auto col_index = colB.index();
                        OffDiagMatrixBlockWellType tmp1;
                        Detail::multMatrixImpl(invDuneD[rowC][rowB], (*colB), tmp1, std::true_type());
                        typename SparseMatrixAdapter::MatrixBlock tmp2;
                        Detail::multMatrixTransposedImpl((*colC), tmp1, tmp2, std::false_type());
                        jacobian.addToBlock(row_index, col_index, tmp2);
                    }
                }
            }
        }
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
    WellSegments::CompPressureDrop
    MultisegmentWell<TypeTag>::
    compPressureDrop() const
    {
        return segmentSet().compPressureDrop();
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

        if (has_wfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
            return primary_variables_evaluation_[seg][WFrac];
        }

        if (has_gfrac_variable && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
            return primary_variables_evaluation_[seg][GFrac];
        }

        // Oil fraction
        EvalWell oil_fraction = 1.0;
        if (has_wfrac_variable) {
            oil_fraction -= primary_variables_evaluation_[seg][WFrac];
        }

        if (has_gfrac_variable) {
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
                            const double Tw,
                            const int seg,
                            const int perf,
                            const EvalWell& segment_pressure,
                            const bool& allow_cf,
                            std::vector<EvalWell>& cq_s,
                            EvalWell& perf_press,
                            double& perf_dis_gas_rate,
                            double& perf_vap_oil_rate,
                            DeferredLogger& deferred_logger) const

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
            if (!allow_cf && this->isInjector()) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell cq_p = - Tw * (mob_perfcells[comp_idx] * drawdown);
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
            if (!allow_cf && this->isProducer()) {
                return;
            }

            // for injecting perforations, we use total mobility
            EvalWell total_mob = mob_perfcells[0];
            for (int comp_idx = 1; comp_idx < num_components_; ++comp_idx) {
                total_mob += mob_perfcells[comp_idx];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - Tw * (total_mob * drawdown);

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
                    OPM_DEFLOG_THROW(NumericalIssue, "Zero d value obtained for well " << name() << " during flux calcuation"
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
        EvalWell saltConcentration;
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
            saltConcentration = extendEval(fs.saltConcentration());
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
            std::vector<EvalWell>& phase_densities = this->segment_phase_densities_[seg];

            const EvalWell seg_pressure = getSegmentPressure(seg);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                b[waterCompIdx] =
                    FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, saltConcentration);
                visc[waterCompIdx] =
                    FluidSystem::waterPvt().viscosity(pvt_region_index, temperature, seg_pressure, saltConcentration);
                // TODO: double check here
                // TODO: should not we use phaseIndex here?
                phase_densities[waterCompIdx] = b[waterCompIdx] * surf_dens[waterCompIdx];
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
                        phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx]
                                                    + rv * b[gasCompIdx] * surf_dens[oilCompIdx];
                    } else { // no oil exists
                        b[gasCompIdx] =
                            FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                        visc[gasCompIdx] =
                            FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                        phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx];
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
                        phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx]
                                                    + rs * b[oilCompIdx] * surf_dens[gasCompIdx];
                    } else { // no oil exists
                        b[oilCompIdx] =
                            FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                        visc[oilCompIdx] =
                            FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                        phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx];
                    }
                } else { // no Liquid phase
                    // it is the same with zero mix_s[Oil]
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                }
            }

            segment_phase_viscosities_[seg] = visc;

            std::vector<EvalWell> mix(mix_s);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

                const EvalWell d = 1.0 - rs * rv;

                if (rs != 0.0) { // rs > 0.0?
                    mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
                }
                if (rv != 0.0) { // rv > 0.0?
                    mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
                }
            }

            EvalWell volrat(0.0);
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                volrat += mix[comp_idx] / b[comp_idx];
            }

            segment_viscosities_[seg] = 0.;
            // calculate the average viscosity
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell fraction =  mix[comp_idx] / b[comp_idx] / volrat;
                // TODO: a little more work needs to be done to handle the negative fractions here
                segment_phase_fractions_[seg][comp_idx] = fraction; // >= 0.0 ? fraction : 0.0;
                segment_viscosities_[seg] += visc[comp_idx] * segment_phase_fractions_[seg][comp_idx];
            }

            EvalWell density(0.0);
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                density += surf_dens[comp_idx] * mix_s[comp_idx];
            }
            segment_densities_[seg] = density / volrat;

            // calculate the mass rates
            segment_mass_rates_[seg] = 0.;
            for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                const EvalWell rate = getSegmentRateUpwinding(seg, comp_idx);
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
    getBhp() const
    {
        return getSegmentPressure(0);
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
    getQs(const int comp_idx) const
    {
        return getSegmentRate(0, comp_idx);
    }





    template <typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getSegmentRateUpwinding(const int seg,
                            const size_t comp_idx) const
    {
        const int seg_upwind = upwinding_segments_[seg];
        // the result will contain the derivative with resepct to GTotal in segment seg,
        // and the derivatives with respect to WFrac GFrac in segment seg_upwind.
        // the derivative with respect to SPres should be zero.
        if (seg == 0 && this->isInjector()) {
            const Well& well = Base::wellEcl();
            auto phase = well.getInjectionProperties().injectorType;

            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                    && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == comp_idx
                    && phase == InjectorType::WATER)
                return primary_variables_evaluation_[seg][GTotal] / scalingFactor(ebosCompIdxToFlowCompIdx(comp_idx));


            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                    && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == comp_idx
                    && phase == InjectorType::OIL)
                return primary_variables_evaluation_[seg][GTotal] / scalingFactor(ebosCompIdxToFlowCompIdx(comp_idx));

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
                    && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == comp_idx
                    && phase == InjectorType::GAS)
                return primary_variables_evaluation_[seg][GTotal] / scalingFactor(ebosCompIdxToFlowCompIdx(comp_idx));

            return 0.0;
        }

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
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    getWQTotal() const
    {
        return getSegmentGTotal(0);
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
    assembleControlEq(const WellState& well_state,
                      const GroupState& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      const Well::InjectionControls& inj_controls,
                      const Well::ProductionControls& prod_controls,
                      DeferredLogger& deferred_logger)
    {

        EvalWell control_eq(0.0);

        const auto& well = well_ecl_;

        auto getRates = [&]() {
            std::vector<EvalWell> rates(3, 0.0);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                rates[Water] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                rates[Oil] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                rates[Gas] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
            }
            return rates;
        };

        if (this->wellIsStopped()) {
            control_eq = getWQTotal();
        } else if (this->isInjector() ) {
            // Find scaling factor to get injection rate,
            const InjectorType injectorType = inj_controls.injector_type;
            double scaling = 1.0;
            const auto& pu = phaseUsage();
            switch (injectorType) {
            case InjectorType::WATER:
            {
                scaling = scalingFactor(pu.phase_pos[BlackoilPhases::Aqua]);
                break;
            }
            case InjectorType::OIL:
            {
                scaling = scalingFactor(pu.phase_pos[BlackoilPhases::Liquid]);
                break;
            }
            case InjectorType::GAS:
            {
                scaling = scalingFactor(pu.phase_pos[BlackoilPhases::Vapour]);
                break;
            }
            default:
                throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
            }
            const EvalWell injection_rate = getWQTotal() / scaling;
            // Setup function for evaluation of BHP from THP (used only if needed).
            auto bhp_from_thp = [&]() {
                const auto rates = getRates();
                return calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
            };
            // Call generic implementation.
            Base::assembleControlEqInj(well_state, group_state, schedule, summaryState, inj_controls, getBhp(), injection_rate, bhp_from_thp, control_eq, deferred_logger);
        } else {
            // Find rates.
            const auto rates = getRates();
            // Setup function for evaluation of BHP from THP (used only if needed).
            auto bhp_from_thp = [&]() {
                return calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
            };
            // Call generic implementation.
            Base::assembleControlEqProd(well_state, group_state, schedule, summaryState, prod_controls, getBhp(), rates, bhp_from_thp, control_eq, deferred_logger);
        }

        // using control_eq to update the matrix and residuals
        resWell_[0][SPres] = control_eq.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[0][0][SPres][pv_idx] = control_eq.derivative(pv_idx + numEq);
        }
    }


    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateThp(WellState& well_state, DeferredLogger& deferred_logger) const
    {
        // When there is no vaild VFP table provided, we set the thp to be zero.
        if (!this->isVFPActive(deferred_logger) || this->wellIsStopped()) {
            well_state.update_thp(index_of_well_, 0.);
            return;
        }

        // the well is under other control types, we calculate the thp based on bhp and rates
        std::vector<double> rates(3, 0.0);

        const PhaseUsage& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            rates[ Water ] = well_state.wellRates(index_of_well_)[pu.phase_pos[ Water ] ];
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            rates[ Oil ] = well_state.wellRates(index_of_well_)[pu.phase_pos[ Oil ] ];
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            rates[ Gas ] = well_state.wellRates(index_of_well_)[pu.phase_pos[ Gas ] ];
        }

        const double bhp = well_state.bhp(index_of_well_);

        well_state.update_thp(index_of_well_, calculateThpFromBhp(rates, bhp, deferred_logger));

    }



    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    calculateThpFromBhp(const std::vector<double>& rates,
                        const double bhp,
                        DeferredLogger& deferred_logger) const
    {
        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        const double aqua = rates[Water];
        const double liquid = rates[Oil];
        const double vapour = rates[Gas];

        // pick the density in the top segment
        const double rho = getRefDensity();

        double thp = 0.0;
        if (this->isInjector()) {
            const int table_id = well_ecl_.vfp_table_number();
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id).getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
        }
        else if (this->isProducer()) {
            const int table_id = well_ecl_.vfp_table_number();
            const double alq = well_ecl_.alq_value();
            const double vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id).getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
        }
        else {
            OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
        }

        return thp;
    }


    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    getRefDensity() const
    {
        return segment_densities_[0].value();
    }

    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleDefaultPressureEq(const int seg, WellState& well_state) const
    {
        assert(seg != 0); // not top segment

        // for top segment, the well control equation will be used.
        EvalWell pressure_equation = getSegmentPressure(seg);

        // we need to handle the pressure difference between the two segments
        // we only consider the hydrostatic pressure loss first
        // TODO: we might be able to add member variables to store these values, then we update well state
        // after converged
        const auto hydro_pressure_drop = getHydroPressureLoss(seg);
        auto& segments = well_state.segments(this->index_of_well_);
        segments.pressure_drop_hydrostatic[seg] = hydro_pressure_drop.value();
        pressure_equation -= hydro_pressure_drop;

        if (frictionalPressureLossConsidered()) {
            const auto friction_pressure_drop = getFrictionPressureLoss(seg);
            pressure_equation -= friction_pressure_drop;
            segments.pressure_drop_friction[seg] = friction_pressure_drop.value();
        }

        resWell_[seg][SPres] = pressure_equation.value();
        const int seg_upwind = upwinding_segments_[seg];
        duneD_[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + numEq);
        duneD_[seg][seg][SPres][GTotal] += pressure_equation.derivative(GTotal + numEq);
        if (has_wfrac_variable) {
            duneD_[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + numEq);
        }
        if (has_gfrac_variable) {
            duneD_[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + numEq);
        }

        // contribution from the outlet segment
        const int outlet_segment_index = segmentNumberToIndex(segmentSet()[seg].outletSegment());
        const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

        resWell_[seg][SPres] -= outlet_pressure.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + numEq);
        }

        if (accelerationalPressureLossConsidered()) {
            handleAccelerationPressureLoss(seg, well_state);
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
        const int seg_upwind = upwinding_segments_[seg];
        EvalWell density = segment_densities_[seg_upwind];
        EvalWell visc = segment_viscosities_[seg_upwind];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (seg != seg_upwind) {
            density.clearDerivatives();
            visc.clearDerivatives();
        }
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
    handleAccelerationPressureLoss(const int seg, WellState& well_state) const
    {
        const double area = segmentSet()[seg].crossArea();
        const EvalWell mass_rate = segment_mass_rates_[seg];
        const int seg_upwind = upwinding_segments_[seg];
        EvalWell density = segment_densities_[seg_upwind];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (seg != seg_upwind) {
            density.clearDerivatives();
        }

        EvalWell accelerationPressureLoss = mswellhelpers::velocityHead(area, mass_rate, density);
        // handling the velocity head of intlet segments
        for (const int inlet : segment_inlets_[seg]) {
            const int seg_upwind_inlet = upwinding_segments_[inlet];
            const double inlet_area = segmentSet()[inlet].crossArea();
            EvalWell inlet_density = segment_densities_[seg_upwind_inlet];
            // WARNING
            // We disregard the derivatives from the upwind density to make sure derivatives
            // wrt. to different segments dont get mixed.
            if (inlet != seg_upwind_inlet) {
                inlet_density.clearDerivatives();
            }
            const EvalWell inlet_mass_rate = segment_mass_rates_[inlet];
            accelerationPressureLoss -= mswellhelpers::velocityHead(std::max(inlet_area, area), inlet_mass_rate, inlet_density);
        }

        // We change the sign of the accelerationPressureLoss for injectors.
        // Is this correct? Testing indicates that this is what the reference simulator does
        const double sign = mass_rate < 0. ? 1.0 : - 1.0;
        accelerationPressureLoss *= sign;

        well_state.segPressDropAcceleration(index_of_well_)[seg] = accelerationPressureLoss.value();

        resWell_[seg][SPres] -= accelerationPressureLoss.value();
        duneD_[seg][seg][SPres][SPres] -= accelerationPressureLoss.derivative(SPres + numEq);
        duneD_[seg][seg][SPres][GTotal] -= accelerationPressureLoss.derivative(GTotal + numEq);
        if (has_wfrac_variable) {
            duneD_[seg][seg_upwind][SPres][WFrac] -= accelerationPressureLoss.derivative(WFrac + numEq);
        }
        if (has_gfrac_variable) {
            duneD_[seg][seg_upwind][SPres][GFrac] -= accelerationPressureLoss.derivative(GFrac + numEq);
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


    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkOperabilityUnderBHPLimitProducer(const WellState& /*well_state*/, const Simulator& ebos_simulator, DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const double bhp_limit = Base::mostStrictBhpFromBhpLimits(summaryState);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints(summaryState) ) {
            // if the BHP limit is not defaulted or the well does not have a THP limit
            // we need to check the BHP limit

            double temp = 0;
            for (int p = 0; p < number_of_phases_; ++p) {
                temp += ipr_a_[p] - ipr_b_[p] * bhp_limit;
            }
            if (temp < 0.) {
                this->operability_status_.operable_under_only_bhp_limit = false;
            }

            // checking whether running under BHP limit will violate THP limit
            if (this->operability_status_.operable_under_only_bhp_limit && this->wellHasTHPConstraints(summaryState)) {
                // option 1: calculate well rates based on the BHP limit.
                // option 2: stick with the above IPR curve
                // we use IPR here
                std::vector<double> well_rates_bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                const double thp = calculateThpFromBhp(well_rates_bhp_limit, bhp_limit, deferred_logger);

                const double thp_limit = this->getTHPConstraint(summaryState);

                if (thp < thp_limit) {
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

        // TODO: it only handles the producers for now
        // the formular for the injectors are not formulated yet
        if (this->isInjector()) {
            return;
        }

        // initialize all the values to be zero to begin with
        std::fill(ipr_a_.begin(), ipr_a_.end(), 0.);
        std::fill(ipr_b_.begin(), ipr_b_.end(), 0.);

        const int nseg = numberOfSegments();
        double seg_bhp_press_diff = 0;
        double ref_depth = ref_depth_;
        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the perforation rate for each perforation that belongs to this segment
            const double segment_depth = segmentSet()[seg].depth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth, segment_depth, segment_densities_[seg].value(), gravity_);
            ref_depth = segment_depth;
            seg_bhp_press_diff += dp;
            for (const int perf : segment_perforations_[seg]) {
            //std::vector<EvalWell> mob(num_components_, {numWellEq_ + numEq, 0.0});
            std::vector<EvalWell> mob(num_components_, 0.0);

            // TODO: mabye we should store the mobility somewhere, so that we only need to calculate it one per iteration
            getMobility(ebos_simulator, perf, mob);

            const int cell_idx = well_cells_[perf];
            const auto& int_quantities = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = int_quantities.fluidState();
            // the pressure of the reservoir grid block the well connection is in
                    // pressure difference between the segment and the perforation
            const double perf_seg_press_diff = gravity_ * segment_densities_[seg].value() * perforation_segment_depth_diffs_[perf];
            // pressure difference between the perforation and the grid cell
            const double cell_perf_press_diff = cell_perforation_pressure_diffs_[perf];
            const double pressure_cell = fs.pressure(FluidSystem::oilPhaseIdx).value();

            // calculating the b for the connection
            std::vector<double> b_perf(num_components_);
            for (size_t phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }
                const unsigned comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phase));
                b_perf[comp_idx] = fs.invB(phase).value();
            }

            // the pressure difference between the connection and BHP
            const double h_perf = cell_perf_press_diff + perf_seg_press_diff + seg_bhp_press_diff;
            const double pressure_diff = pressure_cell - h_perf;

            // Let us add a check, since the pressure is calculated based on zero value BHP
            // it should not be negative anyway. If it is negative, we might need to re-formulate
            // to taking into consideration the crossflow here.
            if (pressure_diff <= 0.) {
                deferred_logger.warning("NON_POSITIVE_DRAWDOWN_IPR",
                                "non-positive drawdown found when updateIPR for well " + name());
            }

            // the well index associated with the connection
            const double tw_perf = well_index_[perf]*ebos_simulator.problem().template rockCompTransMultiplier<double>(int_quantities, cell_idx);

            // TODO: there might be some indices related problems here
            // phases vs components
            // ipr values for the perforation
            std::vector<double> ipr_a_perf(ipr_a_.size());
            std::vector<double> ipr_b_perf(ipr_b_.size());
            for (int p = 0; p < number_of_phases_; ++p) {
                const double tw_mob = tw_perf * mob[p].value() * b_perf[p];
                ipr_a_perf[p] += tw_mob * pressure_diff;
                ipr_b_perf[p] += tw_mob;
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

            for (int p = 0; p < number_of_phases_; ++p) {
                // TODO: double check the indices here
                ipr_a_[ebosCompIdxToFlowCompIdx(p)] += ipr_a_perf[p];
                ipr_b_[ebosCompIdxToFlowCompIdx(p)] += ipr_b_perf[p];
            }
            }
        }
    }

    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkOperabilityUnderTHPLimitProducer(const Simulator& ebos_simulator, const WellState& /*well_state*/, DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto obtain_bhp = computeBhpAtThpLimitProd(ebos_simulator, summaryState, deferred_logger);

        if (obtain_bhp) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const double  bhp_limit = Base::mostStrictBhpFromBhpLimits(summaryState);
            this->operability_status_.obey_bhp_limit_with_thp_limit = (*obtain_bhp >= bhp_limit);

            const double thp_limit = this->getTHPConstraint(summaryState);
            if (*obtain_bhp < thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is SMALLER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a producer for well " + name();
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
                                      + " bar for well " + name() + ", the well might need to be closed ");
            }
        }
    }





    template <typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellStateFromPrimaryVariables(WellState& well_state, DeferredLogger& deferred_logger) const
    {
        const PhaseUsage& pu = phaseUsage();
        assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
        const int oil_pos = pu.phase_pos[Oil];

        auto segment_rates = well_state.segRates(this->index_of_well_);
        auto segment_pressure = well_state.segPress(this->index_of_well_);
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
            for (int p = 0; p < number_of_phases_; ++p) {
                const double phase_rate = g_total * fractions[p];
                segment_rates[seg*this->number_of_phases_ + p] = phase_rate;
                if (seg == 0) { // top segment
                    well_state.wellRates(index_of_well_)[p] = phase_rate;
                }
            }

            // update the segment pressure
            segment_pressure[seg] = primary_variables_[seg][SPres];
            if (seg == 0) { // top segment
                well_state.update_bhp(index_of_well_, segment_pressure[seg]);
            }
        }
        updateThp(well_state, deferred_logger);
    }




    template <typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    frictionalPressureLossConsidered() const
    {
        // HF- and HFA needs to consider frictional pressure loss
        return (segmentSet().compPressureDrop() != WellSegments::CompPressureDrop::H__);
    }





    template <typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    accelerationalPressureLossConsidered() const
    {
        return (segmentSet().compPressureDrop() == WellSegments::CompPressureDrop::HFA);
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
        if (!this->isOperable() && !this->wellIsStopped()) return true;

        const int max_iter_number = param_.max_inner_iter_ms_wells_;
        const WellState well_state0 = well_state;
        const std::vector<Scalar> residuals0 = getWellResiduals(Base::B_avg_, deferred_logger);
        std::vector<std::vector<Scalar> > residual_history;
        std::vector<double> measure_history;
        int it = 0;
        // relaxation factor
        double relaxation_factor = 1.;
        const double min_relaxation_factor = 0.6;
        bool converged = false;
        int stagnate_count = 0;
        bool relax_convergence = false;
        for (; it < max_iter_number; ++it, ++debug_cost_counter_) {

            assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            const BVectorWell dx_well = mswellhelpers::applyUMFPack(duneD_, duneDSolver_, resWell_);

            if (it > param_.strict_inner_iter_ms_wells_)
                relax_convergence = true;

            const auto report = getWellConvergence(well_state, Base::B_avg_, deferred_logger, relax_convergence);
            if (report.converged()) {
                converged = true;
                break;
            }

            residual_history.push_back(getWellResiduals(Base::B_avg_, deferred_logger));
            measure_history.push_back(getResidualMeasureValue(well_state, residual_history[it], deferred_logger) );

            bool is_oscillate = false;
            bool is_stagnate = false;

            detectOscillations(measure_history, it, is_oscillate, is_stagnate);
            // TODO: maybe we should have more sophiscated strategy to recover the relaxation factor,
            // for example, to recover it to be bigger

            if (is_oscillate || is_stagnate) {
                // HACK!
                std::ostringstream sstr;
                if (relaxation_factor == min_relaxation_factor) {
                    // Still stagnating, terminate iterations if 5 iterations pass.
                    ++stagnate_count;
                    if (stagnate_count == 6) {
                        sstr << " well " << name() << " observes severe stagnation and/or oscillation. We relax the tolerance and check for convergence. \n";
                        const auto reportStag = getWellConvergence(well_state, Base::B_avg_, deferred_logger, true);
                        if (reportStag.converged()) {
                            converged = true;
                            sstr << " well " << name() << " manages to get converged with relaxed tolerances in " << it << " inner iterations";
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
                    sstr << " well " << name() << " observes stagnation in inner iteration " << it << "\n";

                }
                if (is_oscillate) {
                    sstr << " well " << name() << " observes oscillation in inner iteration " << it << "\n";
                }
                sstr << " relaxation_factor is " << relaxation_factor << " now\n";
                deferred_logger.debug(sstr.str());
            }
            updateWellState(dx_well, well_state, deferred_logger, relaxation_factor);
            initPrimaryVariablesEvaluation();
        }

        // TODO: we should decide whether to keep the updated well_state, or recover to use the old well_state
        if (converged) {
            std::ostringstream sstr;
            sstr << "     Well " << name() << " converged in " << it << " inner iterations.";
            if (relax_convergence)
                sstr << "      (A relaxed tolerance was used after "<< param_.strict_inner_iter_ms_wells_ << " iterations)";
            deferred_logger.debug(sstr.str());
        } else {
            std::ostringstream sstr;
            sstr << "     Well " << name() << " did not converge in " << it << " inner iterations.";
#define EXTRA_DEBUG_MSW 0
#if EXTRA_DEBUG_MSW
            sstr << "***** Outputting the residual history for well " << name() << " during inner iterations:";
            for (int i = 0; i < it; ++i) {
                const auto& residual = residual_history[i];
                sstr << " residual at " << i << "th iteration ";
                for (const auto& res : residual) {
                    sstr << " " << res;
                }
                sstr << " " << measure_history[i] << " \n";
            }
#endif
            deferred_logger.debug(sstr.str());
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

        if (!this->isOperable() && !this->wellIsStopped()) return;

        // update the upwinding segments
        updateUpwindingSegments();

        // calculate the fluid properties needed.
        computeSegmentFluidProperties(ebosSimulator);

        // clear all entries
        duneB_ = 0.0;
        duneC_ = 0.0;

        duneD_ = 0.0;
        resWell_ = 0.0;

        duneDSolver_.reset();

        well_state.wellVaporizedOilRates(index_of_well_) = 0.;
        well_state.wellDissolvedGasRates(index_of_well_) = 0.;

        // for the black oil cases, there will be four equations,
        // the first three of them are the mass balance equations, the last one is the pressure equations.
        //
        // but for the top segment, the pressure equation will be the well control equation, and the other three will be the same.

        const bool allow_cf = getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);

        const int nseg = numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the accumulation term
            // TODO: without considering the efficiencty factor for now
            {
                const EvalWell segment_surface_volume = getSegmentSurfaceVolume(ebosSimulator, seg);

                // Add a regularization_factor to increase the accumulation term
                // This will make the system less stiff and help convergence for
                // difficult cases
                const Scalar regularization_factor =  param_.regularization_factor_ms_wells_;
                // for each component
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    const EvalWell accumulation_term = regularization_factor * (segment_surface_volume * surfaceVolumeFraction(seg, comp_idx)
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
                    const EvalWell segment_rate = getSegmentRateUpwinding(seg, comp_idx) * well_efficiency_factor_;

                    const int seg_upwind = upwinding_segments_[seg];
                    // segment_rate contains the derivatives with respect to GTotal in seg,
                    // and WFrac and GFrac in seg_upwind
                    resWell_[seg][comp_idx] -= segment_rate.value();
                    duneD_[seg][seg][comp_idx][GTotal] -= segment_rate.derivative(GTotal + numEq);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        duneD_[seg][seg_upwind][comp_idx][WFrac] -= segment_rate.derivative(WFrac + numEq);
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        duneD_[seg][seg_upwind][comp_idx][GFrac] -= segment_rate.derivative(GFrac + numEq);
                    }
                    // pressure derivative should be zero
                }
            }

            // considering the contributions from the inlet segments
            {
                for (const int inlet : segment_inlets_[seg]) {
                    for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                        const EvalWell inlet_rate = getSegmentRateUpwinding(inlet, comp_idx) * well_efficiency_factor_;

                        const int inlet_upwind = upwinding_segments_[inlet];
                        // inlet_rate contains the derivatives with respect to GTotal in inlet,
                        // and WFrac and GFrac in inlet_upwind
                        resWell_[seg][comp_idx] += inlet_rate.value();
                        duneD_[seg][inlet][comp_idx][GTotal] += inlet_rate.derivative(GTotal + numEq);
                        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                            duneD_[seg][inlet_upwind][comp_idx][WFrac] += inlet_rate.derivative(WFrac + numEq);
                        }
                        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                            duneD_[seg][inlet_upwind][comp_idx][GFrac] += inlet_rate.derivative(GFrac + numEq);
                        }
                        // pressure derivative should be zero
                    }
                }
            }

            // calculating the perforation rate for each perforation that belongs to this segment
            const EvalWell seg_pressure = getSegmentPressure(seg);
            const int rate_start_offset = first_perf_ * number_of_phases_;
            auto * perf_rates = &well_state.mutable_perfPhaseRates()[rate_start_offset];
            auto& perf_press_state = well_state.perfPress(this->index_of_well_);
            for (const int perf : segment_perforations_[seg]) {
                const int cell_idx = well_cells_[perf];
                const auto& int_quants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                std::vector<EvalWell> mob(num_components_, 0.0);
                getMobility(ebosSimulator, perf, mob);
                const double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(int_quants, cell_idx);
                const double Tw = well_index_[perf] * trans_mult;
                std::vector<EvalWell> cq_s(num_components_, 0.0);
                EvalWell perf_press;
                double perf_dis_gas_rate = 0.;
                double perf_vap_oil_rate = 0.;
                computePerfRatePressure(int_quants, mob, Tw, seg, perf, seg_pressure, allow_cf, cq_s, perf_press, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);

                // updating the solution gas rate and solution oil rate
                if (this->isProducer()) {
                    well_state.wellDissolvedGasRates(index_of_well_) += perf_dis_gas_rate;
                    well_state.wellVaporizedOilRates(index_of_well_) += perf_vap_oil_rate;
                }

                // store the perf pressure and rates
                for (int comp_idx = 0; comp_idx < num_components_; ++comp_idx) {
                    perf_rates[perf*number_of_phases_ + ebosCompIdxToFlowCompIdx(comp_idx)] = cq_s[comp_idx].value();
                }
                perf_press_state[perf] = perf_press.value();

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
                const auto& summaryState = ebosSimulator.vanguard().summaryState();
                const Schedule& schedule = ebosSimulator.vanguard().schedule();
                assembleControlEq(well_state, group_state, schedule, summaryState, inj_controls, prod_controls, deferred_logger);
            } else {
                const UnitSystem& unit_system = ebosSimulator.vanguard().eclState().getDeckUnitSystem();
                assemblePressureEq(seg, unit_system, well_state, deferred_logger);
            }
        }
    }




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assemblePressureEq(const int seg, const UnitSystem& unit_system,
                       WellState& well_state, DeferredLogger& deferred_logger) const
    {
        switch(segmentSet()[seg].segmentType()) {
            case Segment::SegmentType::SICD :
            case Segment::SegmentType::AICD :
            case Segment::SegmentType::VALVE : {
                assembleICDPressureEq(seg, unit_system, well_state,deferred_logger);
                break;
            }
            default :
                assembleDefaultPressureEq(seg, well_state);
        }
    }



    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const
    {
        return !getAllowCrossFlow() && allDrawDownWrongDirection(ebos_simulator);
    }


    template<typename TypeTag>
    bool
    MultisegmentWell<TypeTag>::
    allDrawDownWrongDirection(const Simulator& ebos_simulator) const
    {
        bool all_drawdown_wrong_direction = true;
        const int nseg = numberOfSegments();

        for (int seg = 0; seg < nseg; ++seg) {
            const EvalWell segment_pressure = getSegmentPressure(seg);
            for (const int perf : segment_perforations_[seg]) {

                const int cell_idx = well_cells_[perf];
                const auto& intQuants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                const auto& fs = intQuants.fluidState();

                // pressure difference between the segment and the perforation
                const EvalWell perf_seg_press_diff = gravity_ * segment_densities_[seg] * perforation_segment_depth_diffs_[perf];
                // pressure difference between the perforation and the grid cell
                const double cell_perf_press_diff = cell_perforation_pressure_diffs_[perf];

                const double pressure_cell = (fs.pressure(FluidSystem::oilPhaseIdx)).value();
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
    updateWaterThroughput(const double dt OPM_UNUSED, WellState& well_state OPM_UNUSED) const
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
            const int cell_idx = well_cells_[0];
            const auto& intQuants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            temperature.setValue(fs.temperature(FluidSystem::oilPhaseIdx).value());
            saltConcentration = extendEval(fs.saltConcentration());
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
                FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, saltConcentration);
        }

        EvalWell rv(0.0);
        // gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index, temperature, seg_pressure);
                if (rvmax < 0.0) { // negative rvmax can happen if the seg_pressure is outside the range of the table
                    rvmax = 0.0;
                }
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
                EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index, temperature, seg_pressure);
                if (rsmax < 0.0) { // negative rsmax can happen if the seg_pressure is outside the range of the table
                    rsmax = 0.0;
                }
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
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << name()
                     << " during conversion to surface volume with rs " << rs
                     << ", rv " << rv << " and pressure " << seg_pressure
                     << " obtaining d " << d;
                OpmLog::debug(sstr.str());
                OPM_THROW_NOLOG(NumericalIssue, sstr.str());
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

        // We increase the segment volume with a factor 10 to stabilize the system.
        const double volume = segmentSet()[seg_idx].volume();

        return volume / vol_ratio;
    }





    template<typename TypeTag>
    std::vector<typename MultisegmentWell<TypeTag>::Scalar>
    MultisegmentWell<TypeTag>::
    getWellResiduals(const std::vector<Scalar>& B_avg,
                     DeferredLogger& deferred_logger) const
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
                    OPM_DEFLOG_THROW(NumericalIssue, "nan or inf value for residal get for well " << name()
                                                    << " segment " << seg << " eq_idx " << eq_idx, deferred_logger);
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
               OPM_DEFLOG_THROW(NumericalIssue, "nan or inf value for control residal get for well " << name(), deferred_logger);
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
    getResidualMeasureValue(const WellState& well_state,
                            const std::vector<double>& residuals,
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

        const double control_tolerance = getControlTolerance(well_state, deferred_logger);
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
    getControlTolerance(const WellState& well_state,
                        DeferredLogger& deferred_logger) const
    {
        double control_tolerance = 0.;

        const int well_index = index_of_well_;
        if (this->isInjector() )
        {
            auto current = well_state.currentInjectionControl(well_index);
            switch(current) {
            case Well::InjectorCMode::THP:
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case Well::InjectorCMode::BHP:
                control_tolerance = param_.tolerance_wells_;
                break;
            case Well::InjectorCMode::RATE:
            case Well::InjectorCMode::RESV:
                control_tolerance = param_.tolerance_wells_;
                break;
            case Well::InjectorCMode::GRUP:
                control_tolerance = param_.tolerance_wells_;
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }

        if (this->isProducer() )
        {
            auto current = well_state.currentProductionControl(well_index);
            switch(current) {
            case Well::ProducerCMode::THP:
                control_tolerance = param_.tolerance_pressure_ms_wells_; // 0.1 bar
                break;
            case Well::ProducerCMode::BHP:
                control_tolerance = param_.tolerance_wells_; // 0.01 bar
                break;
            case Well::ProducerCMode::ORAT:
            case Well::ProducerCMode::WRAT:
            case Well::ProducerCMode::GRAT:
            case Well::ProducerCMode::LRAT:
            case Well::ProducerCMode::RESV:
            case Well::ProducerCMode::CRAT:
                control_tolerance = param_.tolerance_wells_; // smaller tolerance for rate control
                break;
            case Well::ProducerCMode::GRUP:
                control_tolerance = param_.tolerance_wells_; // smaller tolerance for rate control
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }

        return control_tolerance;
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    checkConvergenceControlEq(const WellState& well_state,
                              ConvergenceReport& report,
                              DeferredLogger& deferred_logger) const
    {
        double control_tolerance = 0.;
        using CR = ConvergenceReport;
        CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

        const int well_index = index_of_well_;
        if (this->isInjector() )
        {
            auto current = well_state.currentInjectionControl(well_index);
            switch(current) {
            case Well::InjectorCMode::THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case Well::InjectorCMode::BHP:
                ctrltype = CR::WellFailure::Type::ControlBHP;
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case Well::InjectorCMode::RATE:
            case Well::InjectorCMode::RESV:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = param_.tolerance_wells_;
                break;
            case Well::InjectorCMode::GRUP:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = param_.tolerance_wells_;
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }

        if (this->isProducer() )
        {
            auto current = well_state.currentProductionControl(well_index);
            switch(current) {
            case Well::ProducerCMode::THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case Well::ProducerCMode::BHP:
                ctrltype = CR::WellFailure::Type::ControlBHP;
                control_tolerance = param_.tolerance_pressure_ms_wells_;
                break;
            case Well::ProducerCMode::ORAT:
            case Well::ProducerCMode::WRAT:
            case Well::ProducerCMode::GRAT:
            case Well::ProducerCMode::LRAT:
            case Well::ProducerCMode::RESV:
            case Well::ProducerCMode::CRAT:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = param_.tolerance_wells_;
                break;
            case Well::ProducerCMode::GRUP:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = param_.tolerance_wells_;
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }

        const double well_control_residual = std::abs(resWell_[0][SPres]);
        const int dummy_component = -1;
        const double max_residual_allowed = param_.max_residual_allowed_;
        if (std::isnan(well_control_residual)) {
            report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, name()});
        } else if (well_control_residual > max_residual_allowed * 10.) {
            report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, name()});
        } else if ( well_control_residual > control_tolerance) {
            report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, name()});
        }
    }






    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateUpwindingSegments()
    {
        for (int seg = 0; seg < numberOfSegments(); ++seg) {
            // special treatment is needed for segment 0
            if (seg == 0) {
                // we are not supposed to have injecting producers and producing injectors
                assert( ! (this->isProducer() && primary_variables_evaluation_[seg][GTotal] > 0.) );
                assert( ! (this->isInjector() && primary_variables_evaluation_[seg][GTotal] < 0.) );
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




    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    assembleICDPressureEq(const int seg, const UnitSystem& unit_system,
                          WellState& well_state, DeferredLogger& deferred_logger) const
    {
        // TODO: upwinding needs to be taken care of
        // top segment can not be a spiral ICD device
        assert(seg != 0);

        // the pressure equation is something like
        // p_seg - deltaP - p_outlet = 0.
        // the major part is how to calculate the deltaP

        EvalWell pressure_equation = getSegmentPressure(seg);

        EvalWell icd_pressure_drop;
        switch(segmentSet()[seg].segmentType()) {
            case Segment::SegmentType::SICD :
                icd_pressure_drop = pressureDropSpiralICD(seg);
                break;
            case Segment::SegmentType::AICD :
                icd_pressure_drop = pressureDropAutoICD(seg, unit_system);
                break;
            case Segment::SegmentType::VALVE :
                icd_pressure_drop = pressureDropValve(seg);
                break;
            default: {
                OPM_DEFLOG_THROW(std::runtime_error, "Segment " + std::to_string(segmentSet()[seg].segmentNumber())
                                 + " for well " + name() + " is not of ICD type", deferred_logger);
            }
        }
        pressure_equation = pressure_equation - icd_pressure_drop;
        well_state.segments(this->index_of_well_).pressure_drop_friction[seg] = icd_pressure_drop.value();

        const int seg_upwind = upwinding_segments_[seg];
        resWell_[seg][SPres] = pressure_equation.value();
        duneD_[seg][seg][SPres][SPres] += pressure_equation.derivative(SPres + numEq);
        duneD_[seg][seg][SPres][GTotal] += pressure_equation.derivative(GTotal + numEq);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            duneD_[seg][seg_upwind][SPres][WFrac] += pressure_equation.derivative(WFrac + numEq);
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            duneD_[seg][seg_upwind][SPres][GFrac] += pressure_equation.derivative(GFrac + numEq);
        }

        // contribution from the outlet segment
        const int outlet_segment_index = segmentNumberToIndex(segmentSet()[seg].outletSegment());
        const EvalWell outlet_pressure = getSegmentPressure(outlet_segment_index);

        resWell_[seg][SPres] -= outlet_pressure.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            duneD_[seg][outlet_segment_index][SPres][pv_idx] = -outlet_pressure.derivative(pv_idx + numEq);

        }
    }




    template<typename TypeTag>
    std::optional<double>
    MultisegmentWell<TypeTag>::
    computeBhpAtThpLimitProd(const Simulator& ebos_simulator,
                             const SummaryState& summary_state,
                             DeferredLogger& deferred_logger) const
    {
        // Given a VFP function returning bhp as a function of phase
        // rates and thp:
        //     fbhp(rates, thp),
        // a function extracting the particular flow rate used for VFP
        // lookups:
        //     flo(rates)
        // and the inflow function (assuming the reservoir is fixed):
        //     frates(bhp)
        // we want to solve the equation:
        //     fbhp(frates(bhp, thplimit)) - bhp = 0
        // for bhp.
        //
        // This may result in 0, 1 or 2 solutions. If two solutions,
        // the one corresponding to the lowest bhp (and therefore
        // highest rate) should be returned.

        // Make the fbhp() function.
        const auto& controls = well_ecl_.productionControls(summary_state);
        const auto& table = vfp_properties_->getProd()->getTable(controls.vfp_table_number);
        const double vfp_ref_depth = table.getDatumDepth();
        const double rho = getRefDensity(); // Use the density at the top perforation.
        const double thp_limit = this->getTHPConstraint(summary_state);
        const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
        auto fbhp = [this, &controls, thp_limit, dp](const std::vector<double>& rates) {
            assert(rates.size() == 3);
            return this->vfp_properties_->getProd()
            ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit, controls.alq_value) - dp;
        };

        // Make the flo() function.
        auto flo = [&table](const std::vector<double>& rates) {
            return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
        };

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

        // Find the bhp-point where production becomes nonzero.
        double bhp_max = 0.0;
        {
            auto fflo = [&flo, &frates](double bhp) { return flo(frates(bhp)); };
            double low = controls.bhp_limit;
            double high = maxPerfPress(ebos_simulator) + 1.0 * unit::barsa;
            double f_low = fflo(low);
            double f_high = fflo(high);
            deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + name() +
                                  "  low = " + std::to_string(low) +
                                  "  high = " + std::to_string(high) +
                                  "  f(low) = " + std::to_string(f_low) +
                                  "  f(high) = " + std::to_string(f_high));
            int adjustments = 0;
            const int max_adjustments = 10;
            const double adjust_amount = 5.0 * unit::barsa;
            while (f_low * f_high > 0.0 && adjustments < max_adjustments) {
                // Same sign, adjust high to see if we can flip it.
                high += adjust_amount;
                f_high = fflo(high);
                ++adjustments;
            }
            if (f_low * f_high > 0.0) {
                if (f_low > 0.0) {
                    // Even at the BHP limit, we are injecting.
                    // There will be no solution here, return an
                    // empty optional.
                    deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                            "Robust bhp(thp) solve failed due to inoperability for well " + name());
                    return std::optional<double>();
                } else {
                    // Still producing, even at high bhp.
                    assert(f_high < 0.0);
                    bhp_max = high;
                }
            } else {
                // Bisect to find a bhp point where we produce, but
                // not a large amount ('eps' below).
                const double eps = 0.1 * std::fabs(table.getFloAxis().front());
                const int maxit = 50;
                int it = 0;
                while (std::fabs(f_low) > eps && it < maxit) {
                    const double curr = 0.5*(low + high);
                    const double f_curr = fflo(curr);
                    if (f_curr * f_low > 0.0) {
                        low = curr;
                        f_low = f_curr;
                    } else {
                        high = curr;
                        f_high = f_curr;
                    }
                    ++it;
                }
                bhp_max = low;
            }
            deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + name() +
                                  "  low = " + std::to_string(low) +
                                  "  high = " + std::to_string(high) +
                                  "  f(low) = " + std::to_string(f_low) +
                                  "  f(high) = " + std::to_string(f_high) +
                                  "  bhp_max = " + std::to_string(bhp_max));
        }

        // Define the equation we want to solve.
        auto eq = [&fbhp, &frates](double bhp) {
            return fbhp(frates(bhp)) - bhp;
        };

        // Find appropriate brackets for the solution.
        double low = controls.bhp_limit;
        double high = bhp_max;
        {
            double eq_high = eq(high);
            double eq_low = eq(low);
            const double eq_bhplimit = eq_low;
            deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + name() +
                                  "  low = " + std::to_string(low) +
                                  "  high = " + std::to_string(high) +
                                  "  eq(low) = " + std::to_string(eq_low) +
                                  "  eq(high) = " + std::to_string(eq_high));
            if (eq_low * eq_high > 0.0) {
                // Failed to bracket the zero.
                // If this is due to having two solutions, bisect until bracketed.
                double abs_low = std::fabs(eq_low);
                double abs_high = std::fabs(eq_high);
                int bracket_attempts = 0;
                const int max_bracket_attempts = 20;
                double interval = high - low;
                const double min_interval = 1.0 * unit::barsa;
                while (eq_low * eq_high > 0.0 && bracket_attempts < max_bracket_attempts && interval > min_interval) {
                    if (abs_high < abs_low) {
                        low = 0.5 * (low + high);
                        eq_low = eq(low);
                        abs_low = std::fabs(eq_low);
                    } else {
                        high = 0.5 * (low + high);
                        eq_high = eq(high);
                        abs_high = std::fabs(eq_high);
                    }
                    ++bracket_attempts;
                }
                if (eq_low * eq_high > 0.0) {
                    // Still failed bracketing!
                    const double limit = 3.0 * unit::barsa;
                    if (std::min(abs_low, abs_high) < limit) {
                        // Return the least bad solution if less off than 3 bar.
                        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                                "Robust bhp(thp) not solved precisely for well " + name());
                        return abs_low < abs_high ? low : high;
                    } else {
                        // Return failure.
                        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                                "Robust bhp(thp) solve failed due to bracketing failure for well " + name());
                        return std::optional<double>();
                    }
                }
            }
            // We have a bracket!
            // Now, see if (bhplimit, low) is a bracket in addition to (low, high).
            // If so, that is the bracket we shall use, choosing the solution with the
            // highest flow.
            if (eq_low * eq_bhplimit <= 0.0) {
                high = low;
                low = controls.bhp_limit;
            }
        }

        // Solve for the proper solution in the given interval.
        const int max_iteration = 100;
        const double bhp_tolerance = 0.01 * unit::barsa;
        int iteration = 0;
        try {
            const double solved_bhp = RegulaFalsiBisection<ThrowOnError>::
                solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
            return solved_bhp;
        }
        catch (...) {
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
	}
    }




    template<typename TypeTag>
    std::optional<double>
    MultisegmentWell<TypeTag>::
    computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                            const SummaryState& summary_state,
                            DeferredLogger& deferred_logger) const
    {
        // Given a VFP function returning bhp as a function of phase
        // rates and thp:
        //     fbhp(rates, thp),
        // a function extracting the particular flow rate used for VFP
        // lookups:
        //     flo(rates)
        // and the inflow function (assuming the reservoir is fixed):
        //     frates(bhp)
        // we want to solve the equation:
        //     fbhp(frates(bhp, thplimit)) - bhp = 0
        // for bhp.
        //
        // This may result in 0, 1 or 2 solutions. If two solutions,
        // the one corresponding to the lowest bhp (and therefore
        // highest rate) is returned.
        //
        // In order to detect these situations, we will find piecewise
        // linear approximations both to the inverse of the frates
        // function and to the fbhp function.
        //
        // We first take the FLO sample points of the VFP curve, and
        // find the corresponding bhp values by solving the equation:
        //     flo(frates(bhp)) - flo_sample = 0
        // for bhp, for each flo_sample. The resulting (flo_sample,
        // bhp_sample) values give a piecewise linear approximation to
        // the true inverse inflow function, at the same flo values as
        // the VFP data.
        //
        // Then we extract a piecewise linear approximation from the
        // multilinear fbhp() by evaluating it at the flo_sample
        // points, with fractions given by the frates(bhp_sample)
        // values.
        //
        // When we have both piecewise linear curves defined on the
        // same flo_sample points, it is easy to distinguish between
        // the 0, 1 or 2 solution cases, and obtain the right interval
        // in which to solve for the solution we want (with highest
        // flow in case of 2 solutions).

        // Make the fbhp() function.
        const auto& controls = well_ecl_.injectionControls(summary_state);
        const auto& table = vfp_properties_->getInj()->getTable(controls.vfp_table_number);
        const double vfp_ref_depth = table.getDatumDepth();
        const double rho = getRefDensity(); // Use the density at the top perforation.
        const double thp_limit = this->getTHPConstraint(summary_state);
        const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
        auto fbhp = [this, &controls, thp_limit, dp](const std::vector<double>& rates) {
            assert(rates.size() == 3);
            return this->vfp_properties_->getInj()
                    ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit) - dp;
        };

        // Make the flo() function.
        auto flo = [&table](const std::vector<double>& rates) {
            return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
        };

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

        // Get the flo samples, add extra samples at low rates and bhp
        // limit point if necessary.
        std::vector<double> flo_samples = table.getFloAxis();
        if (flo_samples[0] > 0.0) {
            const double f0 = flo_samples[0];
            flo_samples.insert(flo_samples.begin(), { f0/20.0, f0/10.0, f0/5.0, f0/2.0 });
        }
        const double flo_bhp_limit = flo(frates(controls.bhp_limit));
        if (flo_samples.back() < flo_bhp_limit) {
            flo_samples.push_back(flo_bhp_limit);
        }

        // Find bhp values for inflow relation corresponding to flo samples.
        std::vector<double> bhp_samples;
        for (double flo_sample : flo_samples) {
            if (flo_sample > flo_bhp_limit) {
                // We would have to go over the bhp limit to obtain a
                // flow of this magnitude. We associate all such flows
                // with simply the bhp limit. The first one
                // encountered is considered valid, the rest not. They
                // are therefore skipped.
                bhp_samples.push_back(controls.bhp_limit);
                break;
            }
            auto eq = [&flo, &frates, flo_sample](double bhp) {
                return flo(frates(bhp)) - flo_sample;
            };
            // TODO: replace hardcoded low/high limits.
            const double low = 10.0 * unit::barsa;
            const double high = 800.0 * unit::barsa;
            const int max_iteration = 100;
            const double flo_tolerance = 0.05 * std::fabs(flo_samples.back());
            int iteration = 0;
            try {
                const double solved_bhp = RegulaFalsiBisection<WarnAndContinueOnError>::
                        solve(eq, low, high, max_iteration, flo_tolerance, iteration);
                bhp_samples.push_back(solved_bhp);
            }
            catch (...) {
                // Use previous value (or max value if at start) if we failed.
                bhp_samples.push_back(bhp_samples.empty() ? low : bhp_samples.back());
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                        "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + name());
            }
        }

        // Find bhp values for VFP relation corresponding to flo samples.
        const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
        std::vector<double> fbhp_samples(num_samples);
        for (int ii = 0; ii < num_samples; ++ii) {
            fbhp_samples[ii] = fbhp(frates(bhp_samples[ii]));
        }
// #define EXTRA_THP_DEBUGGING
#ifdef EXTRA_THP_DEBUGGING
        std::string dbgmsg;
        dbgmsg += "flo: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(flo_samples[ii]);
        }
        dbgmsg += "\nbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(bhp_samples[ii]);
        }
        dbgmsg += "\nfbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(fbhp_samples[ii]);
        }
        OpmLog::debug(dbgmsg);
#endif // EXTRA_THP_DEBUGGING

        // Look for sign changes for the (fbhp_samples - bhp_samples) piecewise linear curve.
        // We only look at the valid
        int sign_change_index = -1;
        for (int ii = 0; ii < num_samples - 1; ++ii) {
            const double curr = fbhp_samples[ii] - bhp_samples[ii];
            const double next = fbhp_samples[ii + 1] - bhp_samples[ii + 1];
            if (curr * next < 0.0) {
                // Sign change in the [ii, ii + 1] interval.
                sign_change_index = ii; // May overwrite, thereby choosing the highest-flo solution.
            }
        }

        // Handle the no solution case.
        if (sign_change_index == -1) {
            return std::optional<double>();
        }

        // Solve for the proper solution in the given interval.
        auto eq = [&fbhp, &frates](double bhp) {
            return fbhp(frates(bhp)) - bhp;
        };
        // TODO: replace hardcoded low/high limits.
        const double low = bhp_samples[sign_change_index + 1];
        const double high = bhp_samples[sign_change_index];
        const int max_iteration = 100;
        const double bhp_tolerance = 0.01 * unit::barsa;
        int iteration = 0;
        if (low == high) {
            // We are in the high flow regime where the bhp_samples
            // are all equal to the bhp_limit.
            assert(low == controls.bhp_limit);
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }
        try {
            const double solved_bhp = RegulaFalsiBisection<WarnAndContinueOnError>::
                    solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
#ifdef EXTRA_THP_DEBUGGING
            OpmLog::debug("*****    " + name() + "    solved_bhp = " + std::to_string(solved_bhp)
                          + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
#endif // EXTRA_THP_DEBUGGING
            return solved_bhp;
        }
        catch (...) {
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }

    }





    template<typename TypeTag>
    double
    MultisegmentWell<TypeTag>::
    maxPerfPress(const Simulator& ebos_simulator) const
    {
        double max_pressure = 0.0;
        const int nseg = numberOfSegments();
        for (int seg = 0; seg < nseg; ++seg) {
            for (const int perf : segment_perforations_[seg]) {
                const int cell_idx = well_cells_[perf];
                const auto& int_quants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                const auto& fs = int_quants.fluidState();
                double pressure_cell = fs.pressure(FluidSystem::oilPhaseIdx).value();
                max_pressure = std::max(max_pressure, pressure_cell);
            }
        }
        return max_pressure;
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    pressureDropSpiralICD(const int seg) const
    {
        const SICD& sicd = segmentSet()[seg].spiralICD();

        const int seg_upwind = upwinding_segments_[seg];
        const std::vector<EvalWell>& phase_fractions = segment_phase_fractions_[seg_upwind];
        const std::vector<EvalWell>& phase_viscosities = segment_phase_viscosities_[seg_upwind];

        EvalWell water_fraction = 0.;
        EvalWell water_viscosity = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            water_fraction = phase_fractions[water_pos];
            water_viscosity = phase_viscosities[water_pos];
        }

        EvalWell oil_fraction = 0.;
        EvalWell oil_viscosity = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            oil_fraction = phase_fractions[oil_pos];
            oil_viscosity = phase_viscosities[oil_pos];
        }

        EvalWell gas_fraction = 0.;
        EvalWell gas_viscosity = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            gas_fraction = phase_fractions[gas_pos];
            gas_viscosity = phase_viscosities[gas_pos];
        }

        EvalWell density = segment_densities_[seg_upwind];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (seg != seg_upwind) {
            water_fraction.clearDerivatives();
            water_viscosity.clearDerivatives();
            oil_fraction.clearDerivatives();
            oil_viscosity.clearDerivatives();
            gas_fraction.clearDerivatives();
            gas_viscosity.clearDerivatives();
            density.clearDerivatives();
        }

        const EvalWell liquid_emulsion_viscosity = mswellhelpers::emulsionViscosity(water_fraction, water_viscosity,
                                                     oil_fraction, oil_viscosity, sicd);
        const EvalWell mixture_viscosity = (water_fraction + oil_fraction) * liquid_emulsion_viscosity + gas_fraction * gas_viscosity;

        const EvalWell reservoir_rate = segment_mass_rates_[seg] / density;

        const EvalWell reservoir_rate_icd = reservoir_rate * sicd.scalingFactor();

        const double viscosity_cali = sicd.viscosityCalibration();

        using MathTool = MathToolbox<EvalWell>;

        const double density_cali = sicd.densityCalibration();
        const EvalWell temp_value1 = MathTool::pow(density / density_cali, 0.75);
        const EvalWell temp_value2 = MathTool::pow(mixture_viscosity / viscosity_cali, 0.25);

        // formulation before 2016, base_strength is used
        // const double base_strength = sicd.strength() / density_cali;
        // formulation since 2016, strength is used instead
        const double strength = sicd.strength();

        const double sign = reservoir_rate_icd <= 0. ? 1.0 : -1.0;

        return sign * temp_value1 * temp_value2 * strength * reservoir_rate_icd * reservoir_rate_icd;
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    pressureDropAutoICD(const int seg, const UnitSystem& unit_system) const
    {
        const AutoICD& aicd = this->segmentSet()[seg].autoICD();

        const int seg_upwind = this->upwinding_segments_[seg];
        const std::vector<EvalWell>& phase_fractions = this->segment_phase_fractions_[seg_upwind];
        const std::vector<EvalWell>& phase_viscosities = this->segment_phase_viscosities_[seg_upwind];
        const std::vector<EvalWell>& phase_densities = this->segment_phase_densities_[seg_upwind];

        EvalWell water_fraction = 0.;
        EvalWell water_viscosity = 0.;
        EvalWell water_density = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            water_fraction = phase_fractions[water_pos];
            water_viscosity = phase_viscosities[water_pos];
            water_density = phase_densities[water_pos];
        }

        EvalWell oil_fraction = 0.;
        EvalWell oil_viscosity = 0.;
        EvalWell oil_density = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            oil_fraction = phase_fractions[oil_pos];
            oil_viscosity = phase_viscosities[oil_pos];
            oil_density = phase_densities[oil_pos];
        }

        EvalWell gas_fraction = 0.;
        EvalWell gas_viscosity = 0.;
        EvalWell gas_density = 0.;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            gas_fraction = phase_fractions[gas_pos];
            gas_viscosity = phase_viscosities[gas_pos];
            gas_density = phase_densities[gas_pos];
        }

        EvalWell density = segment_densities_[seg_upwind];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (seg != seg_upwind) {
            water_fraction.clearDerivatives();
            water_viscosity.clearDerivatives();
            water_density.clearDerivatives();
            oil_fraction.clearDerivatives();
            oil_viscosity.clearDerivatives();
            oil_density.clearDerivatives();
            gas_fraction.clearDerivatives();
            gas_viscosity.clearDerivatives();
            gas_density.clearDerivatives();
            density.clearDerivatives();
        }

        using MathTool = MathToolbox<EvalWell>;
        const EvalWell mixture_viscosity = MathTool::pow(water_fraction, aicd.waterViscExponent()) * water_viscosity
                                         + MathTool::pow(oil_fraction, aicd.oilViscExponent()) * oil_viscosity
                                         + MathTool::pow(gas_fraction, aicd.gasViscExponent()) * gas_viscosity;

        const EvalWell mixture_density = MathTool::pow(water_fraction, aicd.waterDensityExponent()) * water_density
                                       + MathTool::pow(oil_fraction, aicd.oilDensityExponent()) * oil_density
                                       + MathTool::pow(gas_fraction, aicd.gasDensityExponent()) * gas_density;

        const double rho_reference = aicd.densityCalibration();
        const double visc_reference = aicd.viscosityCalibration();
        const auto volume_rate_icd = this->segment_mass_rates_[seg] * aicd.scalingFactor() / mixture_density;
        const double sign = volume_rate_icd <= 0. ? 1.0 : -1.0;
        // convert 1 unit volume rate
        using M  = UnitSystem::measure;
        const double unit_volume_rate = unit_system.to_si(M::geometric_volume_rate, 1.);

        // TODO: we did not consider the maximum allowed rate here
        const auto result = sign / rho_reference * mixture_density * mixture_density
                          * MathTool::pow(visc_reference/mixture_viscosity, aicd.viscExponent())
                          * aicd.strength() * MathTool::pow( -sign * volume_rate_icd, aicd.flowRateExponent())
                          * std::pow(unit_volume_rate, (2. - aicd.flowRateExponent())) ;
        return result;
    }




    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::EvalWell
    MultisegmentWell<TypeTag>::
    pressureDropValve(const int seg) const
    {
        const Valve& valve = segmentSet()[seg].valve();

        const EvalWell& mass_rate = segment_mass_rates_[seg];
        const int seg_upwind = upwinding_segments_[seg];
        EvalWell visc = segment_viscosities_[seg_upwind];
        EvalWell density = segment_densities_[seg_upwind];
        // WARNING
        // We disregard the derivatives from the upwind density to make sure derivatives
        // wrt. to different segments dont get mixed.
        if (seg != seg_upwind) {
            visc.clearDerivatives();
            density.clearDerivatives();
        }

        const double additional_length = valve.pipeAdditionalLength();
        const double roughness = valve.pipeRoughness();
        const double diameter = valve.pipeDiameter();
        const double area = valve.pipeCrossArea();

        const EvalWell friction_pressure_loss =
            mswellhelpers::frictionPressureLoss(additional_length, diameter, area, roughness, density, mass_rate, visc);

        const double area_con = valve.conCrossArea();
        const double cv = valve.conFlowCoefficient();

        const EvalWell constriction_pressure_loss =
            mswellhelpers::valveContrictionPressureLoss(mass_rate, density, area_con, cv);

        const double sign = mass_rate <= 0. ? 1.0 : -1.0;
        return sign * (friction_pressure_loss + constriction_pressure_loss);
    }




    template<typename TypeTag>
    std::vector<double>
    MultisegmentWell<TypeTag>::
    computeCurrentWellRates(const Simulator& ebosSimulator,
                            DeferredLogger& deferred_logger) const
    {
        // Calculate the rates that follow from the current primary variables.
        std::vector<EvalWell> well_q_s(num_components_, 0.0);
        const bool allow_cf = getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);
        const int nseg = numberOfSegments();
        for (int seg = 0; seg < nseg; ++seg) {
            // calculating the perforation rate for each perforation that belongs to this segment
            const EvalWell seg_pressure = getSegmentPressure(seg);
            for (const int perf : segment_perforations_[seg]) {
                const int cell_idx = well_cells_[perf];
                const auto& int_quants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
                std::vector<EvalWell> mob(num_components_, 0.0);
                getMobility(ebosSimulator, perf, mob);
                const double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(int_quants, cell_idx);
                const double Tw = well_index_[perf] * trans_mult;
                std::vector<EvalWell> cq_s(num_components_, 0.0);
                EvalWell perf_press;
                double perf_dis_gas_rate = 0.;
                double perf_vap_oil_rate = 0.;
                computePerfRatePressure(int_quants, mob, Tw, seg, perf, seg_pressure, allow_cf, cq_s, perf_press, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);
                for (int comp = 0; comp < num_components_; ++comp) {
                    well_q_s[comp] += cq_s[comp];
                }
            }
        }
        std::vector<double> well_q_s_noderiv(well_q_s.size());
        for (int comp = 0; comp < num_components_; ++comp) {
            well_q_s_noderiv[comp] = well_q_s[comp].value();
        }
        return well_q_s_noderiv;
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeConnLevelProdInd(const typename MultisegmentWell<TypeTag>::FluidState& fs,
                            const std::function<double(const double)>& connPICalc,
                            const std::vector<EvalWell>& mobility,
                            double* connPI) const
    {
        const auto& pu = this->phaseUsage();
        const int   np = this->number_of_phases_;
        for (int p = 0; p < np; ++p) {
            // Note: E100's notion of PI value phase mobility includes
            // the reciprocal FVF.
            const auto connMob =
                mobility[ flowPhaseToEbosCompIdx(p) ].value()
                * fs.invB(flowPhaseToEbosPhaseIdx(p)).value();

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





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    computeConnLevelInjInd(const typename MultisegmentWell<TypeTag>::FluidState& fs,
                           const Phase preferred_phase,
                           const std::function<double(const double)>& connIICalc,
                           const std::vector<EvalWell>& mobility,
                           double* connII,
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
                             "Unsupported Injector Type ("
                             << static_cast<int>(preferred_phase)
                             << ") for well " << this->name()
                             << " during connection I.I. calculation",
                             deferred_logger);
        }

        const auto zero   = EvalWell { 0.0 };
        const auto mt     = std::accumulate(mobility.begin(), mobility.end(), zero);
        connII[phase_pos] = connIICalc(mt.value() * fs.invB(flowPhaseToEbosPhaseIdx(phase_pos)).value());
    }
} // namespace Opm
