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
    :Basse(well, time_step, wells)
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
        // More facility will be required from the opm-parser to make all the situations
        // make sense.
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
        invDuneD_.setSize(numSegment(), numSegment(), 100000);
        duneB_.setSize(numSegment(), num_cells, number_of_perforations_);
        duneC_.setSize(numSegment(), num_cells, number_of_perforations_);

        // we need to add the off diagonal ones
        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (int perf = 0 ; perf < ]; ++perf) { // the segments hold some perforations
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        // make the B^T matrix
        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            // the number of the row corresponds to the segment number now.
            for (int perf = wells().well_connpos[row.index()] ; perf < wells().well_connpos[row.index()+1]; ++perf) {
                const int cell_idx = wells().well_cells[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize( nw );

        // resize temporary class variables
        Cx_.resize( duneC_.N() );
        invDrw_.resize( invDuneD_.N() );
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (int seg = 0; seg < numSegment(); ++seg) {
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
                              WellState& xw) const
    {
        // TODO: it can be challenging, when updating the segment and perforation related
        // well rates will be okay.
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    updateWellControl(WellState& xw,
                      wellhelpers::WellSwitchingLogger& logger) const
    {
        // TODO: it will be very similar to the StandardWell, while updateWellStateWithTarget will be chanlleging.
    }





    template<typename TypeTag>
    typename MultisegmentWell<TypeTag>::ConvergenceReport
    MultisegmentWell<TypeTag>::
    getWellConvergence(Simulator& ebosSimulator,
                       const std::vector<double>& B_avg,
                       const ModelParameters& param) const
    {
        // TODO: it will be very similar
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
                                   const WellState& xw)
    {
        // TODO: the name of the function need to change.
        // it will be calculating the pressure difference between the perforation and grid cells
        // With MS well, the depth of the perforation is not necessarily the center of the grid cells.
    }





    template<typename TypeTag>
    void
    MultisegmentWell<TypeTag>::
    void apply(const BVector& x, BVector& Ax) const
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
        return well_ecl_->getSegmentSet(time_step);
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
