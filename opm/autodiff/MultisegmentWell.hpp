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


#ifndef OPM_MULTISEGMENTWELL_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_HEADER_INCLUDED


#include <opm/autodiff/WellInterface.hpp>

namespace Opm
{

    template<typename TypeTag>
    class MultisegmentWell: public WellInterface<TypeTag>
    {
    public:
        typedef WellInterface<TypeTag> Base;

        // TODO: the WellState does not have any information related to segments
        using typename Base::WellState;

        // TODO: for now, not considering the polymer, solvent and so on to simplify the development process.
        // TODO: should I begin with the old primary variable or the new fraction based variable systems?
        // Let us begin with the new one
        enum WellVariablePositions {
            GTotal = 0,
            WFrac = 1,
            GFrac = 2,
            SPres = 3
        };

        ///  the number of well equations // TODO: it should have a more general strategy for it
        static const int numWellEq = 4;

        using typename Base::Scalar;
        using typename Base::ConvergenceReport;
        /// the number of reservior equations
        using Base::numEq;

        /// the matrix and vector types for the reservoir
        using typename Base::Mat;
        using typename Base::BVector;
        using typename Base::Eval;

        // sparsity pattern for the matrices
        // [A C^T    [x       =  [ res
        //  B  D ]   x_well]      res_well]

        // the vector type for the res_well and x_well
        typedef Dune::FieldVector<Scalar, numWellEq> VectorBlockWellType;
        typedef Dune::BlockVector<VectorBlockWellType> BVectorWell;

        // the matrix type for the diagonal matrix D
        typedef Dune::FieldMatrix<Scalar, numWellEq, numWellEq > DiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix <DiagMatrixBlockWellType> DiagMatWell;

        // the matrix type for the non-diagonal matrix B and C^T
        typedef Dune::FieldMatrix<Scalar, numWellEq, numEq>  OffDiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellType> OffDiagMatWell;

        // TODO: for more efficient implementation, we should have EvalReservoir, EvalWell, and EvalRerservoirAndWell
        //                                                         EvalR (Eval), EvalW, EvalRW
        // TODO: for now, we only use one type to save some implementation efforts, while improve later.
        typedef DenseAd::Evaluation<double, /*size=*/numEq + numWellEq> EvalWell;

        MultisegmentWell(const Well* well, const int time_step, const Wells* wells);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<bool>* active_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells);


        virtual void initPrimaryVariablesEvaluation() const;

        virtual void assembleWellEq(Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    bool only_wells);

        /// updating the well state based the control mode specified with current
        // TODO: later will check wheter we need current
        virtual void updateWellStateWithTarget(const int current,
                                               WellState& xw) const;

        // TODO: this should go to the WellInterface, while updateWellStateWithTarget
        // will need touch different types of well_state, we will see.
        virtual void updateWellControl(WellState& xw,
                                       wellhelpers::WellSwitchingLogger& logger) const;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(Simulator& ebosSimulator,
                                                     const std::vector<double>& B_avg,
                                                     const ModelParameters& param) const;

        /// computing the accumulation term for later use in well mass equations
        virtual void computeAccumWell();

        virtual void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                    const WellState& xw);

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x, const ModelParameters& param,
                                                           WellState& well_state) const;

        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials) const;

        virtual void updatePrimaryVariables(const WellState& well_state) const;

        virtual void solveEqAndUpdateWellState(const ModelParameters& param,
                                               WellState& well_state); // const?

        /// number of segments for this well
        /// int number_of_segments_;
        int numberOfSegments() const;

        int numberOfPerforations() const;

    protected:
        int number_segments_;

        // components of the pressure drop to be included
        WellSegment::CompPressureDropEnum compPressureDrop() const;
        // multi-phase flow model
        WellSegment::MultiPhaseModelEnum multiphaseModel() const;

        // get the SegmentSet from the well_ecl_
        const SegmentSet& segmentSet() const;

        using Base::well_ecl_;
        using Base::number_of_perforations_; // TODO: can use well_ecl_?

        using Base::well_cells_; // TODO: are the perforation orders same with StandardWell or Wells?
        using Base::well_index_;

        using Base::well_controls_;

        // TODO: trying to use the information from the Well opm-parser as much
        // as possible, it will possibly be re-implemented later for efficiency reason.


        // indices of the gird blocks that segments locate at.
        // TODO: the grid cell related to a segment should be calculated based on the location
        //       of the segment node.
        //       As the current temporary solution, the grid cell related to a segment determined by the
        //       first perforation cell related to the segment.
        //       when no perforation is related to the segment, use it outlet segment's cell.
        std::vector<int> segment_cell_;

        // Things are easy to get from SegmentSet
        // segment_volume_, segment_cross_area_, segment_length_(total length), segment_depth_
        // segment_internal_diameter_, segment_roughness_
        // outlet_segment_., in the outlet_segment, we store the ID of the segment, we will need to use numberToLocation to get
        // their location in the segmentSet

        // segment number is an ID of the segment, it is specified in the deck
        // get the loation of the segment with a segment number in the segmentSet
        int numberToLocation(const int segment_number) const;

        // TODO, the following should go to a class for computing purpose
        // two off-diagonal matrices
        mutable OffDiagMatWell duneB_;
        mutable OffDiagMatWell duneC_;
        // diagonal matrix for the well
        mutable DiagMatWell invDuneD_;

        // several vector used in the matrix calculation
        mutable BVectorWell Bx_;
        mutable BVectorWell invDrw_;
        mutable BVector scaleAddRes_;

        // residuals of the well equations
        BVectorWell resWell_;

        // the values for the primary varibles
        // based on different solutioin strategies, the wells can have different primary variables
        // TODO: should we introduce a data structure for segment to simplify this?
        // or std::vector<std::vector<double> >
        mutable std::vector<std::array<double, numWellEq> > primary_variables_;

        // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
        mutable std::vector<std::array<EvalWell, numWellEq> > primary_variables_evaluation_;

        // protected functions
        // EvalWell getBhp(); this one should be something similar to getSegmentPressure();
        // EvalWell getQs(); this one should be something similar to getSegmentRates()
        // EValWell wellVolumeFractionScaled, wellVolumeFraction, wellSurfaceVolumeFraction ... these should have different names, and probably will be needed.
        // bool crossFlowAllowed(const Simulator& ebosSimulator) const; probably will be needed
        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

        // updating the well_state based on well solution dwells
        void updateWellState(const BVectorWell& dwells,
                             const BlackoilModelParameters& param,
                             WellState& well_state) const;

        // void computePerfRate() will be a key function here.
    };

}

#include "MultisegmentWell_impl.hpp"

#endif // OPM_MULTISEGMENTWELL_HEADER_INCLUDED
