/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#ifndef OPM_BLACKOILMULTISEGMENTMODEL_HEADER_INCLUDED
#define OPM_BLACKOILMULTISEGMENTMODEL_HEADER_INCLUDED

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/WellStateMultiSegment.hpp>
#include <opm/autodiff/WellMultiSegment.hpp>

namespace Opm {

    struct BlackoilMultiSegmentSolutionState : public DefaultBlackoilSolutionState
    {
        explicit BlackoilMultiSegmentSolutionState(const int np)
            : DefaultBlackoilSolutionState(np)
            , segp  ( ADB::null())
            , segqs ( ADB::null())
        {
        }
        ADB segp; // the segment pressures
        ADB segqs; // the segment phase rate in surface volume
    };

    /// A model implementation for three-phase black oil with support
    /// for multi-segment wells.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    /// \tparam  Grid            UnstructuredGrid or CpGrid.
    /// \tparam  Implementation  Provides concrete state types.
    template<class Grid>
    class BlackoilMultiSegmentModel : public BlackoilModelBase<Grid, BlackoilMultiSegmentModel<Grid>>
    {
    public:

        typedef BlackoilModelBase<Grid, BlackoilMultiSegmentModel<Grid> > Base; // base class
        typedef typename Base::ReservoirState ReservoirState;
        typedef typename Base::WellState WellState;
        typedef BlackoilMultiSegmentSolutionState SolutionState;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilMultiSegmentModel(const typename Base::ModelParameters&  param,
                          const Grid&                     grid ,
                          const BlackoilPropsAdInterface& fluid,
                          const DerivedGeology&           geo  ,
                          const RockCompressibility*      rock_comp_props,
                          const Wells*                    wells,
                          const NewtonIterationBlackoilInterface& linsolver,
                          Opm::EclipseStateConstPtr eclState,
                          const bool has_disgas,
                          const bool has_vapoil,
                          const bool terminal_output,
                          const std::vector<WellMultiSegmentConstPtr>& wells_multisegment);

        /// Called once before each time step.
        /// \param[in] dt                     time step size
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void prepareStep(const double dt,
                         ReservoirState& reservoir_state,
                         WellState& well_state);

        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] dt                     time step size
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void afterStep(const double dt,
                       ReservoirState& reservoir_state,
                       WellState& well_state) {};


        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        void assemble(const ReservoirState& reservoir_state,
                      WellState& well_state,
                      const bool initial_assembly) {};


        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void updateState(const V& dx,
                         ReservoirState& reservoir_state,
                         WellState& well_state) {};

    protected:
     /*
        // ---------  Types and enums  ---------
        // using Base::DataBlock;
        // using Base::ReservoirResidualQuant;
     */
        // ---------  Data members  ---------

        // For the non-segmented well, it should be the density with AVG or SEG way.
        // while usually SEG way
        using Base::well_perforation_densities_; //Density of each well perforation
        using Base::pvdt_;
        using Base::geo_;
        using Base::active_;


        // Diff to the pressure of the related segment.
        // When the well is a usual well, the bhp will be the pressure of the top segment
        // For mutlti-segmented wells, only AVG is allowed.
        // For non-segmented wells, typically SEG is used. AVG way might not have been
        // implemented yet.

        // Diff to bhp for each well perforation. only for usual wells.
        // For segmented wells, they are zeros.
        using Base::well_perforation_pressure_diffs_; // Diff to bhp for each well perforation.

        // ADB version of the densities, when using AVG way, the calculation of the density and hydrostatic head
        // is implicit
        ADB well_perforation_densities_adb_;

        // ADB version. Eventually, only ADB version will be kept.
        ADB well_perforation_pressure_diffs_adb_;

        // Pressure correction due to the different depth of the perforation
        // and the cell center of the grid block
        // For the non-segmented wells, since the perforation are forced to be
        // at the center of the grid cell, it should be ZERO.
        // It should only apply to the mutli-segmented wells.
        V well_perforation_pressure_cell_diffs_;
        ADB well_perforation_pressure_cell_diffs_adb_;

        // Pressure correction due to the depth differennce between segment depth and perforation depth.
        // TODO: It should be able to be merge as a part of the perforation_pressure_diffs_.
        ADB well_perforations_segment_pressure_diffs_;

        // the average of the fluid densities in the grid block
        // which is used to calculate the hydrostatic head correction due to the depth difference of the perforation
        // and the cell center of the grid block
        V well_perforation_cell_densities_;
        ADB well_perforation_cell_densities_adb_;

        V well_perforatoin_cell_pressure_diffs_;

        const std::vector<WellMultiSegmentConstPtr> wells_multisegment_;

        // return wells object
        // TODO: remove this wells structure
        using Base::wells;
        using Base::updatePrimalVariableFromState;

        const std::vector<WellMultiSegmentConstPtr>& wellsMultiSegment() const { return wells_multisegment_; }

        SolutionState
        variableState(const ReservoirState& x,
                      const WellState& xw) const {};


        void updateWellControls(WellState& xw) const {};

        void updateWellState(const V& dwells,
                             WellState& well_state) {};

        std::vector<V>
        variableStateInitials(const ReservoirState& x,
                              const WellState& xw) const {};

        void
        variableWellStateInitials(const WellState& xw,
                                  std::vector<V>& vars0) const;

        void computeWellConnectionPressures(const SolutionState& state,
                                            const WellState& xw) {};

        void
        computeWellFlux(const SolutionState& state,
                        const std::vector<ADB>& mob_perfcells,
                        const std::vector<ADB>& b_perfcells,
                        V& aliveWells,
                        std::vector<ADB>& cq_s);

        void
        solveWellEq(const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    SolutionState& state,
                    WellState& well_state);

        void
        updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                         const SolutionState& state,
                                         WellState& xw) {};

        void
        addWellFluxEq(const std::vector<ADB>& cq_s,
                      const SolutionState& state) {};

        void
        addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                           const SolutionState& state,
                                           const WellState& xw) {};

        void
        addWellControlEq(const SolutionState& state,
                         const WellState& xw,
                         const V& aliveWells) {};

        void
        makeConstantState(SolutionState& state) const;

        void
        variableStateExtractWellsVars(const std::vector<int>& indices,
                                      std::vector<ADB>& vars,
                                      SolutionState& state) const;


/*

        const Grid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const RockCompressibility*      rock_comp_props_;
        const Wells*                    wells_;
        // FOR TEMPORARY
        // SHOUlD BE A REFERENCE
        VFPProperties                   vfp_properties_;
        const NewtonIterationBlackoilInterface&    linsolver_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active phases. Maps active -> canonical phase indices.
        const std::vector<int>          canph_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const bool has_disgas_;
        const bool has_vapoil_;

        ModelParameters                 param_;
        bool use_threshold_pressure_;
        bool wells_active_;
        V threshold_pressures_by_interior_face_;

        std::vector<ReservoirResidualQuant> rq_;
        std::vector<PhasePresence> phaseCondition_;
        V isRs_;
        V isRv_;
        V isSg_;

        std::vector<int>         primalVariable_;
        V pvdt_;

        // ---------  Protected methods  ---------

        /// Access the most-derived class used for
        /// static polymorphism (CRTP).
        Implementation& asImpl()
        {
            return static_cast<Implementation&>(*this);
        }

        /// Access the most-derived class used for
        /// static polymorphism (CRTP).
        const Implementation& asImpl() const
        {
            return static_cast<const Implementation&>(*this);
        }

        // return true if wells are available in the reservoir
        bool wellsActive() const { return wells_active_; }
        // return true if wells are available on this process
        bool localWellsActive() const { return wells_ ? (wells_->number_of_wells > 0 ) : false; }



        void
        variableReservoirStateInitials(const ReservoirState& x,
                                       std::vector<V>& vars0) const;

        std::vector<int>
        variableStateIndices() const;

        SolutionState
        variableStateExtractVars(const ReservoirState& x,
                                 const std::vector<int>& indices,
                                 std::vector<ADB>& vars) const;

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );


        void
        assembleMassBalanceEq(const SolutionState& state);




        bool getWellConvergence(const int iteration);

        bool isVFPActive() const;

        std::vector<ADB>
        computePressures(const ADB& po,
                         const ADB& sw,
                         const ADB& so,
                         const ADB& sg) const;

        V
        computeGasPressure(const V& po,
                           const V& sw,
                           const V& so,
                           const V& sg) const;

        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;

        void
        computeMassFlux(const int               actph ,
                        const V&                transi,
                        const ADB&              kr    ,
                        const ADB&              p     ,
                        const SolutionState&    state );

        void applyThresholdPressures(ADB& dp);

        ADB
        fluidViscosity(const int               phase,
                       const ADB&              p    ,
                       const ADB&              temp ,
                       const ADB&              rs   ,
                       const ADB&              rv   ,
                       const std::vector<PhasePresence>& cond) const;

        ADB
        fluidReciprocFVF(const int               phase,
                         const ADB&              p    ,
                         const ADB&              temp ,
                         const ADB&              rs   ,
                         const ADB&              rv   ,
                         const std::vector<PhasePresence>& cond) const;

        ADB
        fluidDensity(const int  phase,
                     const ADB& b,
                     const ADB& rs,
                     const ADB& rv) const;

        V
        fluidRsSat(const V&                p,
                   const V&                so,
                   const std::vector<int>& cells) const;

        ADB
        fluidRsSat(const ADB&              p,
                   const ADB&              so,
                   const std::vector<int>& cells) const;

        V
        fluidRvSat(const V&                p,
                   const V&                so,
                   const std::vector<int>& cells) const;

        ADB
        fluidRvSat(const ADB&              p,
                   const ADB&              so,
                   const std::vector<int>& cells) const;

        ADB
        poroMult(const ADB& p) const;

        ADB
        transMult(const ADB& p) const;

        const std::vector<PhasePresence>
        phaseCondition() const {return phaseCondition_;}

        void
        classifyCondition(const ReservoirState& state);


        /// update the primal variable for Sg, Rv or Rs. The Gas phase must
        /// be active to call this method.
        void
        updatePrimalVariableFromState(const ReservoirState& state);

        /// Update the phaseCondition_ member based on the primalVariable_ member.
        /// Also updates isRs_, isRv_ and isSg_;
        void
        updatePhaseCondFromPrimalVariable();

        /// \brief Compute the reduction within the convergence check.
        /// \param[in] B     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[in] tempV A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. tempV.col(i) contains the
        ///                   values
        ///                  for phase i.
        /// \param[in] R     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[out] R_sum An array of size MaxNumPhases where entry i contains the sum
        ///                   of R for the phase i.
        /// \param[out] maxCoeff An array of size MaxNumPhases where entry i contains the
        ///                   maximum of tempV for the phase i.
        /// \param[out] B_avg An array of size MaxNumPhases where entry i contains the average
        ///                   of B for the phase i.
        /// \param[out] maxNormWell The maximum of the well equations for each phase.
        /// \param[in]  nc    The number of cells of the local grid.
        /// \param[in]  nw    The number of wells on the local grid.
        /// \return The total pore volume over all cells.
        double
        convergenceReduction(const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases>& B,
                             const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases>& tempV,
                             const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases>& R,
                             std::array<double,MaxNumPhases>& R_sum,
                             std::array<double,MaxNumPhases>& maxCoeff,
                             std::array<double,MaxNumPhases>& B_avg,
                             std::vector<double>& maxNormWell,
                             int nc,
                             int nw) const;

        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; } */

    };

    /// Providing types by template specialisation of ModelTraits for BlackoilMultiSegmentModel.
    template <class GridT>
    struct ModelTraits< BlackoilMultiSegmentModel<GridT> >
    {
        typedef BlackoilState ReservoirState;
        typedef WellStateMultiSegment WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef BlackoilMultiSegmentSolutionState SolutionState;
    };




} // namespace Opm

#include "BlackoilMultiSegmentModel_impl.hpp"

#endif // OPM_BLACKOILMULTISEGMENTMODEL_HEADER_INCLUDED
