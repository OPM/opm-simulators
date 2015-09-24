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

        friend Base;

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
                      const bool initial_assembly);


        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /* void updateState(const V& dx,
                         ReservoirState& reservoir_state,
                         WellState& well_state) {}; */
        using Base::numPhases;

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
        using Base::rq_;
        using Base::fluid_;
        using Base::terminal_output_;
        using Base::grid_;
        using Base::canph_;


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
        using Base::wellsActive;
        using Base::phaseCondition;
        using Base::fluidRvSat;
        using Base::fluidRsSat;
        using Base::fluidDensity;

        const std::vector<WellMultiSegmentConstPtr>& wellsMultiSegment() const { return wells_multisegment_; }

        void updateWellControls(WellState& xw) const;

        using Base::variableState;

        // void updateWellState(const V& dwells,
        //                      WellState& well_state) {};

        void
        variableWellStateInitials(const WellState& xw,
                                  std::vector<V>& vars0) const;

        void computeWellConnectionPressures(const SolutionState& state,
                                            const WellState& xw);

        void
        computeWellFlux(const SolutionState& state,
                        const std::vector<ADB>& mob_perfcells,
                        const std::vector<ADB>& b_perfcells,
                        V& aliveWells,
                        std::vector<ADB>& cq_s);

        void
        updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                         const SolutionState& state,
                                         WellState& xw);

        /* void
        addWellFluxEq(const std::vector<ADB>& cq_s,
                      const SolutionState& state) {}; */

        /* void
        addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                           const SolutionState& state,
                                           const WellState& xw) {}; */

        /* void
        addWellControlEq(const SolutionState& state,
                         const WellState& xw,
                         const V& aliveWells) {}; */

        void
        makeConstantState(SolutionState& state) const;

        void
        variableStateExtractWellsVars(const std::vector<int>& indices,
                                      std::vector<ADB>& vars,
                                      SolutionState& state) const;


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
