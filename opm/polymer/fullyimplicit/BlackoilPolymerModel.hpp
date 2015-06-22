/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#ifndef OPM_BLACKOILPOLYMERMODEL_HEADER_INCLUDED
#define OPM_BLACKOILPOLYMERMODEL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/fullyimplicit/WellStateFullyImplicitBlackoilPolymer.hpp>

namespace Opm {

    /// A model implementation for three-phase black oil with polymer.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa, with polymer
    /// in the water phase. It uses an industry-standard TPFA
    /// discretization with per-phase upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    template<class Grid>
    class BlackoilPolymerModel : public BlackoilModelBase<Grid, BlackoilPolymerModel<Grid> >
    {
    public:

        // ---------  Types and enums  ---------

        typedef BlackoilModelBase<Grid, BlackoilPolymerModel<Grid> > Base;
        typedef typename Base::ReservoirState ReservoirState;
        typedef typename Base::WellState WellState;
        // The next line requires C++11 support available in g++ 4.7.
        // friend Base;
        friend class BlackoilModelBase<Grid, BlackoilPolymerModel<Grid> >;

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] linsolver        linear solver
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] has_polymer      turn on polymer feature
        /// \param[in] has_plyshlog     true when PLYSHLOG keyword available
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilPolymerModel(const typename Base::ModelParameters&   param,
                             const Grid&                             grid,
                             const BlackoilPropsAdInterface&         fluid,
                             const DerivedGeology&                   geo,
                             const RockCompressibility*              rock_comp_props,
                             const PolymerPropsAd&                   polymer_props_ad,
                             const Wells*                            wells,
                             const NewtonIterationBlackoilInterface& linsolver,
                             const bool                              has_disgas,
                             const bool                              has_vapoil,
                             const bool                              has_polymer,
                             const bool                              has_plyshlog,
                             const std::vector<double>&              wells_rep_radius,
                             const std::vector<double>&              wells_perf_length,
                             const bool                              terminal_output);

        /// Called once before each time step.
        /// \param[in] dt                     time step size
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void prepareStep(const double dt,
                         ReservoirState& reservoir_state,
                         WellState& well_state);

        /// Called once after each time step.
        /// \param[in] dt                     time step size
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void afterStep(const double dt,
                       ReservoirState& reservoir_state,
                       WellState& well_state);

        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void updateState(const V& dx,
                         ReservoirState& reservoir_state,
                         WellState& well_state);

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   dt          timestep length
        /// \param[in]   iteration   current iteration number
        bool getConvergence(const double dt, const int iteration);

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        void assemble(const ReservoirState& reservoir_state,
                      WellState& well_state,
                      const bool initial_assembly);


    protected:

        // ---------  Types and enums  ---------

        typedef typename Base::SolutionState SolutionState;
        typedef typename Base::DataBlock DataBlock;
        enum { Concentration = CanonicalVariablePositions::Next };

        // ---------  Data members  ---------

        const PolymerPropsAd& polymer_props_ad_;
        const bool has_polymer_;
        const bool has_plyshlog_;
        const int  poly_pos_;
        V cmax_;

        // representative radius and perforation length of well perforations
        // to be used in shear-thinning computation.
        std::vector<double> wells_rep_radius_;
        std::vector<double> wells_perf_length_;

        // shear-thinning factor for cell faces
        std::vector<double> shear_mult_faces_;
        // shear-thinning factor for well perforations
        std::vector<double> shear_mult_wells_;

        // Need to declare Base members we want to use here.
        using Base::grid_;
        using Base::fluid_;
        using Base::geo_;
        using Base::rock_comp_props_;
        using Base::wells_;
        using Base::linsolver_;
        using Base::active_;
        using Base::canph_;
        using Base::cells_;
        using Base::ops_;
        using Base::wops_;
        using Base::has_disgas_;
        using Base::has_vapoil_;
        using Base::param_;
        using Base::use_threshold_pressure_;
        using Base::threshold_pressures_by_interior_face_;
        using Base::rq_;
        using Base::phaseCondition_;
        using Base::well_perforation_pressure_diffs_;
        using Base::residual_;
        using Base::terminal_output_;
        using Base::primalVariable_;
        using Base::pvdt_;

        // ---------  Protected methods  ---------

        // Need to declare Base members we want to use here.
        using Base::wellsActive;
        using Base::wells;
        using Base::variableState;
        using Base::computePressures;
        using Base::computeGasPressure;
        using Base::applyThresholdPressures;
        using Base::fluidViscosity;
        using Base::fluidReciprocFVF;
        using Base::fluidDensity;
        using Base::fluidRsSat;
        using Base::fluidRvSat;
        using Base::poroMult;
        using Base::transMult;
        using Base::updatePrimalVariableFromState;
        using Base::updatePhaseCondFromPrimalVariable;
        using Base::dpMaxRel;
        using Base::dsMax;
        using Base::drMaxRel;
        using Base::maxResidualAllowed;

        using Base::updateWellControls;
        using Base::computeWellConnectionPressures;
        using Base::addWellControlEq;
        using Base::computeRelPerm;


        void
        makeConstantState(SolutionState& state) const;

        std::vector<V>
        variableStateInitials(const ReservoirState& x,
                              const WellState& xw) const;

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

        void
        addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                           const SolutionState& state,
                                           WellState& xw);

        void
        computeMassFlux(const int               actph ,
                        const V&                transi,
                        const ADB&              kr    ,
                        const ADB&              p     ,
                        const SolutionState&    state );

        void
        computeCmax(ReservoirState& state);

        ADB
        computeMc(const SolutionState& state) const;

        const std::vector<PhasePresence>
        phaseCondition() const {return this->phaseCondition_;}

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
        /// \param[in]  nc    The number of cells of the local grid.
        /// \return The total pore volume over all cells.
        double
        convergenceReduction(const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& B,
                             const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& tempV,
                             const Eigen::Array<double, Eigen::Dynamic, MaxNumPhases+1>& R,
                             std::array<double,MaxNumPhases+1>& R_sum,
                             std::array<double,MaxNumPhases+1>& maxCoeff,
                             std::array<double,MaxNumPhases+1>& B_avg,
                             std::vector<double>& maxNormWell,
                             int nc,
                             int nw) const;

        /// Computing the water velocity without shear-thinning for the cell faces.
        /// The water velocity will be used for shear-thinning calculation.
        void computeWaterShearVelocityFaces(const V& transi, const std::vector<ADB>& kr,
                                            const std::vector<ADB>& phasePressure, const SolutionState& state,
                                            std::vector<double>& water_vel, std::vector<double>& visc_mult);

        /// Computing the water velocity without shear-thinning for the well perforations based on the water flux rate.
        /// The water velocity will be used for shear-thinning calculation.
        void computeWaterShearVelocityWells(const SolutionState& state, WellState& xw, const ADB& cq_sw,
                                            std::vector<double>& water_vel_wells, std::vector<double>& visc_mult_wells);

    };



    /// Need to include concentration in our state variables, otherwise all is as
    /// the default blackoil model.
    struct BlackoilPolymerSolutionState : public DefaultBlackoilSolutionState
    {
        explicit BlackoilPolymerSolutionState(const int np)
            : DefaultBlackoilSolutionState(np),
              concentration( ADB::null())
        {
        }
        ADB concentration;
    };



    /// Providing types by template specialisation of ModelTraits for BlackoilPolymerModel.
    template <class Grid>
    struct ModelTraits< BlackoilPolymerModel<Grid> >
    {
        typedef PolymerBlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoilPolymer WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef BlackoilPolymerSolutionState SolutionState;
    };

} // namespace Opm

#include "BlackoilPolymerModel_impl.hpp"


#endif // OPM_BLACKOILPOLYMERMODEL_HEADER_INCLUDED
