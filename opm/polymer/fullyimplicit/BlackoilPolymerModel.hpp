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

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        void assemble(const ReservoirState& reservoir_state,
                      WellState& well_state,
                      const bool initial_assembly);
        // void assemble(const PolymerBlackoilState& reservoir_state,
        //               WellStateFullyImplicitBlackoilPolymer& well_state,
        //               const bool initial_assembly);
        /// \brief Compute the residual norms of the mass balance for each phase,
        /// the well flux, and the well equation.
        /// \return a vector that contains for each phase the norm of the mass balance
        /// and afterwards the norm of the residual of the well flux and the well equation.
        std::vector<double> computeResidualNorms() const;

        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        V solveJacobianSystem() const;

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

    protected:

        // ---------  Types and enums  ---------

        typedef typename Base::SolutionState SolutionState;
        typedef typename Base::DataBlock DataBlock;

        // ---------  Data members  ---------

        const PolymerPropsAd& polymer_props_ad_;
        const bool has_polymer_;
        const int  poly_pos_;
        V cmax_;

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

        SolutionState
        constantState(const ReservoirState& x,
                      const WellState& xw) const;

        void
        makeConstantState(SolutionState& state) const;

        SolutionState
        variableState(const ReservoirState& x,
                      const WellState& xw) const;

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );

        void computeWellConnectionPressures(const SolutionState& state,
                                            const WellState& xw);

        void
        addWellControlEq(const SolutionState& state,
                         const WellState& xw,
                         const V& aliveWells);

        void
        addWellEq(const SolutionState& state,
                  WellState& xw,
                  V& aliveWells);

        void updateWellControls(WellState& xw) const;

        std::vector<ADB>
        computePressures(const SolutionState& state) const;

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

        void
        computeCmax(ReservoirState& state);

        ADB
        computeMc(const SolutionState& state) const;

        void applyThresholdPressures(ADB& dp);

        ADB
        fluidViscosity(const int               phase,
                       const ADB&              p    ,
                       const ADB&              temp ,
                       const ADB&              rs   ,
                       const ADB&              rv   ,
                       const std::vector<PhasePresence>& cond,
                       const std::vector<int>& cells) const;

        ADB
        fluidReciprocFVF(const int               phase,
                         const ADB&              p    ,
                         const ADB&              temp ,
                         const ADB&              rs   ,
                         const ADB&              rv   ,
                         const std::vector<PhasePresence>& cond,
                         const std::vector<int>& cells) const;

        ADB
        fluidDensity(const int               phase,
                     const ADB&              p    ,
                     const ADB&              temp ,
                     const ADB&              rs   ,
                     const ADB&              rv   ,
                     const std::vector<PhasePresence>& cond,
                     const std::vector<int>& cells) const;

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

        void
        classifyCondition(const SolutionState&        state,
                          std::vector<PhasePresence>& cond ) const;

        const std::vector<PhasePresence>
        phaseCondition() const {return this->phaseCondition_;}

        void
        classifyCondition(const ReservoirState&        state);


        /// update the primal variable for Sg, Rv or Rs. The Gas phase must
        /// be active to call this method.
        void
        updatePrimalVariableFromState(const ReservoirState&        state);

        /// Update the phaseCondition_ member based on the primalVariable_ member.
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

        double dpMaxRel() const { return this->param_.dp_max_rel_; }
        double dsMax() const { return this->param_.ds_max_; }
        double drMaxRel() const { return this->param_.dr_max_rel_; }
        double maxResidualAllowed() const { return this->param_.max_residual_allowed_; }

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
