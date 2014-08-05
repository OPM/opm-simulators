/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/LinearisedBlackoilResidual.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>

struct UnstructuredGrid;
struct Wells;

namespace Opm {

    namespace parameter { class ParameterGroup; }
    class DerivedGeology;
    class RockCompressibility;
    class NewtonIterationBlackoilInterface;
    class BlackoilState;
    class WellStateFullyImplicitBlackoil;


    /// A fully implicit solver for the black-oil problem.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil (but not vice versa). It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    template<class T>
    class FullyImplicitBlackoilSolver
    {
    public:
        /// \brief The type of the grid that we use.
        typedef T Grid;
        /// Construct a solver. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param       parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] linsolver        linear solver
        FullyImplicitBlackoilSolver(const parameter::ParameterGroup& param,
                                    const Grid&                     grid ,
                                    const BlackoilPropsAdInterface& fluid,
                                    const DerivedGeology&           geo  ,
                                    const RockCompressibility*      rock_comp_props,
                                    const Wells&                    wells,
                                    const NewtonIterationBlackoilInterface& linsolver,
                                    const bool has_disgas,
                                    const bool has_vapoil );

        /// Take a single forward step, modifiying
        ///   state.pressure()
        ///   state.faceflux()
        ///   state.saturation()
        ///   state.gasoilratio()
        ///   wstate.bhp()
        /// \param[in] dt        time step size
        /// \param[in] state     reservoir state
        /// \param[in] wstate    well state
        void
        step(const double   dt    ,
             BlackoilState& state ,
             WellStateFullyImplicitBlackoil&     wstate);

    private:
        // Types and enums
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;

        struct ReservoirResidualQuant {
            ReservoirResidualQuant();
            std::vector<ADB> accum; // Accumulations
            ADB              mflux; // Mass flux (surface conditions)
            ADB              b;     // Reciprocal FVF
            ADB              head;  // Pressure drop across int. interfaces
            ADB              mob;   // Phase mobility (per cell)
        };

        struct SolutionState {
            SolutionState(const int np);
            ADB              pressure;
            std::vector<ADB> saturation;
            ADB              rs;
            ADB              rv;
            ADB              qs;
            ADB              bhp;       
        };

        struct WellOps {
            WellOps(const Wells& wells);
            M w2p;              // well -> perf (scatter)
            M p2w;              // perf -> well (gather)
        };

        enum { Water = BlackoilPropsAdInterface::Water,
               Oil   = BlackoilPropsAdInterface::Oil  ,
               Gas   = BlackoilPropsAdInterface::Gas  };

        // the Newton relaxation type
        enum RelaxType { DAMPEN, SOR };
        enum PrimalVariables { Sg = 0, RS = 1, RV = 2 };

        // Member data
        const Grid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const RockCompressibility*      rock_comp_props_;
        const Wells&                    wells_;
        const NewtonIterationBlackoilInterface&    linsolver_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active faces. Maps active -> canonical phase indices.
        const std::vector<int>          canph_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const WellOps                   wops_;
        const bool has_disgas_;
        const bool has_vapoil_;
        double                          dp_max_rel_;
        double                          ds_max_;
        double                          drs_max_rel_;
        enum RelaxType                  relax_type_;
        double                          relax_max_;
        double                          relax_increment_;
        double                          relax_rel_tol_;
        int                             max_iter_;

        std::vector<ReservoirResidualQuant> rq_;
        std::vector<PhasePresence> phaseCondition_;
        V well_perforation_pressure_diffs_; // Diff to bhp for each well perforation.

        LinearisedBlackoilResidual residual_;

        std::vector<int>         primalVariable_;

        // Private methods.
        SolutionState
        constantState(const BlackoilState& x,
                      const WellStateFullyImplicitBlackoil& xw);

        SolutionState
        variableState(const BlackoilState& x,
                      const WellStateFullyImplicitBlackoil& xw);

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );

        void computeWellConnectionPressures(const SolutionState& state,
                                            const WellStateFullyImplicitBlackoil& xw);

        void
        addWellControlEq(const SolutionState& state,
                         const WellStateFullyImplicitBlackoil& xw,
                         const V& aliveWells);

        void
        addWellEq(const SolutionState& state,
                  WellStateFullyImplicitBlackoil& xw,
                  V& aliveWells);

        void updateWellControls(ADB& bhp,
                                ADB& well_phase_flow_rate,
                                WellStateFullyImplicitBlackoil& xw) const;

        void
        assemble(const V&             dtpv,
                 const BlackoilState& x,
                 WellStateFullyImplicitBlackoil& xw);

        V solveJacobianSystem() const;

        void updateState(const V& dx,
                         BlackoilState& state,
                         WellStateFullyImplicitBlackoil& well_state);

        std::vector<ADB>
        computePressures(const SolutionState& state) const;

        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;

        std::vector<ADB>
        computeRelPermWells(const SolutionState& state,
                            const DataBlock& well_s,
                            const std::vector<int>& well_cells) const;

        void
        computeMassFlux(const int               actph ,
                        const V&                transi,
                        const ADB&              kr    ,
                        const ADB&              p     ,
                        const SolutionState&    state );

        double
        residualNorm() const;

        std::vector<double> residuals() const;

        ADB
        fluidViscosity(const int               phase,
                       const ADB&              p    ,
                       const ADB&              rs   ,
                       const ADB&              rv   ,
                       const std::vector<PhasePresence>& cond,
                       const std::vector<int>& cells) const;

        ADB
        fluidReciprocFVF(const int               phase,
                         const ADB&              p    ,
                         const ADB&              rs   ,
                         const ADB&              rv   ,
                         const std::vector<PhasePresence>& cond,
                         const std::vector<int>& cells) const;

        ADB
        fluidDensity(const int               phase,
                     const ADB&              p    ,
                     const ADB&              rs   ,
                     const ADB&              rv   ,
                     const std::vector<PhasePresence>& cond,
                     const std::vector<int>& cells) const;

        V
        fluidRsSat(const V&                p,
                   const std::vector<int>& cells) const;

        ADB
        fluidRsSat(const ADB&              p,
                   const std::vector<int>& cells) const;

        V
        fluidRvSat(const V&                p,
                   const std::vector<int>& cells) const;

        ADB
        fluidRvSat(const ADB&              p,
                   const std::vector<int>& cells) const;

        ADB
        poroMult(const ADB& p) const;

        ADB
        transMult(const ADB& p) const;

        void
        classifyCondition(const SolutionState&        state,
                          std::vector<PhasePresence>& cond ) const;

        const std::vector<PhasePresence>
        phaseCondition() const {return phaseCondition_;}

        void
        classifyCondition(const BlackoilState&        state);


        /// update the primal variable for Sg, Rv or Rs. The Gas phase must
        /// be active to call this method.
        void
        updatePrimalVariableFromState(const BlackoilState&        state);

        /// Update the phaseCondition_ member based on the primalVariable_ member.
        void
        updatePhaseCondFromPrimalVariable();

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        bool getConvergence(const double dt);

        void detectNewtonOscillations(const std::vector<std::vector<double>>& residual_history,
                                      const int it, const double relaxRelTol,
                                      bool& oscillate, bool& stagnate) const;

        void stablizeNewton(V& dx, V& dxOld, const double omega, const RelaxType relax_type) const;

        double dpMaxRel() const { return dp_max_rel_; }
        double dsMax() const { return ds_max_; }
        double drsMaxRel() const { return drs_max_rel_; }
        enum RelaxType relaxType() const { return relax_type_; }
        double relaxMax() const { return relax_max_; };
        double relaxIncrement() const { return relax_increment_; };
        double relaxRelTol() const { return relax_rel_tol_; };
        double maxIter() const { return max_iter_; }

    };
} // namespace Opm

#include "FullyImplicitBlackoilSolver_impl.hpp"

#endif // OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED
