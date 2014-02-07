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

struct UnstructuredGrid;
struct Wells;

namespace Opm {

    class DerivedGeology;
    class RockCompressibility;
    class LinearSolverInterface;
    class BlackoilState;
    class WellState;


    /// A fully implicit solver for the black-oil problem.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil (but not vice versa). It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    class FullyImplicitBlackoilSolver
    {
    public:
        /// Construct a solver. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] linsolver        linear solver
        FullyImplicitBlackoilSolver(const UnstructuredGrid&         grid ,
                                    const BlackoilPropsAdInterface& fluid,
                                    const DerivedGeology&           geo  ,
                                    const RockCompressibility*      rock_comp_props,
                                    const Wells&                    wells,
                                    const LinearSolverInterface&    linsolver);

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
             WellState&     wstate);

        /// Typedef used throughout the solver, in the public section
        /// since it is used in the definition of Residual.
        typedef AutoDiffBlock<double> ADB;

        /// Residual structure of the fully implicit solver.
        /// All equations are given as AD types, with multiple
        /// jacobian blocks corresponding to the primary unknowns. The
        /// primary unknowns are for a three-phase simulation, in order:
        ///    p    (pressure)
        ///    sw   (water saturation)
        ///    xvar (gas saturation, gas-oil ratio or oil-gas ratio)
        ///    qs   (well outflows by well and phase)
        ///    bhp  (bottom hole pressures)
        /// In the above, the xvar variable will have a different
        /// meaning from cell to cell, corresponding to the state in
        /// that cell (saturated, undersaturated oil or undersaturated
        /// gas). In a two-phase simulation, either sw or xvar is not
        /// used, depending on which face is missing.
        struct Residual {
            /// The mass_balance vector has one element for each
            /// active phase, each of which has size equal to the
            /// number of cells. Each mass balance equation is given
            /// in terms of surface volumes.
            std::vector<ADB> mass_balance;
            /// The well_flux_eq has size equal to the number of wells
            /// times the number of phases. It contains the well flow
            /// equations, relating the total well flows to
            /// bottom-hole pressures and reservoir conditions.
            ADB well_flux_eq;
            /// The well_eq has size equal to the number of wells. It
            /// contains the well control equations, that is for each
            /// well either a rate specification or bottom hole
            /// pressure specification.
            ADB well_eq;
        };


    private:
        // Types and enums
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

        // Member data
        const UnstructuredGrid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const RockCompressibility*      rock_comp_props_;
        const Wells&                    wells_;
        const LinearSolverInterface&    linsolver_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active faces. Maps active -> canonical phase indices.
        const std::vector<int>          canph_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const WellOps                   wops_;
        const M                         grav_;

        std::vector<ReservoirResidualQuant> rq_;
        std::vector<PhasePresence> phaseCondition_;

        Residual residual_;

        // Private methods.
        SolutionState
        constantState(const BlackoilState& x,
                      const WellState&     xw);

        SolutionState
        variableState(const BlackoilState& x,
                      const WellState&     xw);

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );

        void
        assemble(const V&             dtpv,
                 const BlackoilState& x   ,
                 const WellState&     xw  );

        V solveJacobianSystem() const;

        void updateState(const V& dx,
                         BlackoilState& state,
                         WellState& well_state);

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


    };
} // namespace Opm


#endif // OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED
