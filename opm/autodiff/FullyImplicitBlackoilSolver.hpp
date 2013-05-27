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
    class LinearSolverInterface;
    class BlackoilState;
    class WellState;


    /// A fully implicit TPFA-based solver for the black-oil problem.
    class FullyImplicitBlackoilSolver
    {
    public:
        FullyImplicitBlackoilSolver(const UnstructuredGrid&         grid ,
                                    const BlackoilPropsAdInterface& fluid,
                                    const DerivedGeology&           geo  ,
                                    const Wells&                    wells,
                                    const LinearSolverInterface&    linsolver);

        /// Take a single forward step, modifiying
        ///   state.pressure()
        ///   state.faceflux()
        ///   state.saturation()
        ///   state.surfacevol()
        void
        step(const double   dt    ,
             BlackoilState& state ,
             WellState&     wstate);

    private:
        // Types and enums
        typedef AutoDiff::ForwardBlock<double> ADB;
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
            ADB              Rs;
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

        // The mass_balance vector has one element for each active phase,
        // each of which has size equal to the number of cells.
        // The well_eq has size equal to the number of wells.
        struct {
            std::vector<ADB> mass_balance;
            ADB well_eq;
        } residual_;

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

        void solveJacobianSystem(BlackoilState& x,
                                 WellState& xw) const;

        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;

        std::vector<ADB>
        computeRelPermWells(const SolutionState& state,
                            const DataBlock& well_s,
                            const std::vector<int>& well_cells) const;

        void
        computeMassFlux(const int               actph ,
                        const V&                transi,
                        const std::vector<ADB>& kr    ,
                        const SolutionState&    state );

        double
        residualNorm() const;

        ADB
        fluidViscosity(const int               phase,
                       const ADB&              p    ,
                       const std::vector<int>& cells) const;

        ADB
        fluidReciprocFVF(const int               phase,
                         const ADB&              p    ,
                         const std::vector<int>& cells) const;

        ADB
        fluidDensity(const int               phase,
                     const ADB&              p    ,
                     const std::vector<int>& cells) const;

        ADB
        fluidRsMax(const ADB&              p,
                   const std::vector<int>& cells) const;
    };
} // namespace Opm


#endif // OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED
