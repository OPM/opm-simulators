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
                                    const Wells&                    wells);

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
        };

        enum { Water = BlackoilPropsAdInterface::Water,
               Oil   = BlackoilPropsAdInterface::Oil  ,
               Gas   = BlackoilPropsAdInterface::Gas  };

        // Member data
        const UnstructuredGrid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const Wells&                    wells_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active faces. Maps active -> canonical phase indices.
        const std::vector<int>          canph_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const M                         grav_;

        std::vector<ReservoirResidualQuant> rq_;

        struct {
            std::vector<ADB> reservoir;
        } residual_;

        // Private methods.
        void
        allocateResidual();

        SolutionState
        constantState(const BlackoilState& x);

        SolutionState
        variableState(const BlackoilState& x);

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );

        void
        assemble(const V&             dtpv,
                 const BlackoilState& x   ,
                 const WellState&     xw  );

        std::vector<ADB>
        computeRelPerm(const SolutionState& state);

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
    };
} // namespace Opm


#endif // OPM_FULLYIMPLICITBLACKOILSOLVER_HEADER_INCLUDED
