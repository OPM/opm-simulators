/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
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


#ifndef OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/IncompPropsAdInterface.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


struct UnstructuredGrid;
struct Wells;

namespace Opm {
    class LinearSolverInterface;
    class TwophaseState;
    class WellState;
    
    /// A fully implicit solver for incompressible oil-water problem.
    ///
    /// The simulator is capable of handling incompressible oil-water
    /// problems.It uses an industry-standard TPFA discretization 
    /// with per-phase upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    class FullyImplicitTwoPhaseSolver
    {
    public:
        /// Construct a solver. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] linsolver        linear solver
        /// \param[in] wells            well structure
        /// \param[in] gravity			gravity
        FullyImplicitTwoPhaseSolver(const UnstructuredGrid&        grid,
                                    const IncompPropsAdInterface&  fluid,
                                    const LinearSolverInterface&   linsolver,
                                    const Wells&                   wells,
                                    const double*                  gravity);

        /// Take a single forward step, modifiying
        ///   state.pressure()
        ///   state.faceflux()
        ///   state.saturation()
        ///   wstate.bhp()
        /// \param[in] dt        time step size
        /// \param[in] state     reservoir state
        /// \param[in] well_tate well state
        void step(const double   dt,
                  TwophaseState& state,
                  WellState&     well_state);
    private:
		// Types
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;

        struct SolutionState {
            SolutionState(const int np);
            ADB              pressure;
            std::vector<ADB> saturation;
            ADB              qs;
            ADB              bhp;
        };
        struct WellOps {
            WellOps(const Wells& wells);
            M w2p;              // well -> perf (scatter)
            M p2w;              // perf -> well (gather)
        };

        const UnstructuredGrid&         grid_;
        const IncompPropsAdInterface&   fluid_;
        const LinearSolverInterface&    linsolver_;
        const Wells&                    wells_;
        const double*                   gravity_;
        const std::vector<int>          cells_;
        HelperOps                       ops_;
        const WellOps                   wops_;
        std::vector<ADB>                mob_;
       
        struct {
            std::vector<ADB>    mass_balance;
            ADB                 well_eq;
            ADB                 well_flux_eq;
        } residual_;
       
        SolutionState
        constantState(const TwophaseState& x,
                      const WellState&     xw);
        SolutionState
        variableState(const TwophaseState& x,
                      const WellState&     xw);
        void
        assemble(const V&               pvdt,
                 const SolutionState&   old_state,
                 const TwophaseState&   x,
                 const WellState&       xw);

        V solveJacobianSystem() const;
        void updateState(const V&             dx,
                         TwophaseState&       x,
                         WellState&           xw) const;

        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;

        V
        transmissibility() const;

        ADB
        computeFracFlow(const int phase);

        ADB 
        accumSource(const int phase,
                    const std::vector<ADB>& kr,
                    const std::vector<double>& src) const;

		std::vector<ADB>
		computePressures(const SolutionState& state) const;

        ADB
        computeMassFlux(const int               phase,
                        const V&                trans,
                        const std::vector<ADB>& kr,
						const ADB&				phasePress,
                        const SolutionState&    state);

        double
        residualNorm() const;

        ADB
        fluidDensity(const int phase,
                     const ADB p) const;

        ADB
        rockPorosity(const ADB& p) const;

        ADB
        rockPermeability(const ADB& p) const;

        ADB
        transMult(const ADB& p) const;
    };
} // namespace Opm

#endif// OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
