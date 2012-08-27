/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_INCOMPTPFAPOLYMER_HEADER_INCLUDED
#define OPM_INCOMPTPFAPOLYMER_HEADER_INCLUDED


#include <opm/core/pressure/CompressibleTpfa.hpp>

#include <vector>

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class BlackoilState;
    class PolymerBlackoilState;
    class RockCompressibility;
    class PolymerProperties;
    class LinearSolverInterface;
    class PolymerState;
    class WellState;

    /// Encapsulating a tpfa pressure solver for the compressible-fluid case  with polymer.
    /// Supports gravity, wells controlled by bhp or reservoir rates,
    /// boundary conditions and simple sources as driving forces.
    /// Rock compressibility can be included, and necessary nonlinear
    /// iterations are handled.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class CompressibleTpfaPolymer : public CompressibleTpfa
    {
    public:

	/// Construct solver, possibly with rock compressibility.
        /// \param[in] grid             A 2d or 3d grid.
        /// \param[in] props            Rock and fluid properties.
        /// \param[in] rock_comp_props  Rock compressibility properties. May be null.
        /// \param[in] poly_props       Polymer properties.
        /// \param[in] linsolver        Linear solver to use.
        /// \param[in] residual_tol     Solution accepted if inf-norm of residual is smaller.
        /// \param[in] change_tol       Solution accepted if inf-norm of change in pressure is smaller.
        /// \param[in] maxiter          Maximum acceptable number of iterations.
        /// \param[in] gravity          Gravity vector. If non-null, the array should
        ///                             have D elements.
        /// \param[in] wells            The wells argument. Will be used in solution,
        ///                             is ignored if NULL.
        ///                             Note: this class observes the well object, and
        ///                                   makes the assumption that the well topology
        ///                                   and completions does not change during the
        ///                                   run. However, controls (only) are allowed
        ///                                   to change.
	CompressibleTpfaPolymer(const UnstructuredGrid& grid,
                                const BlackoilPropertiesInterface& props,
                                const RockCompressibility* rock_comp_props,
                                const PolymerProperties& poly_props,
                                const LinearSolverInterface& linsolver,
                                const double residual_tol,
                                const double change_tol,
                                const int maxiter,
                                const double* gravity,
                                const Wells* wells);

        /// Solve the pressure equation. The nonlinear equations ares solved by a
        /// Newton-Raphson scheme.
        /// May throw an exception if the number of iterations
        /// exceed maxiter (set in constructor).
        void solve(const double dt,
                   PolymerBlackoilState& state,
                   WellState& well_state);

    private:
        virtual void computeCellDynamicData(const double dt,
                                            const BlackoilState& state,
                                            const WellState& well_state);

        virtual void computePerSolveDynamicData(const double dt,
                                                const BlackoilState& state,
                                                const WellState& well_state);


    private:
        // ------ Data that will remain unmodified after construction. ------
        const PolymerProperties& poly_props_;
        // ------ Data that will be updated every solve() call. ------
        const std::vector<double>* c_;
        const std::vector<double>* cmax_;
        std::vector<double> cell_eff_viscosity_;
        std::vector<double> cell_relperm_;
        std::vector<double> cell_eff_relperm_;
    };

} // namespace Opm

#endif // OPM_INCOMPTPFAPOLYMER_HEADER_INCLUDED
