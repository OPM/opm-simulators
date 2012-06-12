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

#ifndef OPM_INCOMPTPFA_HEADER_INCLUDED
#define OPM_INCOMPTPFA_HEADER_INCLUDED


#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <vector>

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class IncompPropertiesInterface;
    class RockCompressibility;
    class LinearSolverInterface;
    class TwophaseState;
    class WellState;

    /// Encapsulating a tpfa pressure solver for the incompressible-fluid case.
    /// Supports gravity, wells controlled by bhp or reservoir rates,
    /// boundary conditions and simple sources as driving forces.
    /// Rock compressibility can be included, and necessary nonlinear
    /// iterations are handled.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class IncompTpfa
    {
    public:
	/// Construct solver.
        /// \param[in] grid             A 2d or 3d grid.
        /// \param[in] props            Rock and fluid properties.
        /// \param[in] rock_comp_props  Rock compressibility properties.
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
        /// \param[in] src              Source terms. May be empty().
        /// \param[in] bcs              Boundary conditions, treat as all noflow if null.
	IncompTpfa(const UnstructuredGrid& grid,
                   const IncompPropertiesInterface& props,
                   const RockCompressibility* rock_comp_props,
                   LinearSolverInterface& linsolver,
                   const double residual_tol,
                   const double change_tol,
                   const int maxiter,
                   const double* gravity,
                   const Wells* wells,
		   const std::vector<double>& src,
		   const FlowBoundaryConditions* bcs);

	/// Destructor.
	~IncompTpfa();

        /// Solve the pressure equation. If there is no pressure
        /// dependency introduced by rock compressibility effects,
        /// the equation is linear, and it is solved directly.
        /// Otherwise, the nonlinear equations ares solved by a
        /// Newton-Raphson scheme.
        /// May throw an exception if the number of iterations
        /// exceed maxiter (set in constructor).
        void solve(const double dt,
                   TwophaseState& state,
                   WellState& well_state);


        /// Expose read-only reference to internal half-transmissibility.
        const std::vector<double>& getHalfTrans() const { return htrans_; }

    private:
        // Solve with no rock compressibility (linear eqn).
        void solveIncomp(const double dt,
                         TwophaseState& state,
                         WellState& well_state);
        // Solve with rock compressibility (nonlinear eqn).
        void solveRockComp(const double dt,
                           TwophaseState& state,
                           WellState& well_state);
        // Helper functions.
        void computePerSolveDynamicData(const double dt,
                                        const TwophaseState& state,
                                        const WellState& well_state);
        void computePerIterationDynamicData(const double dt,
                                            const TwophaseState& state,
                                            const WellState& well_state);
        void assemble(const double dt,
                      const TwophaseState& state,
                      const WellState& well_state);
        void solveIncrement();
        double residualNorm() const;
        double incrementNorm() const;
	void computeResults(TwophaseState& state,
                            WellState& well_state) const;

    private:
        // ------ Data that will remain unmodified after construction. ------
	const UnstructuredGrid& grid_;
        const IncompPropertiesInterface& props_;
        const RockCompressibility* rock_comp_props_;
        const LinearSolverInterface& linsolver_;
        const double residual_tol_;
        const double change_tol_;
        const int maxiter_;
        const double* gravity_; // May be NULL
        const Wells* wells_;    // May be NULL, outside may modify controls (only) between calls to solve().
        const std::vector<double>& src_;
        const FlowBoundaryConditions* bcs_;
	std::vector<double> htrans_;
	std::vector<double> gpress_;
        std::vector<int> allcells_;

        // ------ Data that will be modified for every solve. ------
	std::vector<double> trans_ ;
        std::vector<double> wdp_;
        std::vector<double> totmob_;
        std::vector<double> omega_;
	std::vector<double> gpress_omegaweighted_;
        std::vector<double> initial_porevol_;
        struct ifs_tpfa_forces forces_;

        // ------ Data that will be modified for every solver iteration. ------
        std::vector<double> porevol_;
        std::vector<double> rock_comp_;
        std::vector<double> pressures_;

        // ------ Internal data for the ifs_tpfa solver. ------
	struct ifs_tpfa_data* h_;
    };

} // namespace Opm

#endif // OPM_INCOMPTPFA_HEADER_INCLUDED
