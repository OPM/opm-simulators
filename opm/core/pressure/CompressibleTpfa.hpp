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

#ifndef OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
#define OPM_COMPRESSIBLETPFA_HEADER_INCLUDED


#include <vector>

struct UnstructuredGrid;
struct cfs_tpfa_res_data;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class BlackoilState;
    class BlackoilPropertiesInterface;
    class RockCompressibility;
    class LinearSolverInterface;
    class WellState;

    /// Encapsulating a tpfa pressure solver for the compressible-fluid case.
    /// Supports gravity, wells and simple sources as driving forces.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class CompressibleTpfa
    {
    public:
        /// Construct solver.
        /// \param[in] grid             A 2d or 3d grid.
        /// \param[in] props            Rock and fluid properties.
        /// \param[in] rock_comp_props  Rock compressibility properties. May be null.
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
        CompressibleTpfa(const UnstructuredGrid& grid,
                         const BlackoilPropertiesInterface& props,
                         const RockCompressibility* rock_comp_props,
                         const LinearSolverInterface& linsolver,
                         const double residual_tol,
                         const double change_tol,
                         const int maxiter,
                         const double* gravity,
                         const Wells* wells);

        /// Destructor.
        ~CompressibleTpfa();

        /// Solve the pressure equation by Newton-Raphson scheme.
        /// May throw an exception if the number of iterations
        /// exceed maxiter (set in constructor).
        void solve(const double dt,
                   BlackoilState& state,
                   WellState& well_state);

        /// @brief After solve(), was the resulting pressure singular.
        /// Returns true if the pressure is singular in the following
        /// sense: if everything is incompressible and there are no
        /// pressure conditions, the absolute values of the pressure
        /// solution are arbitrary. (But the differences in pressure
        /// are significant.)
        bool singularPressure() const;

    private:
        virtual void computePerSolveDynamicData(const double dt,
                                                const BlackoilState& state,
                                                const WellState& well_state);
        void computePerIterationDynamicData(const double dt,
                                            const BlackoilState& state,
                                            const WellState& well_state);
        virtual void computeCellDynamicData(const double dt,
                                            const BlackoilState& state,
                                            const WellState& well_state);
        void computeFaceDynamicData(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state);
        void computeWellDynamicData(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state);
        void assemble(const double dt,
                      const BlackoilState& state,
                      const WellState& well_state);
        void solveIncrement();
        double residualNorm() const;
        double incrementNorm() const;
        void computeResults(BlackoilState& state,
                            WellState& well_state) const;
    protected:
        void computeWellPotentials(const BlackoilState& state);

        // ------ Data that will remain unmodified after construction. ------
        const UnstructuredGrid& grid_;
        const BlackoilPropertiesInterface& props_;
        const RockCompressibility* rock_comp_props_;
        const LinearSolverInterface& linsolver_;
        const double residual_tol_;
        const double change_tol_;
        const int maxiter_;
        const double* gravity_; // May be NULL
        const Wells* wells_;    // May be NULL, outside may modify controls (only) between calls to solve().
        std::vector<double> htrans_;
        std::vector<double> trans_ ;
        std::vector<int> allcells_;

        // ------ Internal data for the cfs_tpfa_res solver. ------
        struct cfs_tpfa_res_data* h_;

        // ------ Data that will be modified for every solve. ------
        std::vector<double> wellperf_gpot_;
        std::vector<double> initial_porevol_;

        // ------ Data that will be modified for every solver iteration. ------
        std::vector<double> cell_A_;
        std::vector<double> cell_dA_;
        std::vector<double> cell_viscosity_;
        std::vector<double> cell_phasemob_;
        std::vector<double> cell_voldisc_;
        std::vector<double> face_A_;
        std::vector<double> face_phasemob_;
        std::vector<double> face_gravcap_;
        std::vector<double> wellperf_A_;
        std::vector<double> wellperf_phasemob_;
        std::vector<double> porevol_;   // Only modified if rock_comp_props_ is non-null.
        std::vector<double> rock_comp_; // Empty unless rock_comp_props_ is non-null.
        // The update to be applied to the pressures (cell and bhp).
        std::vector<double> pressure_increment_;
        // True if the matrix assembled would be singular but for the
        // adjustment made in the cfs_*_assemble() calls. This happens
        // if everything is incompressible and there are no pressure
        // conditions.
        bool singular_;
    };

} // namespace Opm


#endif // OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
