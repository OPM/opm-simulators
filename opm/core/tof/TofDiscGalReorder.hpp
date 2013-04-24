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

#ifndef OPM_TOFDISCGALREORDER_HEADER_INCLUDED
#define OPM_TOFDISCGALREORDER_HEADER_INCLUDED

#include <opm/core/transport/reorder/ReorderSolverInterface.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <map>
#include <ostream>

struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;
    class VelocityInterpolationInterface;
    class DGBasisInterface;
    namespace parameter { class ParameterGroup; }
    template <typename T> class SparseTable;

    /// Implements a discontinuous Galerkin solver for
    /// (single-phase) time-of-flight using reordering.
    /// The equation solved is:
    ///     \f[v \cdot \nabla\tau = \phi\f]
    /// in which \f$ v \f$ is the fluid velocity, \f$ \tau \f$ is time-of-flight and
    /// \f$ \phi \f$ is the porosity. This is a boundary value problem, and
    /// \f$ \tau \f$ is specified to be zero on all inflow boundaries.
    /// The user may specify the polynomial degree of the basis function space
    /// used, but only degrees 0 and 1 are supported so far.
    class TofDiscGalReorder : public ReorderSolverInterface
    {
    public:
        /// Construct solver.
        /// \param[in] grid      A 2d or 3d grid.
        /// \param[in] param     Parameters for the solver.
        ///                      The following parameters are accepted (defaults):\n
        ///   - \c dg_degree (0)                           -- Polynomial degree of basis functions.
        ///   - \c use_tensorial_basis (false)             -- Use tensor-product basis, interpreting dg_degree as
        ///                                                   bi/tri-degree not total degree.
        ///   - \c use_cvi (false)                         -- Use ECVI velocity interpolation.
        ///   - \c use_limiter (false)                     -- Use a slope limiter. If true, the next three parameters are used.
        ///   - \c limiter_relative_flux_threshold (1e-3)  -- Ignore upstream fluxes below this threshold,
        ///                                                   relative to total cell flux.
        ///   - \c limiter_method ("MinUpwindFace")        -- Limiter method used. Accepted methods are:
        ///             - MinUpwindFace              -- Limit cell tof to >= inflow face tofs.
        ///             - MinUpwindAverage           -- Limit cell tof to >= inflow cell average tofs.
        ///   - \c limiter_usage ("DuringComputations")    -- Usage pattern for limiter. Accepted choices are:
        ///             - DuringComputations         -- Apply limiter to cells as they are computed,
        ///                                             so downstream cells' solutions may be affected
        ///                                             by limiting in upstream cells.
        ///             - AsPostProcess              -- Apply in dependency order, but only after
        ///                                             computing (unlimited) solution.
        ///             - AsSimultaneousPostProcess  -- Apply to each cell independently, using un-
        ///                                             limited solution in neighbouring cells.
        TofDiscGalReorder(const UnstructuredGrid& grid,
                          const parameter::ParameterGroup& param);


        /// Solve for time-of-flight.
        /// \param[in]  darcyflux         Array of signed face fluxes.
        /// \param[in]  porevolume        Array of pore volumes.
        /// \param[in]  source            Source term. Sign convention is:
        ///                                 (+) inflow flux,
        ///                                 (-) outflow flux.
        /// \param[out] tof_coeff         Array of time-of-flight solution coefficients.
        ///                               The values are ordered by cell, meaning that
        ///                               the K coefficients corresponding to the first
        ///                               cell come before the K coefficients corresponding
        ///                               to the second cell etc.
        ///                               K depends on degree and grid dimension.
        void solveTof(const double* darcyflux,
                      const double* porevolume,
                      const double* source,
                      std::vector<double>& tof_coeff);

        /// Solve for time-of-flight and a number of tracers.
        /// \param[in]  darcyflux         Array of signed face fluxes.
        /// \param[in]  porevolume        Array of pore volumes.
        /// \param[in]  source            Source term. Sign convention is:
        ///                                 (+) inflow flux,
        ///                                 (-) outflow flux.
        /// \param[in]  tracerheads       Table containing one row per tracer, and each
        ///                               row contains the source cells for that tracer.
        /// \param[out] tof_coeff         Array of time-of-flight solution coefficients.
        ///                               The values are ordered by cell, meaning that
        ///                               the K coefficients corresponding to the first
        ///                               cell comes before the K coefficients corresponding
        ///                               to the second cell etc.
        ///                               K depends on degree and grid dimension.
        /// \param[out] tracer_coeff      Array of tracer solution coefficients. N*K per cell,
        ///                               where N is equal to tracerheads.size(). All K coefs
        ///                               for a tracer are consecutive, and all tracers' coefs
        ///                               for a cell come before those for the next cell.
        void solveTofTracer(const double* darcyflux,
                            const double* porevolume,
                            const double* source,
                            const SparseTable<int>& tracerheads,
                            std::vector<double>& tof_coeff,
                            std::vector<double>& tracer_coeff);

    private:
        virtual void solveSingleCell(const int cell);
        virtual void solveMultiCell(const int num_cells, const int* cells);

    private:
        // Disable copying and assignment.
        TofDiscGalReorder(const TofDiscGalReorder&);
        TofDiscGalReorder& operator=(const TofDiscGalReorder&);

        // Data members
        const UnstructuredGrid& grid_;
        boost::shared_ptr<VelocityInterpolationInterface> velocity_interpolation_;
        bool use_cvi_;
        bool use_limiter_;
        double limiter_relative_flux_threshold_;
        enum LimiterMethod { MinUpwindFace, MinUpwindAverage };
        LimiterMethod limiter_method_;
        enum LimiterUsage { DuringComputations, AsPostProcess, AsSimultaneousPostProcess };
        LimiterUsage limiter_usage_;
        const double* darcyflux_;   // one flux per grid face
        const double* porevolume_;  // one volume per cell
        const double* source_;      // one volumetric source term per cell
        boost::shared_ptr<DGBasisInterface> basis_func_;
        double* tof_coeff_;
        // For tracers.
        double* tracer_coeff_;
        int num_tracers_;
        enum { NoTracerHead = -1 };
        std::vector<int> tracerhead_by_cell_;
        // Used by solveSingleCell().
        std::vector<double> rhs_;   // single-cell right-hand-sides
        std::vector<double> jac_;   // single-cell jacobian
        std::vector<double> orig_rhs_;   // single-cell right-hand-sides (copy)
        std::vector<double> orig_jac_;   // single-cell jacobian (copy)
        std::vector<double> coord_;
        mutable std::vector<double> basis_;
        mutable std::vector<double> basis_nb_;
        std::vector<double> grad_basis_;
        std::vector<double> velocity_;
        int num_singlesolves_;
        // Used by solveMultiCell():
        double gauss_seidel_tol_;
        int num_multicell_;
        int max_size_multicell_;
        int max_iter_multicell_;

        // Private methods

        // Apply some limiter, writing to array tof
        // (will read data from tof_coeff_, it is ok to call
        //  with tof_coeff as tof argument.
        void applyLimiter(const int cell, double* tof);
        void applyMinUpwindLimiter(const int cell, const bool face_min, double* tof);
        void applyLimiterAsPostProcess();
        void applyLimiterAsSimultaneousPostProcess();
        double totalFlux(const int cell) const;
        double minCornerVal(const int cell, const int face) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTRACERTOFDISCGAL_HEADER_INCLUDED
