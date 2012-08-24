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

#ifndef OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
#define OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED

#include <opm/core/transport/reorder/TransportModelInterface.hpp>
#include <vector>
#include <map>
#include <ostream>
struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    /// Implements a reordering transport solver for incompressible two-phase flow.
    class TransportModelTwophase : public TransportModelInterface
    {
    public:
        /// Construct solver.
        /// \param[in] grid      A 2d or 3d grid.
        /// \param[in] props     Rock and fluid properties.
        /// \param[in] tol       Tolerance used in the solver.
        /// \param[in] maxit     Maximum number of non-linear iterations used.
        TransportModelTwophase(const UnstructuredGrid& grid,
                               const Opm::IncompPropertiesInterface& props,
                               const double tol,
                               const int maxit);

        /// Solve for saturation at next timestep.
        /// \param[in] darcyflux         Array of signed face fluxes.
        /// \param[in] porevolume        Array of pore volumes.
        /// \param[in] source            Transport source term.
        /// \param[in] dt                Time step.
        /// \param[in, out] saturation   Phase saturations.
        void solve(const double* darcyflux,
                   const double* porevolume,
                   const double* source,
                   const double dt,
                   std::vector<double>& saturation);

        /// Initialise quantities needed by gravity solver.
        /// \param[in] grav    gravity vector
        void initGravity(const double* grav);

        /// Solve for gravity segregation.
        /// This uses a column-wise nonlinear Gauss-Seidel approach.
        /// It assumes that the input columns contain cells in a single
        /// vertical stack, that do not interact with other columns (for
        /// gravity segregation.
        /// \param[in] columns           Vector of cell-columns.
        /// \param[in] porevolume        Array of pore volumes.
        /// \param[in] dt                Time step.
        /// \param[in, out] saturation   Phase saturations.
        void solveGravity(const std::vector<std::vector<int> >& columns,
                          const double* porevolume,
                          const double dt,
                          std::vector<double>& saturation);

        //// Return the number of iterations used by the reordering solver.
        //// \param[out] vector of iteration per cell
        const std::vector<int>& getReorderIterations() const;

    private:
        virtual void solveSingleCell(const int cell);
        virtual void solveMultiCell(const int num_cells, const int* cells);

        void solveSingleCellGravity(const std::vector<int>& cells,
                                    const int pos,
                                    const double* gravflux);
        int solveGravityColumn(const std::vector<int>& cells);
    private:
        const UnstructuredGrid& grid_;
        const IncompPropertiesInterface& props_;
        const double* visc_;
        std::vector<double> smin_;
        std::vector<double> smax_;
        double tol_;
        double maxit_;

        const double* darcyflux_;   // one flux per grid face
        const double* porevolume_;  // one volume per cell
        const double* source_;      // one source per cell
        double dt_;
        std::vector<double> saturation_;        // one per cell, only water saturation!
        std::vector<double> fractionalflow_;  // = m[0]/(m[0] + m[1]) per cell
        std::vector<int> reorder_iterations_;
        //std::vector<double> reorder_fval_;
        // For gravity segregation.
        std::vector<double> gravflux_;
        std::vector<double> mob_;
        std::vector<double> s0_;

        // Storing the upwind and downwind graphs for experiments.
        std::vector<int> ia_upw_;
        std::vector<int> ja_upw_;
        std::vector<int> ia_downw_;
        std::vector<int> ja_downw_;

        struct Residual;
        double fracFlow(double s, int cell) const;

        struct GravityResidual;
        void mobility(double s, int cell, double* mob) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTMODELTWOPHASE_HEADER_INCLUDED
