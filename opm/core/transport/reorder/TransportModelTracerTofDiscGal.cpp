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

#include <opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Opm
{


    /// A class providing discontinuous Galerkin basis functions.
    struct DGBasis
    {
        static int numBasisFunc(const int dimensions,
                                const int degree)
        {
            switch (dimensions) {
            case 1:
                return degree + 1;
            case 2:
                return (degree + 2)*(degree + 1)/2;
            case 3:
                return (degree + 3)*(degree + 2)*(degree + 1)/6;
            default:
                THROW("Dimensions must be 1, 2 or 3.");
            }
        }

        /// Evaluate all nonzero basis functions at x,
        /// writing to f_x. The array f_x must have
        /// size numBasisFunc(grid.dimensions, degree).
        ///
        /// The basis functions are the following
        ///     Degree 0: 1.
        ///     Degree 1: x - xc, y - yc, z - zc etc.
        /// Further degrees await development.
        static void eval(const UnstructuredGrid& grid,
                         const int cell,
                         const int degree,
                         const double* x,
                         double* f_x)
        {
            const int dim = grid.dimensions;
            const double* cc = grid.cell_centroids + dim*cell;
            // Note intentional fallthrough in this switch statement!
            switch (degree) {
            case 1:
                for (int ix = 0; ix < dim; ++ix) {
                    f_x[1 + ix] = x[ix] - cc[ix];
                }
            case 0:
                f_x[0] = 1;
                break;
            default:
                THROW("Maximum degree is 1 for now.");
            }
        }

        /// Evaluate gradients of all nonzero basis functions at x,
        /// writing to grad_f_x. The array grad_f_x must have size
        /// numBasisFunc(grid.dimensions, degree) * grid.dimensions.
        /// The <grid.dimensions> components of the first basis function
        /// gradient come before the components of the second etc.
        static void evalGrad(const UnstructuredGrid& grid,
                             const int /*cell*/,
                             const int degree,
                             const double* /*x*/,
                             double* grad_f_x)
        {
            const int dim = grid.dimensions;
            const int num_basis = numBasisFunc(dim, degree);
            std::fill(grad_f_x, grad_f_x + num_basis*dim, 0.0);
            if (degree > 1) {
                THROW("Maximum degree is 1 for now.");
            } else if (degree == 1) {
                for (int ix = 0; ix < dim; ++ix) {
                    grad_f_x[dim*(ix + 1) + ix] = 1.0;
                }
            }
        }
    };




    /// A class providing numerical quadrature functions.
    struct Quadrature
    {
        static void x(){}
    };





    /// Construct solver.
    /// \param[in] grid      A 2d or 3d grid.
    TransportModelTracerTofDiscGal::TransportModelTracerTofDiscGal(const UnstructuredGrid& grid)
        : grid_(grid)
    {
    }




    /// Solve for time-of-flight at next timestep.
    /// \param[in]  darcyflux         Array of signed face fluxes.
    /// \param[in]  porevolume        Array of pore volumes.
    /// \param[in]  source            Transport source term.
    /// \param[in]  degree            Polynomial degree of DG basis functions used.
    /// \param[out] tof_coeff         Array of time-of-flight solution coefficients.
    ///                               The values are ordered by cell, meaning that
    ///                               the K coefficients corresponding to the first
    ///                               cell comes before the K coefficients corresponding
    ///                               to the second cell etc.
    ///                               K depends on degree and grid dimension.
    void TransportModelTracerTofDiscGal::solveTof(const double* darcyflux,
                                                  const double* porevolume,
                                                  const double* source,
                                                  const int degree,
                                                  std::vector<double>& tof_coeff)
    {
        darcyflux_ = darcyflux;
        porevolume_ = porevolume;
        source_ = source;
#ifndef NDEBUG
        // Sanity check for sources.
        const double cum_src = std::accumulate(source, source + grid_.number_of_cells, 0.0);
        if (std::fabs(cum_src) > *std::max_element(source, source + grid_.number_of_cells)*1e-2) {
            THROW("Sources do not sum to zero: " << cum_src);
        }
#endif
        degree_ = degree;
        tof_coeff.resize(grid_.number_of_cells);
        std::fill(tof_coeff.begin(), tof_coeff.end(), 0.0);
        tof_coeff_ = &tof_coeff[0];
        reorderAndTransport(grid_, darcyflux);
    }




    void TransportModelTracerTofDiscGal::solveSingleCell(const int cell)
    {
        // Residual:
        // For each cell K, basis function b_j (spanning V_h),
        // writing the solution u_h|K = \sum_i c_i b_i
        //  Res = - \int_K \sum_i c_i b_i v(x) \cdot \grad b_j dx
        //        + \int_{\partial K} F(u_h, u_h^{ext}, v(x) \cdot n) b_j ds
        //        - \int_K \phi b_j
        THROW("Not implemented yet!");
    }




    void TransportModelTracerTofDiscGal::solveMultiCell(const int num_cells, const int* cells)
    {
        std::cout << "Pretending to solve multi-cell dependent equation with " << num_cells << " cells." << std::endl;
        for (int i = 0; i < num_cells; ++i) {
            solveSingleCell(cells[i]);
        }
    }




} // namespace Opm
