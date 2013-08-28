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

#include "config.h"
#include <opm/core/tof/DGBasis.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <numeric>

namespace Opm
{

    // ----------------  Methods for class DGBasisInterface ----------------


    /// Virtual destructor.
    DGBasisInterface::~DGBasisInterface()
    {
    }

    /// Evaluate function f = sum_i c_i b_i at the point x.
    /// Note that this function is not virtual, but implemented in
    /// terms of the virtual functions of the class.
    /// \param[in] cell          Cell index
    /// \param[in] coefficients  Coefficients {c_i} for a single cell.
    /// \param[in] x             Point at which to compute f(x).
    double DGBasisInterface::evalFunc(const int cell,
                                      const double* coefficients,
                                      const double* x) const
    {
        bvals_.resize(numBasisFunc());
        eval(cell, x, &bvals_[0]);
        return std::inner_product(bvals_.begin(), bvals_.end(), coefficients, 0.0);
    }



    // ----------------  Methods for class DGBasisBoundedTotalDegree ----------------


    /// Constructor.
    /// \param[in]  grid    grid on which basis is used (cell-wise)
    /// \param[in]  degree  polynomial degree of basis
    DGBasisBoundedTotalDegree::DGBasisBoundedTotalDegree(const UnstructuredGrid& grid,
                                                         const int degree_arg)
        : grid_(grid),
          degree_(degree_arg)
    {
        if (grid_.dimensions > 3) {
            OPM_THROW(std::runtime_error, "Grid dimension must be 1, 2 or 3.");
        }
        if (degree_ > 1 || degree_ < 0) {
            OPM_THROW(std::runtime_error, "Degree must be 0 or 1.");
        }
    }

    /// Destructor.
    DGBasisBoundedTotalDegree::~DGBasisBoundedTotalDegree()
    {
    }

    /// The number of basis functions per cell.
    int DGBasisBoundedTotalDegree::numBasisFunc() const
    {
        switch (dimensions()) {
        case 1:
            return degree_ + 1;
        case 2:
            return (degree_ + 2)*(degree_ + 1)/2;
        case 3:
            return (degree_ + 3)*(degree_ + 2)*(degree_ + 1)/6;
        default:
            OPM_THROW(std::runtime_error, "Dimensions must be 1, 2 or 3.");
        }
    }

    /// The number of space dimensions.
    int DGBasisBoundedTotalDegree::dimensions() const
    {
        return grid_.dimensions;
    }

    /// The polynomial degree of the basis functions.
    int DGBasisBoundedTotalDegree::degree() const
    {
        return degree_;
    }

    /// Evaluate all basis functions associated with cell at x,
    /// writing to f_x. The array f_x must have size equal to
    /// numBasisFunc().
    void DGBasisBoundedTotalDegree::eval(const int cell,
                                         const double* x,
                                         double* f_x) const
    {
        const int dim = dimensions();
        const double* cc = grid_.cell_centroids + dim*cell;
        // Note intentional fallthrough in this switch statement!
        switch (degree_) {
        case 1:
            for (int ix = 0; ix < dim; ++ix) {
                f_x[1 + ix] = x[ix] - cc[ix];
            }
        case 0:
            f_x[0] = 1;
            break;
        default:
            OPM_THROW(std::runtime_error, "Maximum degree is 1 for now.");
        }
    }


    /// Evaluate gradients of all basis functions associated with
    /// cell at x, writing to grad_f_x. The array grad_f_x must
    /// have size numBasisFunc() * dimensions().  The dimensions()
    /// components of the first basis function gradient come
    /// before the components of the second etc.
    void DGBasisBoundedTotalDegree::evalGrad(const int /*cell*/,
                                             const double* /*x*/,
                                             double* grad_f_x) const
    {
        const int dim = dimensions();
        const int num_basis = numBasisFunc();
        std::fill(grad_f_x, grad_f_x + num_basis*dim, 0.0);
        if (degree_ == 1) {
            for (int ix = 0; ix < dim; ++ix) {
                grad_f_x[dim*(ix + 1) + ix] = 1.0;
            }
        }
    }

    /// Modify basis coefficients to add to the function value.
    /// A function f = sum_i c_i b_i is assumed, and we change
    /// it to (f + increment) by modifying the c_i. This is done without
    /// modifying its gradient.
    /// \param[in]  increment     Add this value to the function.
    /// \param[out] coefficients  Coefficients {c_i} for a single cell.
    void DGBasisBoundedTotalDegree::addConstant(const double increment,
                                                double* coefficients) const
    {
        coefficients[0] += increment;
    }

    /// Modify basis coefficients to change the function's slope.
    /// A function f = sum_i c_i b_i is assumed, and we change
    /// it to a function g with the property that grad g = factor * grad f
    /// by modifying the c_i. This is done without modifying the average,
    /// i.e. the integrals of g and f over the cell are the same.
    /// \param[in]  factor        Multiply gradient by this factor.
    /// \param[out] coefficients  Coefficients {c_i} for a single cell.
    void DGBasisBoundedTotalDegree::multiplyGradient(const double factor,
                                                     double* coefficients) const
    {
        const int nb = numBasisFunc();
        for (int ix = 1; ix < nb; ++ix) {
            coefficients[ix] *= factor;
        }
    }

    /// Compute the average of the function f = sum_i c_i b_i.
    /// \param[in] coefficients  Coefficients {c_i} for a single cell.
    double DGBasisBoundedTotalDegree::functionAverage(const double* coefficients) const
    {
        return coefficients[0];
    }



    // ----------------  Methods for class DGBasisMultilin ----------------


    /// Constructor.
    /// \param[in]  grid    grid on which basis is used (cell-wise)
    /// \param[in]  degree  polynomial degree of basis
    DGBasisMultilin::DGBasisMultilin(const UnstructuredGrid& grid,
                                     const int degree_arg)
        : grid_(grid),
          degree_(degree_arg)
    {
        if (grid_.dimensions > 3) {
            OPM_THROW(std::runtime_error, "Grid dimension must be 1, 2 or 3.");
        }
        if (degree_ > 1 || degree_ < 0) {
            OPM_THROW(std::runtime_error, "Degree must be 0 or 1.");
        }
    }

    /// Destructor.
    DGBasisMultilin::~DGBasisMultilin()
    {
    }

    /// The number of basis functions per cell.
    int DGBasisMultilin::numBasisFunc() const
    {
        switch (dimensions()) {
        case 1:
            return degree_ + 1;
        case 2:
            return (degree_ + 1)*(degree_ + 1);
        case 3:
            return (degree_ + 1)*(degree_ + 1)*(degree_ + 1);
        default:
            OPM_THROW(std::runtime_error, "Dimensions must be 1, 2 or 3.");
        }
    }

    /// The number of space dimensions.
    int DGBasisMultilin::dimensions() const
    {
        return grid_.dimensions;
    }

    /// The polynomial degree of the basis functions.
    int DGBasisMultilin::degree() const
    {
        return degree_;
    }

    /// Evaluate all basis functions associated with cell at x,
    /// writing to f_x. The array f_x must have size equal to
    /// numBasisFunc().
    void DGBasisMultilin::eval(const int cell,
                               const double* x,
                               double* f_x) const
    {
        const int dim = dimensions();
        const int num_basis = numBasisFunc();
        const double* cc = grid_.cell_centroids + dim*cell;
        switch (degree_) {
        case 0:
            f_x[0] = 1;
            break;
        case 1:
            std::fill(f_x, f_x + num_basis, 1.0);
            for (int dd = 0; dd < dim; ++dd) {
                const double f[2] = { 0.5 - x[dd] + cc[dd],  0.5 + x[dd] - cc[dd] };
                const int divi = 1 << (dim - dd - 1); // { 4, 2, 1 } for 3d, for example.
                for (int ix = 0; ix < num_basis; ++ix) {
                    f_x[ix] *= f[(ix/divi) % 2];
                }
            }
            break;
        default:
            OPM_THROW(std::runtime_error, "Maximum degree is 1 for now.");
        }
    }


    /// Evaluate gradients of all basis functions associated with
    /// cell at x, writing to grad_f_x. The array grad_f_x must
    /// have size numBasisFunc() * dimensions().  The dimensions()
    /// components of the first basis function gradient come
    /// before the components of the second etc.
    void DGBasisMultilin::evalGrad(const int cell,
                                   const double* x,
                                   double* grad_f_x) const
    {
        const int dim = dimensions();
        const int num_basis = numBasisFunc();
            const double* cc = grid_.cell_centroids + dim*cell;
            switch (degree_) {
            case 0:
                std::fill(grad_f_x, grad_f_x + num_basis*dim, 0.0);
                break;
            case 1:
                std::fill(grad_f_x, grad_f_x + num_basis*dim, 1.0);
                for (int dd = 0; dd < dim; ++dd) {
                    const double f[2] = { 0.5 - x[dd] + cc[dd],  0.5 + x[dd] - cc[dd] };
                    const double fder[2] = { -1.0,  1.0 };
                    const int divi = 1 << (dim - dd - 1); // { 4, 2, 1 } for 3d, for example.
                    for (int ix = 0; ix < num_basis; ++ix) {
                        const int ind = (ix/divi) % 2;
                        for (int dder = 0; dder < dim; ++dder) {
                            grad_f_x[ix*dim + dder] *= (dder == dd ? fder[ind] : f[ind]);
                        }
                    }
                }
                break;
            default:
                OPM_THROW(std::runtime_error, "Maximum degree is 1 for now.");
            }
    }

    /// Modify basis coefficients to add to the function value.
    /// A function f = sum_i c_i b_i is assumed, and we change
    /// it to (f + increment) by modifying the c_i. This is done without
    /// modifying its gradient.
    /// \param[in]  increment     Add this value to the function.
    /// \param[out] coefficients  Coefficients {c_i} for a single cell.
    void DGBasisMultilin::addConstant(const double increment,
                                      double* coefficients) const
    {
        const int nb = numBasisFunc();
        const double term = increment/double(nb);
        for (int ix = 0; ix < nb; ++ix) {
            coefficients[ix] += term;
        }
    }

    /// Modify basis coefficients to change the function's slope.
    /// A function f = sum_i c_i b_i is assumed, and we change
    /// it to a function g with the property that grad g = factor * grad f
    /// by modifying the c_i. This is done without modifying the average,
    /// i.e. the integrals of g and f over the cell are the same.
    /// \param[in]  factor        Multiply gradient by this factor.
    /// \param[out] coefficients  Coefficients {c_i} for a single cell.
    void DGBasisMultilin::multiplyGradient(const double factor,
                                           double* coefficients) const
    {
        const int nb = numBasisFunc();
        const double aver = functionAverage(coefficients);
        for (int ix = 0; ix < nb; ++ix) {
            coefficients[ix] = factor*(coefficients[ix] - aver) + aver;
        }
    }

    /// Compute the average of the function f = sum_i c_i b_i.
    /// \param[in] coefficients  Coefficients {c_i} for a single cell.
    double DGBasisMultilin::functionAverage(const double* coefficients) const
    {
        const int nb = numBasisFunc();
        return std::accumulate(coefficients, coefficients + nb, 0.0)/double(nb);
    }

} // namespace Opm
