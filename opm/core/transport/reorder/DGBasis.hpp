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

#ifndef OPM_DGBASIS_HEADER_INCLUDED
#define OPM_DGBASIS_HEADER_INCLUDED

struct UnstructuredGrid;

namespace Opm
{

    /// Base class for Discontinuous Galerkin bases, intended for time-of-flight computations.
    class DGBasisInterface
    {
    public:
        /// Virtual destructor.
        virtual ~DGBasisInterface();

        /// The number of basis functions per cell.
        virtual int numBasisFunc() const = 0;

        /// The number of space dimensions.
        virtual int dimensions() const = 0;

        /// The polynomial degree of the basis functions.
        virtual int degree() const = 0;

        /// Evaluate all basis functions associated with cell at x,
        /// writing to f_x. The array f_x must have size equal to
        /// numBasisFunc().
        virtual void eval(const int cell,
                          const double* x,
                          double* f_x) const = 0;

        /// Evaluate gradients of all basis functions associated with
        /// cell at x, writing to grad_f_x. The array grad_f_x must
        /// have size numBasisFunc() * dimensions().  The dimensions()
        /// components of the first basis function gradient come
        /// before the components of the second etc.
        virtual void evalGrad(const int cell,
                              const double* x,
                              double* grad_f_x) const = 0;

        /// Modify basis coefficients to add to the function value.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to (f + increment) by modifying the c_i. This is done without
        /// modifying its gradient.
        /// \param[in]  increment     Add this value to the function.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void addConstant(const double increment,
                                 double* coefficients) const = 0;

        /// Modify basis coefficients to change the function's slope.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to a function g with the property that grad g = factor * grad f
        /// by modifying the c_i. This is done without modifying the average,
        /// i.e. the integrals of g and f over the cell are the same.
        /// \param[in]  factor        Multiply gradient by this factor.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void multiplyGradient(const double factor,
                                      double* coefficients) const = 0;

        /// Compute the average of the function f = sum_i c_i b_i.
        /// \param[in] coefficients  Coefficients {c_i} for a single cell.
        virtual double functionAverage(const double* coefficients) const = 0;
    };





    /// A class providing discontinuous Galerkin basis functions
    /// of bounded total degree.
    ///
    /// The basis functions are the following for each cell (example for 3d):
    ///     Degree 0: 1.
    ///     Degree 1: 1, x - xc, y - yc, z - zc
    /// where (xc, yc, zc) are the coordinates of the cell centroid.
    /// Further degrees await development.
    class DGBasisBoundedTotalDegree : public DGBasisInterface
    {
    public:
        /// Constructor.
        /// \param[in]  grid    grid on which basis is used (cell-wise)
        /// \param[in]  degree  polynomial degree of basis
        DGBasisBoundedTotalDegree(const UnstructuredGrid& grid, const int degree);

        /// Destructor.
        virtual ~DGBasisBoundedTotalDegree();

        /// The number of basis functions per cell.
        virtual int numBasisFunc() const;

        /// The number of space dimensions.
        virtual int dimensions() const;

        /// The polynomial degree of the basis functions.
        virtual int degree() const;

        /// Evaluate all basis functions associated with cell at x,
        /// writing to f_x. The array f_x must have size equal to
        /// numBasisFunc().
        virtual void eval(const int cell,
                          const double* x,
                          double* f_x) const;

        /// Evaluate gradients of all basis functions associated with
        /// cell at x, writing to grad_f_x. The array grad_f_x must
        /// have size numBasisFunc() * dimensions().  The dimensions()
        /// components of the first basis function gradient come
        /// before the components of the second etc.
        virtual void evalGrad(const int cell,
                              const double* x,
                              double* grad_f_x) const;

        /// Modify basis coefficients to add to the function value.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to (f + increment) by modifying the c_i. This is done without
        /// modifying its gradient.
        /// \param[in]  increment     Add this value to the function.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void addConstant(const double increment,
                                 double* coefficients) const;

        /// Modify basis coefficients to change the function's slope.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to a function g with the property that grad g = factor * grad f
        /// by modifying the c_i. This is done without modifying the average,
        /// i.e. the integrals of g and f over the cell are the same.
        /// \param[in]  factor        Multiply gradient by this factor.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void multiplyGradient(const double factor,
                                      double* coefficients) const;

        /// Compute the average of the function f = sum_i c_i b_i.
        /// \param[in] coefficients  Coefficients {c_i} for a single cell.
        virtual double functionAverage(const double* coefficients) const;

    private:
        const UnstructuredGrid& grid_;
        const int degree_;
    };




    /// A class providing discontinuous Galerkin basis functions of
    /// multi-degree 1 (bilinear or trilinear functions).
    ///
    /// The basis functions for a cell are the following
    ///     Degree 0: 1.
    /// (for 2 dims:)
    ///         (Bi)degree 1: (x-)(y-), (x-)(y+), (x+)(y-), (x+)(y+)
    ///         where (x-) = (1/2 - x + xc), (x+) = (1/2 + x - xc)
    ///         and xc is the x-coordinate of the cell centroid.
    ///         Similar for (y-), (y+).
    class DGBasisMultilin : public DGBasisInterface
    {
    public:
        /// Constructor.
        /// \param[in]  grid    grid on which basis is used (cell-wise)
        /// \param[in]  degree  polynomial degree of basis (in each coordinate)
        DGBasisMultilin(const UnstructuredGrid& grid, const int degree);

        /// Destructor.
        virtual ~DGBasisMultilin();

        /// The number of basis functions per cell.
        virtual int numBasisFunc() const;

        /// The number of space dimensions.
        virtual int dimensions() const;

        /// The polynomial degree of the basis functions.
        virtual int degree() const;

        /// Evaluate all basis functions associated with cell at x,
        /// writing to f_x. The array f_x must have size equal to
        /// numBasisFunc().
        virtual void eval(const int cell,
                          const double* x,
                          double* f_x) const;

        /// Evaluate gradients of all basis functions associated with
        /// cell at x, writing to grad_f_x. The array grad_f_x must
        /// have size numBasisFunc() * dimensions().  The dimensions()
        /// components of the first basis function gradient come
        /// before the components of the second etc.
        virtual void evalGrad(const int cell,
                              const double* x,
                              double* grad_f_x) const;

        /// Modify basis coefficients to add to the function value.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to (f + increment) by modifying the c_i. This is done without
        /// modifying its gradient.
        /// \param[in]  increment     Add this value to the function.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void addConstant(const double increment,
                                 double* coefficients) const;

        /// Modify basis coefficients to change the function's slope.
        /// A function f = sum_i c_i b_i is assumed, and we change
        /// it to a function g with the property that grad g = factor * grad f
        /// by modifying the c_i. This is done without modifying the average,
        /// i.e. the integrals of g and f over the cell are the same.
        /// \param[in]  factor        Multiply gradient by this factor.
        /// \param[out] coefficients  Coefficients {c_i} for a single cell.
        virtual void multiplyGradient(const double factor,
                                      double* coefficients) const;

        /// Compute the average of the function f = sum_i c_i b_i.
        /// \param[in] coefficients  Coefficients {c_i} for a single cell.
        virtual double functionAverage(const double* coefficients) const;

    private:
        const UnstructuredGrid& grid_;
        const int degree_;

    };




} // namespace Opm


#endif // OPM_DGBASIS_HEADER_INCLUDED
