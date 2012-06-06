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

#ifndef OPM_LINEARSOLVERINTERFACE_HEADER_INCLUDED
#define OPM_LINEARSOLVERINTERFACE_HEADER_INCLUDED


struct CSRMatrix;

namespace Opm
{


    /// Abstract interface for linear solvers.
    class LinearSolverInterface
    {
    public:
        /// Virtual destructor.
        virtual ~LinearSolverInterface();

        /// Struct for reporting data about the solution process back
        /// to the caller. The only field that is mandatory to set is
        /// 'converged' (even for direct solvers) to indicate success.
        struct LinearSolverReport
        {
            bool converged;
            int iterations;
            double residual_reduction;
        };

        /// Solve a linear system, with a matrix given in compressed sparse row format.
        /// \param[in] A           matrix in CSR format
        /// \param[in] rhs         array of length A->m containing the right hand side
        /// \param[inout] solution array of length A->m to which the solution will be written, may also be used
        ///                        as initial guess by iterative solvers.
        /// Note: this method is a convenience method that calls the virtual solve() method.
        LinearSolverReport solve(const CSRMatrix* A,
                                 const double* rhs,
                                 double* solution) const;

        /// Solve a linear system, with a matrix given in compressed sparse row format.
        /// \param[in] size        # of rows in matrix
        /// \param[in] nonzeros    # of nonzeros elements in matrix
        /// \param[in] ia          array of length (size + 1) containing start and end indices for each row
        /// \param[in] ja          array of length nonzeros containing column numbers for the nonzero elements
        /// \param[in] sa          array of length nonzeros containing the values of the nonzero elements
        /// \param[in] rhs         array of length size containing the right hand side
        /// \param[inout] solution array of length size to which the solution will be written, may also be used
        ///                        as initial guess by iterative solvers.
        virtual LinearSolverReport solve(const int size,
                                         const int nonzeros,
                                         const int* ia,
                                         const int* ja,
                                         const double* sa,
                                         const double* rhs,
                                         double* solution) const = 0;

        /// Set tolerance for the linear solver.
        /// \param[in] tol         tolerance value
        virtual void setTolerance(const double tol) = 0;

        /// Get tolerance for the linear solver.
        /// \param[out] tolerance value
        virtual double getTolerance() const = 0;
 
    };


} // namespace Opm



#endif // OPM_LINEARSOLVERINTERFACE_HEADER_INCLUDED
