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

#ifndef OPM_LINEARSOLVERPETSC_HEADER_INCLUDED
#define OPM_LINEARSOLVERPETSC_HEADER_INCLUDED
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <string>

namespace Opm
{


    /// Concrete class encapsulating some Petsc linear solvers.
    class LinearSolverPetsc : public LinearSolverInterface
    {
    public:
        /// Default constructor.
        /// All parameters controlling the solver are defaulted:
        ///   ksp_type                      gmres
        ///   pc_type                       jacobi
        ///   ksp_view                      0
        ///   ksp_rtol                      1e-5
        ///   ksp_atol                      1e-50
        ///   ksp_dtol                      1e5
        ///   ksp_max_it                    1e5
        LinearSolverPetsc();

        /// Construct from parameters
        /// Accepted parameters are, with defaults, listed in the
        /// default constructor.
        LinearSolverPetsc(const parameter::ParameterGroup& param);

        /// Destructor.
        virtual ~LinearSolverPetsc();

        using LinearSolverInterface::solve;

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
                                         double* solution,
                                         const boost::any&) const;

        /// Set tolerance for the residual in dune istl linear solver.
        /// \param[in] tol         tolerance value
        virtual void setTolerance(const double tol);

        /// Get tolerance ofthe linear solver.
        /// \param[out] tolerance value
        virtual double getTolerance() const;
    private:
        std::string     ksp_type_;
        std::string     pc_type_;
        int             ksp_view_;
        double          rtol_;
        double          atol_;
        double          dtol_;
        int             maxits_;
    };


} // namespace Opm



#endif // OPM_LINEARSOLVERPETSC_HEADER_INCLUDED
