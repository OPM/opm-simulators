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

#ifndef OPM_LINEARSOLVERISTL_HEADER_INCLUDED
#define OPM_LINEARSOLVERISTL_HEADER_INCLUDED


#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <string>


namespace Opm
{


    /// Abstract interface for linear solvers.
    class LinearSolverIstl : public LinearSolverInterface
    {
    public:
        /// Default constructor.
        /// All parameters controlling the solver are defaulted:
        ///   linsolver_residual_tolerance  1e-8
        ///   linsolver_verbosity           0
        ///   linsolver_type                1 ( = CG_AMG), alternatives are:
        ///                                 CG_ILU0 = 0, CG_AMG = 1, BiCGStab_ILU0 = 2
        ///   linsolver_save_system         false
        ///   linsolver_save_filename       <empty string>
        ///   linsolver_max_iterations      0 (unlimited)
        LinearSolverIstl();

        /// Construct from parameters
        /// Accepted parameters are, with defaults, listed in the
        /// default constructor.
        LinearSolverIstl(const parameter::ParameterGroup& param);

        /// Destructor.
        virtual ~LinearSolverIstl();

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
                                         double* solution) const;
    private:
        double linsolver_residual_tolerance_;
        int linsolver_verbosity_;
        enum LinsolverType { CG_ILU0 = 0, CG_AMG = 1, BiCGStab_ILU0 = 2 };
        LinsolverType linsolver_type_;
        bool linsolver_save_system_;
        std::string linsolver_save_filename_;
        int linsolver_max_iterations_;
    };


} // namespace Opm



#endif // OPM_LINEARSOLVERISTL_HEADER_INCLUDED
