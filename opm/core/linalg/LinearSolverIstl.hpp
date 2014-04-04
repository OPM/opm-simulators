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


    /// Concrete class encapsulating some dune-istl linear solvers.
    class LinearSolverIstl : public LinearSolverInterface
    {
    public:
        /// Default constructor.
        /// All parameters controlling the solver are defaulted:
        ///   linsolver_residual_tolerance  1e-8
        ///   linsolver_verbosity           0
        ///   linsolver_type                1 ( = CG_AMG), alternatives are:
        ///                                 CG_ILU0 = 0, CG_AMG = 1, BiCGStab_ILU0 = 2
        ///                                 FastAMG=3, KAMG=4 };
        ///   linsolver_save_system         false
        ///   linsolver_save_filename       <empty string>
        ///   linsolver_max_iterations      0 (unlimited=5000)
        ///   linsolver_residual_tolerance  1e-8
        ///   linsolver_smooth_steps        2
        ///   linsolver_prolongate_factor   1.6
        ///   linsolver_verbosity           0
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

        /// Set tolerance for the residual in dune istl linear solver.
        /// \param[in] tol         tolerance value
        virtual void setTolerance(const double tol);

        /// Get tolerance ofthe linear solver.
        /// \param[out] tolerance value
        virtual double getTolerance() const;

    private:
        /// \brief Solve the linear system using ISTL
        /// \param[in] opA The linear operator of the system to solve.
        /// \param[out]    solution C array for storing the solution vector.
        /// \param[in]     rhs C array containing the right hand side.
        /// \param[in]     sp The scalar product to use.
        /// \param[in]     comm The information about the parallel domain decomposition.
        /// \param[in]     maxit The maximum number of iterations allowed.
        template<class O, class S, class C>
        LinearSolverReport solveSystem(O& opA, double* solution, const double *rhs,
                                       S& sp, const C& comm, int maxit) const;

        double linsolver_residual_tolerance_;
        int linsolver_verbosity_;
        enum LinsolverType { CG_ILU0 = 0, CG_AMG = 1, BiCGStab_ILU0 = 2, FastAMG=3, KAMG=4 };
        LinsolverType linsolver_type_;
        bool linsolver_save_system_;
        std::string linsolver_save_filename_;
        int linsolver_max_iterations_;
        /** \brief The number smoothing steps to apply in AMG. */
        int linsolver_smooth_steps_;
        /** \brief The factor to scale the coarse grid correction with. */
        double linsolver_prolongate_factor_;

    };


} // namespace Opm



#endif // OPM_LINEARSOLVERISTL_HEADER_INCLUDED
