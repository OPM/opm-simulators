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

#ifndef OPM_LINEARSOLVERFACTORY_HEADER_INCLUDED
#define OPM_LINEARSOLVERFACTORY_HEADER_INCLUDED


#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <boost/shared_ptr.hpp>

namespace Opm
{

    namespace parameter { class ParameterGroup; }


    /// Concrete class encapsulating any available linear solver.
    /// For the moment, this means UMFPACK and dune-istl.
    /// Since both are optional dependencies, either or both
    /// may be unavailable, depending on configuration.
    class LinearSolverFactory : public LinearSolverInterface
    {
    public:
        /// Default constructor.
        LinearSolverFactory();

        /// Construct from parameters.
        /// The accepted parameters are (default) (allowed values):
        ///    linsolver ("umfpack")   ("umfpack", "istl")
        /// For the umfpack solver to be available, this class must be
        /// compiled with UMFPACK support, as indicated by the
        /// variable HAVE_SUITESPARSE_UMFPACK_H in config.h.
        /// For the istl solver to be available, this class must be
        /// compiled with dune-istl support, as indicated by the
        /// variable HAVE_DUNE_ISTL in config.h.
        /// Any further parameters are passed on to the constructors
        /// of the actual solver used, see LinearSolverUmfpack
        /// and LinearSolverIstl for details.
        LinearSolverFactory(const parameter::ParameterGroup& param);

        /// Destructor.
        virtual ~LinearSolverFactory();

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

        /// Set tolerance for the linear solver.
        /// \param[in] tol         tolerance value
        /// Not used for LinearSolverFactory
        virtual void setTolerance(const double tol);

        /// Get tolerance for the linear solver.
        /// \param[out] tolerance value
        /// Not used for LinearSolverFactory. Returns -1.
        virtual double getTolerance() const;

    private:
        boost::shared_ptr<LinearSolverInterface> solver_;
    };


} // namespace Opm


#endif // OPM_LINEARSOLVERFACTORY_HEADER_INCLUDED
