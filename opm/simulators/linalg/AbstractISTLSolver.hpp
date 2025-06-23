/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_ABSTRACTISTLSOLVER_HEADER_INCLUDED
#define OPM_ABSTRACTISTLSOLVER_HEADER_INCLUDED

#include <opm/common/Exceptions.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm
{

/**
 * \brief Abstract interface for ISTL solvers.
 *
 * This class defines the interface for ISTL solvers used in OPM.
 * It provides methods for preparing the solver, setting and getting
 * residuals, solving the system, and managing communication.
 *
 * \note This class is used in the ISTLSolverRuntimeOptionProxy which we
 *       where we can set the solver type at runtime, and this proxy holds a
 *       pointer to an instance of this class.
 *
 * \todo We should remove the use of setResidual() and setMatrix() methods
 *       and instead use prepare() with the SparseMatrixAdapter and Vector
 *       directly. This would simplify the interface and reduce the number of
 *       required method calls before solving the system.
 */
template <class TypeTag>
class AbstractISTLSolver
{
public:
#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;

    virtual ~AbstractISTLSolver() = default;

    /**
     * \brief Signals that the memory for the matrix internally in the solver could be erased
     *
     * \note This call could be ignored by the solver, but it is a hint that the
     *       solver does not need the matrix anymore.
     */
    virtual void eraseMatrix() = 0;

    /**
     * \brief Set the active solver by its index.
     *
     * \param num The index of the solver to set as active.
     *
     * \note The index corresponds to the order in which solvers are registered.
     */
    virtual void setActiveSolver(int num) = 0;

    /**
     * \brief Get the number of available solvers.
     *
     * \return The number of solvers that can be used.
     */
    virtual int numAvailableSolvers() const = 0;

    /**
     * \brief Prepare the solver with the given matrix and right-hand side vector.
     *
     * This method initializes the solver with the provided matrix and vector,
     * preparing it for solving the system of equations.
     *
     * \param M The matrix representing the system of equations.
     * \param b The right-hand side vector.
     */
    virtual void prepare(const Matrix& M, Vector& b) = 0;

    /**
     * \brief Prepare the solver with the given sparse matrix and right-hand side vector.
     *
     * This method initializes the solver with the provided sparse matrix and vector,
     * preparing it for solving the system of equations.
     *
     * \param M The sparse matrix representing the system of equations.
     * \param b The right-hand side vector.
     *
     * \note This method should be called *in addition* to setResidual() and setMatrix() before calling solve().
     */
    virtual void prepare(const SparseMatrixAdapter& M, Vector& b) = 0;

    /**
     * \brief Set the residual vector.
     *
     * This method sets the residual vector for the solver.
     *
     * \param b The residual vector to set.
     *
     * \note This method should be called *in addition* to prepare() and setMatrix() before calling solve().
     */
    virtual void setResidual(Vector& b) = 0;

    /**
     * \brief Get the residual vector.
     *
     * This method retrieves the current residual vector from the solver.
     *
     * \param b The vector to store the residual.
     */
    virtual void getResidual(Vector& b) const = 0;

    /**
     * \brief Set the matrix for the solver.
     *
     * This method sets the matrix that the solver will use to solve the system of equations.
     *
     * \param M The sparse matrix adapter containing the matrix data.
     *
     * \note This method should be called *in addition* to prepare() and setResidual() before calling solve().
     */
    virtual void setMatrix(const SparseMatrixAdapter& M) = 0;

    /**
     * \brief Solve the system of equations Ax = b.
     *
     * This method solves the linear system represented by the matrix A and the right-hand side vector b,
     * storing the solution in vector x.
     *
     * \param x The vector to store the solution.
     * \return true if the solver converged, false otherwise.
     *
     * Before this function is called, the following function calls should have been made:
     * - prepare(const Matrix& M, Vector& b) or prepare(const SparseMatrixAdapter& M, Vector& b)
     * - setResidual(Vector& b) or setResidual(const Vector& b)
     * - setMatrix(const SparseMatrixAdapter& M)
     */
    virtual bool solve(Vector& x) = 0;

    /**
     * \brief Get the number of iterations used in the last solve.
     *
     * This method returns the number of iterations that the solver performed during the last call to solve().
     *
     * \return The number of iterations.
     *
     * \note This value is only valid after a call to solve().
     */
    virtual int iterations() const = 0;

    /**
     * \brief Get the communication object used by the solver.
     *
     * This method returns a pointer to the communication object used by the solver.
     *
     * \return A pointer to the communication object.
     */
    virtual const CommunicationType* comm() const = 0;

    /**
     * \brief Get the count of how many times the solver has been called.
     *
     * This method returns the number of times the solve() method has been called.
     *
     * \return The count of solve calls.
     */
    virtual int getSolveCount() const = 0;

protected:

    /**
     * \brief Check the convergence of the linear solver.
     *
     * This method checks if the linear solver has converged based on the result and parameters.
     *
     * \param result The result of the linear solver.
     * \param parameters The parameters used for the linear solver.
     * \return true if the solver has converged, false otherwise.
     */
    static bool checkConvergence(const Dune::InverseOperatorResult& result,
                                 const FlowLinearSolverParameters& parameters)
    {
        if (!result.converged) {
            if (result.reduction < parameters.relaxed_linear_solver_reduction_) {
                std::stringstream ss;
                ss << "Full linear solver tolerance not achieved. The reduction is:" << result.reduction << " after "
                   << result.iterations << " iterations ";
                OpmLog::warning(ss.str());
                return true;
            }
        }
        // Check for failure of linear solver.
        if (!parameters.ignoreConvergenceFailure_ && !result.converged) {
            const std::string msg("Convergence failure for linear solver.");
            OPM_THROW_NOLOG(NumericalProblem, msg);
        }

        return result.converged;
    }
};
} // namespace Opm

#endif
