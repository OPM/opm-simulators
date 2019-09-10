// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::ParallelIstlSolverBackend
 */
#ifndef EWOMS_PARALLEL_ISTL_BACKEND_HH
#define EWOMS_PARALLEL_ISTL_BACKEND_HH

#include "parallelbasebackend.hh"
#include "istlsolverwrappers.hh"
#include "istlsparsematrixadapter.hh"

#include <dune/common/version.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(ParallelIstlLinearSolver, INHERITS_FROM(ParallelBaseLinearSolver));

NEW_PROP_TAG(LinearSolverWrapper);
NEW_PROP_TAG(SparseMatrixAdapter);

//! number of iterations between solver restarts for the GMRES solver
NEW_PROP_TAG(GMResRestart);

END_PROPERTIES

namespace Opm {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides all unmodified linear solvers from dune-istl
 *
 * To set the linear solver, use
 * \code
 * SET_TYPE_PROP(YourTypeTag, LinearSolverWrapper,
 *               Opm::Linear::SolverWrapper$SOLVER<TypeTag>);
 * \endcode
 *
 * The possible choices for '\c $SOLVER' are:
 * - \c Richardson: A fixpoint solver using the Richardson iteration
 * - \c SteepestDescent: The steepest descent solver
 * - \c ConjugatedGradients: A conjugated gradients solver
 * - \c BiCGStab: A stabilized bi-conjugated gradients solver
 * - \c MinRes: A solver based on the  minimized residual algorithm
 * - \c RestartedGMRes: A restarted GMRES solver
 *
 * Chosing the preconditioner works in an analogous way:
 * \code
 * SET_TYPE_PROP(YourTypeTag, PreconditionerWrapper,
 *               Opm::Linear::PreconditionerWrapper$PRECONDITIONER<TypeTag>);
 * \endcode
 *
 * Where the choices possible for '\c $PRECONDITIONER' are:
 * - \c Jacobi: A Jacobi preconditioner
 * - \c GaussSeidel: A Gauss-Seidel preconditioner
 * - \c SSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - \c SOR: A successive overrelaxation (SOR) preconditioner
 * - \c ILUn: An ILU(n) preconditioner
 * - \c ILU0: A specialized (and optimized) ILU(0) preconditioner
 */
template <class TypeTag>
class ParallelIstlSolverBackend : public ParallelBaseBackend<TypeTag>
{
    typedef ParallelBaseBackend<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverWrapper) LinearSolverWrapper;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

    typedef typename ParentType::ParallelOperator ParallelOperator;
    typedef typename ParentType::OverlappingVector OverlappingVector;
    typedef typename ParentType::ParallelPreconditioner ParallelPreconditioner;
    typedef typename ParentType::ParallelScalarProduct ParallelScalarProduct;

    typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlock;
    typedef typename LinearSolverWrapper::RawSolver RawLinearSolver;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock> >::value,
                  "The ParallelIstlSolverBackend linear solver backend requires the IstlSparseMatrixAdapter");

public:
    ParallelIstlSolverBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        LinearSolverWrapper::registerParameters();
    }

protected:
    friend ParentType;

    std::shared_ptr<RawLinearSolver> prepareSolver_(ParallelOperator& parOperator,
                                                    ParallelScalarProduct& parScalarProduct,
                                                    ParallelPreconditioner& parPreCond)
    {
        return solverWrapper_.get(parOperator,
                                  parScalarProduct,
                                  parPreCond);
    }

    void cleanupSolver_()
    {
        solverWrapper_.cleanup();
    }

    std::pair<bool, int> runSolver_(std::shared_ptr<RawLinearSolver> solver)
    {
        Dune::InverseOperatorResult result;
        solver->apply(*this->overlappingx_, *this->overlappingb_, result);
        return std::make_pair(result.converged, result.iterations);
    }

    LinearSolverWrapper solverWrapper_;
};

}} // namespace Linear, Ewoms

BEGIN_PROPERTIES

SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverBackend,
              Opm::Linear::ParallelIstlSolverBackend<TypeTag>);

SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverWrapper,
              Opm::Linear::SolverWrapperBiCGStab<TypeTag>);

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
SET_TYPE_PROP(ParallelIstlLinearSolver,
              PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILU<TypeTag>);
#else
SET_TYPE_PROP(ParallelIstlLinearSolver,
              PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILU0<TypeTag>);
#endif

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(ParallelIstlLinearSolver, GMResRestart, 10);

END_PROPERTIES

#endif
