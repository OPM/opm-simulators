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

#include "linalgproperties.hh"
#include "parallelbasebackend.hh"
#include "istlsolverwrappers.hh"
#include "istlsparsematrixadapter.hh"

#include <dune/common/version.hh>

namespace Opm::Properties::TTag {

// Create new type tag
struct ParallelIstlLinearSolver { using InheritsFrom = std::tuple<ParallelBaseLinearSolver>; };

} // namespace Opm::Properties::TTag

namespace Opm::Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides all unmodified linear solvers from dune-istl
 *
 * To set the linear solver, use
 * \code
 * template<class TypeTag>
 * struct LinearSolverWrapper<TypeTag, TTag::YourTypeTag>
 * { using type = Opm::Linear::SolverWrapper$SOLVER<TypeTag>; };
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
 * template<class TypeTag>
 * struct PreconditionerWrapper<TypeTag, TTag::YourTypeTag>
 * { using type = Opm::Linear::PreconditionerWrapper$PRECONDITIONER<TypeTag>; };
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
    using ParentType = ParallelBaseBackend<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using LinearSolverWrapper = GetPropType<TypeTag, Properties::LinearSolverWrapper>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

    using ParallelOperator = typename ParentType::ParallelOperator;
    using OverlappingVector = typename ParentType::OverlappingVector;
    using ParallelPreconditioner = typename ParentType::ParallelPreconditioner;
    using ParallelScalarProduct = typename ParentType::ParallelScalarProduct;

    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;
    using RawLinearSolver = typename LinearSolverWrapper::RawSolver;

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

} // namespace Opm::Linear

namespace Opm::Properties {

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::ParallelIstlLinearSolver>
{ using type = Opm::Linear::ParallelIstlSolverBackend<TypeTag>; };

template<class TypeTag>
struct LinearSolverWrapper<TypeTag, TTag::ParallelIstlLinearSolver>
{ using type = Opm::Linear::SolverWrapperBiCGStab<TypeTag>; };

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
template<class TypeTag>
struct PreconditionerWrapper<TypeTag, TTag::ParallelIstlLinearSolver>
{ using type = Opm::Linear::PreconditionerWrapperILU<TypeTag>; };
#else
template<class TypeTag>
struct PreconditionerWrapper<TypeTag, TTag::ParallelIstlLinearSolver>
{ using type = Opm::Linear::PreconditionerWrapperILU0<TypeTag>; };
#endif

//! set the GMRes restart parameter to 10 by default
template<class TypeTag>
struct GMResRestart<TypeTag, TTag::ParallelIstlLinearSolver> { static constexpr int value = 10; };

} // namespace Opm::Properties

#endif
