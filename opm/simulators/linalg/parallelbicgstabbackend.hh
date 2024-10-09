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
 * \copydoc Opm::Linear::ParallelBiCGStabSolverBackend
 */
#ifndef EWOMS_PARALLEL_BICGSTAB_BACKEND_HH
#define EWOMS_PARALLEL_BICGSTAB_BACKEND_HH

#include <opm/simulators/linalg/bicgstabsolver.hh>
#include <opm/simulators/linalg/combinedcriterion.hh>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/linalgparameters.hh>
#include <opm/simulators/linalg/linalgproperties.hh>
#include <opm/simulators/linalg/parallelbasebackend.hh>

#include <memory>

namespace Opm::Linear {

template <class TypeTag>
class ParallelBiCGStabSolverBackend;

} // namespace Opm::Linear

namespace Opm::Properties {

// Create new type tags
namespace TTag {

struct ParallelBiCGStabLinearSolver
{ using InheritsFrom = std::tuple<ParallelBaseLinearSolver>; };

} // end namespace TTag

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::ParallelBiCGStabLinearSolver>
{ using type = Opm::Linear::ParallelBiCGStabSolverBackend<TypeTag>; };

} // namespace Opm::Properties

namespace Opm::Linear {

/*!
 * \ingroup Linear
 *
 * \brief Implements a generic linear solver abstraction.
 *
 * Chosing the preconditioner works by setting the "PreconditionerWrapper" property:
 *
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
 * - \c ILU0: An ILU(0) preconditioner. The results of this
 *            preconditioner are the same as setting the
 *            PreconditionerOrder property to 0 and using the ILU(n)
 *            preconditioner. The reason for the existence of ILU0 is
 *            that it is computationally cheaper because it does not
 *            need to consider things which are only required for
 *            higher orders
 */
template <class TypeTag>
class ParallelBiCGStabSolverBackend : public ParallelBaseBackend<TypeTag>
{
    using ParentType = ParallelBaseBackend<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

    using ParallelOperator = typename ParentType::ParallelOperator;
    using OverlappingVector = typename ParentType::OverlappingVector;
    using ParallelPreconditioner = typename ParentType::ParallelPreconditioner;
    using ParallelScalarProduct = typename ParentType::ParallelScalarProduct;

    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;

    using RawLinearSolver = BiCGStabSolver<ParallelOperator,
                                           OverlappingVector,
                                           ParallelPreconditioner>;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock> >::value,
                  "The ParallelIstlSolverBackend linear solver backend requires the IstlSparseMatrixAdapter");

public:
    ParallelBiCGStabSolverBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    static void registerParameters()
    {
        ParentType::registerParameters();

        Parameters::Register<Parameters::LinearSolverMaxError<Scalar>>
            ("The maximum residual error which the linear solver tolerates"
             " without giving up");
    }

protected:
    friend ParentType;

    std::shared_ptr<RawLinearSolver> prepareSolver_(ParallelOperator& parOperator,
                                                    ParallelScalarProduct& parScalarProduct,
                                                    ParallelPreconditioner& parPreCond)
    {
        const auto& gridView = this->simulator_.gridView();
        using CCC = CombinedCriterion<OverlappingVector, decltype(gridView.comm())>;

        Scalar linearSolverTolerance = Parameters::Get<Parameters::LinearSolverTolerance<Scalar>>();
        Scalar linearSolverAbsTolerance = Parameters::Get<Parameters::LinearSolverAbsTolerance<Scalar>>();
        if(linearSolverAbsTolerance < 0.0)
            linearSolverAbsTolerance = this->simulator_.model().newtonMethod().tolerance() / 100.0;

        convCrit_.reset(new CCC(gridView.comm(),
                                /*residualReductionTolerance=*/linearSolverTolerance,
                                /*absoluteResidualTolerance=*/linearSolverAbsTolerance,
                                Parameters::Get<Parameters::LinearSolverMaxError<Scalar>>()));

        auto bicgstabSolver =
            std::make_shared<RawLinearSolver>(parPreCond, *convCrit_, parScalarProduct);

        int verbosity = 0;
        if (parOperator.overlap().myRank() == 0)
            verbosity = Parameters::Get<Parameters::LinearSolverVerbosity>();
        bicgstabSolver->setVerbosity(verbosity);
        bicgstabSolver->setMaxIterations(Parameters::Get<Parameters::LinearSolverMaxIterations>());
        bicgstabSolver->setLinearOperator(&parOperator);
        bicgstabSolver->setRhs(this->overlappingb_);

        return bicgstabSolver;
    }

    std::pair<bool,int> runSolver_(std::shared_ptr<RawLinearSolver> solver)
    {
        bool converged = solver->apply(*this->overlappingx_);
        return std::make_pair(converged, int(solver->report().iterations()));
    }

    void cleanupSolver_()
    { /* nothing to do */ }

    std::unique_ptr<ConvergenceCriterion<OverlappingVector> > convCrit_;
};

} // namespace Opm::Linear

#endif
