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
 * \copydoc Opm::Linear::ParallelBaseBackend
 */
#ifndef EWOMS_PARALLEL_BASE_BACKEND_HH
#define EWOMS_PARALLEL_BASE_BACKEND_HH

#include <opm/common/Exceptions.hpp>

#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/overlappingbcrsmatrix.hh>
#include <opm/simulators/linalg/overlappingblockvector.hh>
#include <opm/simulators/linalg/overlappingpreconditioner.hh>
#include <opm/simulators/linalg/overlappingscalarproduct.hh>
#include <opm/simulators/linalg/overlappingoperator.hh>
#include <opm/simulators/linalg/parallelbasebackend.hh>
#include <opm/simulators/linalg/istlpreconditionerwrappers.hh>
//#include <MatrixMarketSpecializations.hpp>
//#include <dune/istl/matrixmarket.hh>

#include <opm/models/utils/genericguard.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/linalgproperties.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <sstream>
#include <memory>
#include <iostream>

namespace Opm::Properties {

namespace TTag {
struct ParallelBaseLinearSolver {};
}

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
template<class TypeTag>
struct SparseMatrixAdapter<TypeTag, TTag::ParallelBaseLinearSolver>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = Opm::MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Opm::Linear::IstlSparseMatrixAdapter<Block>;
};

} // namespace Opm::Properties

namespace Opm {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides the common code which is required by most linear solvers.
 *
 * This class provides access to all preconditioners offered by dune-istl using the
 * PreconditionerWrapper property:
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
class ParallelBaseBackend
{
protected:
    using Implementation = GetPropType<TypeTag, Properties::LinearSolverBackend>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LinearSolverScalar = GetPropType<TypeTag, Properties::LinearSolverScalar>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using BorderListCreator = GetPropType<TypeTag, Properties::BorderListCreator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using Overlap = GetPropType<TypeTag, Properties::Overlap>;
    using OverlappingVector = GetPropType<TypeTag, Properties::OverlappingVector>;
    using OverlappingMatrix = GetPropType<TypeTag, Properties::OverlappingMatrix>;

    using PreconditionerWrapper = GetPropType<TypeTag, Properties::PreconditionerWrapper>;
    using SequentialPreconditioner = typename PreconditionerWrapper::SequentialPreconditioner;

    using ParallelPreconditioner = Opm::Linear::OverlappingPreconditioner<SequentialPreconditioner, Overlap>;
    using ParallelScalarProduct = Opm::Linear::OverlappingScalarProduct<OverlappingVector, Overlap>;
    using ParallelOperator = Opm::Linear::OverlappingOperator<OverlappingMatrix,
                                                              OverlappingVector,
                                                              OverlappingVector>;

    enum { dimWorld = GridView::dimensionworld };

public:
    ParallelBaseBackend(const Simulator& simulator)
        : simulator_(simulator)
        , gridSequenceNumber_( -1 )
        , lastIterations_( -1 )
    {
        overlappingMatrix_ = nullptr;
        overlappingb_ = nullptr;
        overlappingx_ = nullptr;
    }

    ~ParallelBaseBackend()
    { cleanup_(); }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverTolerance,
                             "The maximum allowed error between of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance,
                             "The maximum accepted error of the norm of the residual");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, LinearSolverOverlapSize,
                             "The size of the algebraic overlap for the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIterations,
                             "The maximum number of iterations of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");

        PreconditionerWrapper::registerParameters();
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { cleanup_(); }

    /*!
     * \brief Set up the internal data structures required for the linear solver.
     *
     * This only specified the topology of the linear system of equations; it does does
     * *not* assign the values of the residual vector and its Jacobian matrix.
     */
    void prepare(const SparseMatrixAdapter& M, const Vector& )
    {
        // if grid has changed the sequence number has changed too
        int curSeqNum = simulator_.vanguard().gridSequenceNumber();
        if (gridSequenceNumber_ == curSeqNum && overlappingMatrix_)
            // the grid has not changed since the overlappingMatrix_has been created, so
            // there's noting to do
            return;

        asImp_().cleanup_();
        gridSequenceNumber_ = curSeqNum;

        BorderListCreator borderListCreator(simulator_.gridView(),
                                            simulator_.model().dofMapper());

        // create the overlapping Jacobian matrix
        unsigned overlapSize = EWOMS_GET_PARAM(TypeTag, unsigned, LinearSolverOverlapSize);
        overlappingMatrix_ = new OverlappingMatrix(M.istlMatrix(),
                                                   borderListCreator.borderList(),
                                                   borderListCreator.blackList(),
                                                   overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

        // writeOverlapToVTK_();
    }

    /*!
     * \brief Assign values to the internal data structure for the residual vector.
     *
     * This method also cares about synchronizing that vector with the peer processes.
     */
    void setResidual(const Vector& b)
    {
        // copy the interior values of the non-overlapping residual vector to the
        // overlapping one
        overlappingb_->assignAddBorder(b);
    }

    /*!
     * \brief Retrieve the synchronized internal residual vector.
     *
     * This only deals with entries which are local to the current process.
     */
    void getResidual(Vector& b) const
    {
        // update the non-overlapping vector with the overlapping one
        overlappingb_->assignTo(b);
    }

    /*!
     * \brief Sets the values of the residual's Jacobian matrix.
     *
     * This method also synchronizes the data structure across the processes which are
     * involved in the simulation run.
     */
    void setMatrix(const SparseMatrixAdapter& M)
    {
        overlappingMatrix_->assignFromNative(M.istlMatrix());
        overlappingMatrix_->syncAdd();
    }

    /*!
     * \brief Actually solve the linear system of equations.
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    bool solve(Vector& x)
    {
        (*overlappingx_) = 0.0;

        auto parPreCond = asImp_().preparePreconditioner_();
        auto precondCleanupFn = [this]() -> void
                                { this->asImp_().cleanupPreconditioner_(); };
        auto precondCleanupGuard = Opm::make_guard(precondCleanupFn);
        // create the parallel scalar product and the parallel operator
        ParallelScalarProduct parScalarProduct(overlappingMatrix_->overlap());
        ParallelOperator parOperator(*overlappingMatrix_);

        //const auto& matrix = parOperator.getMatrix();
        //Dune::storeMatrixMarket(*overlappingMatrix_, std::string("mymatrix.mm"));

        // retrieve the linear solver
        auto solver = asImp_().prepareSolver_(parOperator,
                                              parScalarProduct,
                                              *parPreCond);

        auto cleanupSolverFn =
            [this]() -> void
            { this->asImp_().cleanupSolver_(); };
        GenericGuard<decltype(cleanupSolverFn)> solverGuard(cleanupSolverFn);

        // run the linear solver and have some fun
        auto result = asImp_().runSolver_(solver);
        // store number of iterations used
        lastIterations_ = result.second;

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);
        overlappingx_->print();

        // return the result of the solver
        return result.first;
    }

    /*!
     * \brief Return number of iterations used during last solve.
     */
    size_t iterations () const
    { return lastIterations_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void cleanup_()
    {
        // create the overlapping Jacobian matrix and vectors
        delete overlappingMatrix_;
        delete overlappingb_;
        delete overlappingx_;

        overlappingMatrix_ = 0;
        overlappingb_ = 0;
        overlappingx_ = 0;
    }

    std::shared_ptr<ParallelPreconditioner> preparePreconditioner_()
    {
        int preconditionerIsReady = 1;
        try {
            // update sequential preconditioner
            precWrapper_.prepare(*overlappingMatrix_);
        }
        catch (const Dune::Exception& e) {
            std::cout << "Preconditioner threw exception \"" << e.what()
                      << " on rank " << overlappingMatrix_->overlap().myRank()
                      << "\n"  << std::flush;
            preconditionerIsReady = 0;
        }

        // make sure that the preconditioner is also ready on all peer
        // ranks.
        preconditionerIsReady = simulator_.gridView().comm().min(preconditionerIsReady);
        if (!preconditionerIsReady)
            throw NumericalProblem("Creating the preconditioner failed");

        // create the parallel preconditioner
        return std::make_shared<ParallelPreconditioner>(precWrapper_.get(), overlappingMatrix_->overlap());
    }

    void cleanupPreconditioner_()
    {
        precWrapper_.cleanup();
    }

    void writeOverlapToVTK_()
    {
        for (int lookedAtRank = 0;
             lookedAtRank < simulator_.gridView().comm().size(); ++lookedAtRank) {
            std::cout << "writing overlap for rank " << lookedAtRank << "\n"  << std::flush;
            using VtkField = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
            int n = simulator_.gridView().size(/*codim=*/dimWorld);
            VtkField isInOverlap(n);
            VtkField rankField(n);
            isInOverlap = 0.0;
            rankField = 0.0;
            assert(rankField.two_norm() == 0.0);
            assert(isInOverlap.two_norm() == 0.0);
            auto vIt = simulator_.gridView().template begin</*codim=*/dimWorld>();
            const auto& vEndIt = simulator_.gridView().template end</*codim=*/dimWorld>();
            const auto& overlap = overlappingMatrix_->overlap();
            for (; vIt != vEndIt; ++vIt) {
                int nativeIdx = simulator_.model().vertexMapper().map(*vIt);
                int localIdx = overlap.foreignOverlap().nativeToLocal(nativeIdx);
                if (localIdx < 0)
                    continue;
                rankField[nativeIdx] = simulator_.gridView().comm().rank();
                if (overlap.peerHasIndex(lookedAtRank, localIdx))
                    isInOverlap[nativeIdx] = 1.0;
            }

            using VtkWriter = Dune::VTKWriter<GridView>;
            VtkWriter writer(simulator_.gridView(), Dune::VTK::conforming);
            writer.addVertexData(isInOverlap, "overlap");
            writer.addVertexData(rankField, "rank");

            std::ostringstream oss;
            oss << "overlap_rank=" << lookedAtRank;
            writer.write(oss.str().c_str(), Dune::VTK::ascii);
        }
    }

    const Simulator& simulator_;
    int gridSequenceNumber_;
    size_t lastIterations_;

    OverlappingMatrix *overlappingMatrix_;
    OverlappingVector *overlappingb_;
    OverlappingVector *overlappingx_;

    PreconditionerWrapper precWrapper_;
};
}} // namespace Linear, Opm

namespace Opm::Properties {

//! make the linear solver shut up by default
template<class TypeTag>
struct LinearSolverVerbosity<TypeTag, TTag::ParallelBaseLinearSolver> { static constexpr int value = 0; };

//! set the preconditioner relaxation parameter to 1.0 by default
template<class TypeTag>
struct PreconditionerRelaxation<TypeTag, TTag::ParallelBaseLinearSolver>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

//! set the preconditioner order to 0 by default
template<class TypeTag>
struct PreconditionerOrder<TypeTag, TTag::ParallelBaseLinearSolver> { static constexpr int value = 0; };

//! by default use the same kind of floating point values for the linearization and for
//! the linear solve
template<class TypeTag>
struct LinearSolverScalar<TypeTag, TTag::ParallelBaseLinearSolver>
{ using type = GetPropType<TypeTag, Properties::Scalar>; };

template<class TypeTag>
struct OverlappingMatrix<TypeTag, TTag::ParallelBaseLinearSolver>
{
private:
    static constexpr int numEq = getPropValue<TypeTag, Properties::NumEq>();
    using LinearSolverScalar = GetPropType<TypeTag, Properties::LinearSolverScalar>;
    using MatrixBlock = Opm::MatrixBlock<LinearSolverScalar, numEq, numEq>;
    using NonOverlappingMatrix = Dune::BCRSMatrix<MatrixBlock>;

public:
    using type = Opm::Linear::OverlappingBCRSMatrix<NonOverlappingMatrix>;
};

template<class TypeTag>
struct Overlap<TypeTag, TTag::ParallelBaseLinearSolver>
{ using type = typename GetPropType<TypeTag, Properties::OverlappingMatrix>::Overlap; };

template<class TypeTag>
struct OverlappingVector<TypeTag, TTag::ParallelBaseLinearSolver>
{
    static constexpr int numEq = getPropValue<TypeTag, Properties::NumEq>();
    using LinearSolverScalar = GetPropType<TypeTag, Properties::LinearSolverScalar>;
    using VectorBlock = Dune::FieldVector<LinearSolverScalar, numEq>;
    using Overlap = GetPropType<TypeTag, Properties::Overlap>;
    using type = Opm::Linear::OverlappingBlockVector<VectorBlock, Overlap>;
};

template<class TypeTag>
struct OverlappingScalarProduct<TypeTag, TTag::ParallelBaseLinearSolver>
{
    using OverlappingVector = GetPropType<TypeTag, Properties::OverlappingVector>;
    using Overlap = GetPropType<TypeTag, Properties::Overlap>;
    using type = Opm::Linear::OverlappingScalarProduct<OverlappingVector, Overlap>;
};

template<class TypeTag>
struct OverlappingLinearOperator<TypeTag, TTag::ParallelBaseLinearSolver>
{
    using OverlappingMatrix = GetPropType<TypeTag, Properties::OverlappingMatrix>;
    using OverlappingVector = GetPropType<TypeTag, Properties::OverlappingVector>;
    using type = Opm::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector,
                                                  OverlappingVector>;
};

template<class TypeTag>
struct PreconditionerWrapper<TypeTag, TTag::ParallelBaseLinearSolver>
{ using type = Opm::Linear::PreconditionerWrapperILU<TypeTag>; };

//! set the default overlap size to 2
template<class TypeTag>
struct LinearSolverOverlapSize<TypeTag, TTag::ParallelBaseLinearSolver> { static constexpr int value = 2; };

//! set the default number of maximum iterations for the linear solver
template<class TypeTag>
struct LinearSolverMaxIterations<TypeTag, TTag::ParallelBaseLinearSolver> { static constexpr int value = 1000; };

} // namespace Opm::Properties

#endif
