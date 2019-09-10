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

#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/overlappingbcrsmatrix.hh>
#include <opm/simulators/linalg/overlappingblockvector.hh>
#include <opm/simulators/linalg/overlappingpreconditioner.hh>
#include <opm/simulators/linalg/overlappingscalarproduct.hh>
#include <opm/simulators/linalg/overlappingoperator.hh>
#include <opm/simulators/linalg/parallelbasebackend.hh>
#include <opm/simulators/linalg/istlpreconditionerwrappers.hh>

#include <ewoms/common/genericguard.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <sstream>
#include <memory>
#include <iostream>

BEGIN_PROPERTIES
NEW_TYPE_TAG(ParallelBaseLinearSolver);

// forward declaration of the required property tags
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(GridView);

NEW_PROP_TAG(BorderListCreator);
NEW_PROP_TAG(Overlap);
NEW_PROP_TAG(OverlappingVector);
NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingScalarProduct);
NEW_PROP_TAG(OverlappingLinearOperator);

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);

//! the preconditioner used by the linear solver
NEW_PROP_TAG(PreconditionerWrapper);


//! The floating point type used internally by the linear solver
NEW_PROP_TAG(LinearSolverScalar);

/*!
 * \brief The size of the algebraic overlap of the linear solver.
 *
 * Algebraic overlaps can be thought as being the same as the overlap
 * of a grid, but it is only existant for the linear system of
 * equations.
 */
NEW_PROP_TAG(LinearSolverOverlapSize);

/*!
 * \brief Maximum accepted error of the solution of the linear solver.
 */
NEW_PROP_TAG(LinearSolverTolerance);

/*!
 * \brief Maximum accepted error of the norm of the residual.
 */
NEW_PROP_TAG(LinearSolverAbsTolerance);

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LinearSolverVerbosity);

//! Maximum number of iterations eyecuted by the linear solver
NEW_PROP_TAG(LinearSolverMaxIterations);

//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
SET_PROP(ParallelBaseLinearSolver, SparseMatrixAdapter)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Opm::MatrixBlock<Scalar, numEq, numEq> Block;

public:
    typedef typename Opm::Linear::IstlSparseMatrixAdapter<Block> type;
};

END_PROPERTIES

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
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, BorderListCreator) BorderListCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PreconditionerWrapper) PreconditionerWrapper;
    typedef typename PreconditionerWrapper::SequentialPreconditioner SequentialPreconditioner;

    typedef Opm::Linear::OverlappingPreconditioner<SequentialPreconditioner, Overlap> ParallelPreconditioner;
    typedef Opm::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> ParallelScalarProduct;
    typedef Opm::Linear::OverlappingOperator<OverlappingMatrix,
                                             OverlappingVector,
                                             OverlappingVector> ParallelOperator;

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
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
        Dune::FMatrixPrecision<LinearSolverScalar>::set_singular_limit(1.e-30);
        Dune::FMatrixPrecision<LinearSolverScalar>::set_absolute_limit(1.e-30);
#endif

        (*overlappingx_) = 0.0;

        auto parPreCond = asImp_().preparePreconditioner_();
        auto precondCleanupFn = [this]() -> void
                                { this->asImp_().cleanupPreconditioner_(); };
        auto precondCleanupGuard = Opm::make_guard(precondCleanupFn);
        // create the parallel scalar product and the parallel operator
        ParallelScalarProduct parScalarProduct(overlappingMatrix_->overlap());
        ParallelOperator parOperator(*overlappingMatrix_);

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
            throw Opm::NumericalIssue("Creating the preconditioner failed");

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
            typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > VtkField;
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

            typedef Dune::VTKWriter<GridView> VtkWriter;
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
}} // namespace Linear, Ewoms

BEGIN_PROPERTIES

//! make the linear solver shut up by default
SET_INT_PROP(ParallelBaseLinearSolver, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(ParallelBaseLinearSolver, PreconditionerRelaxation, 1.0);

//! set the preconditioner order to 0 by default
SET_INT_PROP(ParallelBaseLinearSolver, PreconditionerOrder, 0);

//! by default use the same kind of floating point values for the linearization and for
//! the linear solve
SET_TYPE_PROP(ParallelBaseLinearSolver,
              LinearSolverScalar,
              typename GET_PROP_TYPE(TypeTag, Scalar));

SET_PROP(ParallelBaseLinearSolver, OverlappingMatrix)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef Opm::MatrixBlock<LinearSolverScalar, numEq, numEq> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> NonOverlappingMatrix;

public:
    typedef Opm::Linear::OverlappingBCRSMatrix<NonOverlappingMatrix> type;
};

SET_TYPE_PROP(ParallelBaseLinearSolver,
              Overlap,
              typename GET_PROP_TYPE(TypeTag, OverlappingMatrix)::Overlap);

SET_PROP(ParallelBaseLinearSolver, OverlappingVector)
{
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef Dune::FieldVector<LinearSolverScalar, numEq> VectorBlock;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Opm::Linear::OverlappingBlockVector<VectorBlock, Overlap> type;
};

SET_PROP(ParallelBaseLinearSolver, OverlappingScalarProduct)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Opm::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> type;
};

SET_PROP(ParallelBaseLinearSolver, OverlappingLinearOperator)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef Opm::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector,
                                             OverlappingVector> type;
};

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
SET_TYPE_PROP(ParallelBaseLinearSolver,
              PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILU<TypeTag>);
#else
SET_TYPE_PROP(ParallelBaseLinearSolver,
              PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILU0<TypeTag>);
#endif

//! set the default overlap size to 2
SET_INT_PROP(ParallelBaseLinearSolver, LinearSolverOverlapSize, 2);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(ParallelBaseLinearSolver, LinearSolverMaxIterations, 1000);

END_PROPERTIES

#endif
