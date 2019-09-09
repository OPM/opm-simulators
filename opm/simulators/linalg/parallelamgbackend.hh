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
 * \copydoc Opm::Linear::ParallelAmgBackend
 */
#ifndef EWOMS_PARALLEL_AMG_BACKEND_HH
#define EWOMS_PARALLEL_AMG_BACKEND_HH

#include "parallelbasebackend.hh"
#include "bicgstabsolver.hh"
#include "combinedcriterion.hh"
#include "istlsparsematrixadapter.hh"

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/common/version.hh>

#include <iostream>

namespace Opm {
namespace Linear {
template <class TypeTag>
class ParallelAmgBackend;
}}

BEGIN_PROPERTIES

NEW_TYPE_TAG(ParallelAmgLinearSolver, INHERITS_FROM(ParallelBaseLinearSolver));

NEW_PROP_TAG(AmgCoarsenTarget);
NEW_PROP_TAG(LinearSolverMaxError);

//! The target number of DOFs per processor for the parallel algebraic
//! multi-grid solver
SET_INT_PROP(ParallelAmgLinearSolver, AmgCoarsenTarget, 5000);

SET_SCALAR_PROP(ParallelAmgLinearSolver, LinearSolverMaxError, 1e7);

SET_TYPE_PROP(ParallelAmgLinearSolver, LinearSolverBackend,
              Opm::Linear::ParallelAmgBackend<TypeTag>);

END_PROPERTIES

namespace Opm {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides a linear solver backend using the parallel
 *        algebraic multi-grid (AMG) linear solver from DUNE-ISTL.
 */
template <class TypeTag>
class ParallelAmgBackend : public ParallelBaseBackend<TypeTag>
{
    typedef ParallelBaseBackend<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

    typedef typename ParentType::ParallelOperator ParallelOperator;
    typedef typename ParentType::OverlappingVector OverlappingVector;
    typedef typename ParentType::ParallelPreconditioner ParallelPreconditioner;
    typedef typename ParentType::ParallelScalarProduct ParallelScalarProduct;

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef Dune::FieldVector<LinearSolverScalar, numEq> VectorBlock;
    typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlock;
    typedef typename SparseMatrixAdapter::IstlMatrix IstlMatrix;

    typedef Dune::BlockVector<VectorBlock> Vector;

    // define the smoother used for the AMG and specify its
    // arguments
    typedef Dune::SeqSOR<IstlMatrix, Vector, Vector> SequentialSmoother;
// typedef Dune::SeqSSOR<IstlMatrix,Vector,Vector> SequentialSmoother;
// typedef Dune::SeqJac<IstlMatrix,Vector,Vector> SequentialSmoother;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
// typedef Dune::SeqILU<IstlMatrix,Vector,Vector> SequentialSmoother;
#else
// typedef Dune::SeqILU0<IstlMatrix,Vector,Vector> SequentialSmoother;
// typedef Dune::SeqILUn<IstlMatrix,Vector,Vector> SequentialSmoother;
#endif

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Opm::Linear::Index>
    OwnerOverlapCopyCommunication;
    typedef Dune::OverlappingSchwarzOperator<IstlMatrix,
                                             Vector,
                                             Vector,
                                             OwnerOverlapCopyCommunication> FineOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<Vector,
                                                  OwnerOverlapCopyCommunication> FineScalarProduct;
    typedef Dune::BlockPreconditioner<Vector,
                                      Vector,
                                      OwnerOverlapCopyCommunication,
                                      SequentialSmoother> ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator,
                           Vector,
                           ParallelSmoother,
                           OwnerOverlapCopyCommunication> AMG;
#else
    typedef Dune::MatrixAdapter<IstlMatrix, Vector, Vector> FineOperator;
    typedef Dune::SeqScalarProduct<Vector> FineScalarProduct;
    typedef SequentialSmoother ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother> AMG;
#endif

    typedef BiCGStabSolver<ParallelOperator,
                           OverlappingVector,
                           AMG> RawLinearSolver;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock> >::value,
                  "The ParallelAmgBackend linear solver backend requires the IstlSparseMatrixAdapter");

public:
    ParallelAmgBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverMaxError,
                             "The maximum residual error which the linear solver tolerates"
                             " without giving up");
        EWOMS_REGISTER_PARAM(TypeTag, int, AmgCoarsenTarget,
                             "The coarsening target for the agglomerations of "
                             "the AMG preconditioner");
    }

protected:
    friend ParentType;

    std::shared_ptr<AMG> preparePreconditioner_()
    {
#if HAVE_MPI
        // create and initialize DUNE's OwnerOverlapCopyCommunication
        // using the domestic overlap
        istlComm_ = std::make_shared<OwnerOverlapCopyCommunication>(MPI_COMM_WORLD);
        setupAmgIndexSet_(this->overlappingMatrix_->overlap(), istlComm_->indexSet());
        istlComm_->remoteIndices().template rebuild<false>();
#endif

        // create the parallel scalar product and the parallel operator
#if HAVE_MPI
        fineOperator_ = std::make_shared<FineOperator>(*this->overlappingMatrix_, *istlComm_);
#else
        fineOperator_ = std::make_shared<FineOperator>(*this->overlappingMatrix_);
#endif

        setupAmg_();

        return amg_;
    }

    void cleanupPreconditioner_()
    { /* nothing to do */ }

    std::shared_ptr<RawLinearSolver> prepareSolver_(ParallelOperator& parOperator,
                                                    ParallelScalarProduct& parScalarProduct,
                                                    AMG& parPreCond)
    {
        const auto& gridView = this->simulator_.gridView();
        typedef CombinedCriterion<OverlappingVector, decltype(gridView.comm())> CCC;

        Scalar linearSolverTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        Scalar linearSolverAbsTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance);
        if(linearSolverAbsTolerance < 0.0)
            linearSolverAbsTolerance = this->simulator_.model().newtonMethod().tolerance()/100.0;

        convCrit_.reset(new CCC(gridView.comm(),
                                /*residualReductionTolerance=*/linearSolverTolerance,
                                /*absoluteResidualTolerance=*/linearSolverAbsTolerance,
                                EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverMaxError)));

        auto bicgstabSolver =
            std::make_shared<RawLinearSolver>(parPreCond, *convCrit_, parScalarProduct);

        int verbosity = 0;
        if (parOperator.overlap().myRank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        bicgstabSolver->setVerbosity(verbosity);
        bicgstabSolver->setMaxIterations(EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations));
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

#if HAVE_MPI
    template <class ParallelIndexSet>
    void setupAmgIndexSet_(const Overlap& overlap, ParallelIndexSet& istlIndices)
    {
        typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
        typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet GridAttributeSet;

        // create DUNE's ParallelIndexSet from a domestic overlap
        istlIndices.beginResize();
        for (Index curIdx = 0; static_cast<size_t>(curIdx) < overlap.numDomestic(); ++curIdx) {
            GridAttributeSet gridFlag =
                overlap.iAmMasterOf(curIdx)
                ? GridAttributes::owner
                : GridAttributes::copy;

            // an index is used by other processes if it is in the
            // domestic or in the foreign overlap.
            bool isShared = overlap.isInOverlap(curIdx);

            assert(curIdx == overlap.globalToDomestic(overlap.domesticToGlobal(curIdx)));
            istlIndices.add(/*globalIdx=*/overlap.domesticToGlobal(curIdx),
                            Dune::ParallelLocalIndex<GridAttributeSet>(static_cast<size_t>(curIdx),
                                                                       gridFlag,
                                                                       isShared));
        }
        istlIndices.endResize();
    }
#endif

    void setupAmg_()
    {
        if (amg_)
            amg_.reset();

        int verbosity = 0;
        if (this->simulator_.vanguard().gridView().comm().rank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        typedef typename Dune::Amg::SmootherTraits<ParallelSmoother>::Arguments SmootherArgs;

        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1.0;

        // specify the coarsen criterion:
        //
        // typedef
        // Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<IstlMatrix,
        //                             Dune::Amg::FirstDiagonal>>
        typedef Dune::Amg::
            CoarsenCriterion<Dune::Amg::SymmetricCriterion<IstlMatrix, Dune::Amg::FrobeniusNorm> >
            CoarsenCriterion;
        int coarsenTarget = EWOMS_GET_PARAM(TypeTag, int, AmgCoarsenTarget);
        CoarsenCriterion coarsenCriterion(/*maxLevel=*/15, coarsenTarget);
        coarsenCriterion.setDefaultValuesAnisotropic(GridView::dimension,
                                                     /*aggregateSizePerDim=*/3);
        if (verbosity > 0)
            coarsenCriterion.setDebugLevel(1);
        else
            coarsenCriterion.setDebugLevel(0); // make the AMG shut up

        // reduce the minium coarsen rate (default is 1.2)
        coarsenCriterion.setMinCoarsenRate(1.05);
        // coarsenCriterion.setAccumulate(Dune::Amg::noAccu);
        coarsenCriterion.setAccumulate(Dune::Amg::atOnceAccu);
        coarsenCriterion.setSkipIsolated(false);

// instantiate the AMG preconditioner
#if HAVE_MPI
        amg_ = std::make_shared<AMG>(*fineOperator_, coarsenCriterion, smootherArgs, *istlComm_);
#else
        amg_ = std::make_shared<AMG>(*fineOperator_, coarsenCriterion, smootherArgs);
#endif
    }

    std::unique_ptr<ConvergenceCriterion<OverlappingVector> > convCrit_;

    std::shared_ptr<FineOperator> fineOperator_;
    std::shared_ptr<AMG> amg_;

#if HAVE_MPI
    std::shared_ptr<OwnerOverlapCopyCommunication> istlComm_;
#endif
};

} // namespace Linear
} // namespace Opm

#endif
