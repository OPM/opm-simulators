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

#include "linalgproperties.hh"
#include "parallelbasebackend.hh"
#include "bicgstabsolver.hh"
#include "combinedcriterion.hh"
#include "istlsparsematrixadapter.hh"

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <memory>
#include <tuple>
#include <utility>

namespace Opm::Linear {

template <class TypeTag>
class ParallelAmgBackend;

} // namespace Opm::Linear

namespace Opm::Properties {

// Create new type tags
namespace TTag {

struct ParallelAmgLinearSolver
{ using InheritsFrom = std::tuple<ParallelBaseLinearSolver>; };

} // end namespace TTag

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::ParallelAmgLinearSolver>
{ using type = Opm::Linear::ParallelAmgBackend<TypeTag>; };

} // namespace Opm::Properties

namespace Opm::Parameters {

//! The target number of DOFs per processor for the parallel algebraic
//! multi-grid solver
struct AmgCoarsenTarget { static constexpr int value = 5000; };

}

namespace Opm::Linear {

/*!
 * \ingroup Linear
 *
 * \brief Provides a linear solver backend using the parallel
 *        algebraic multi-grid (AMG) linear solver from DUNE-ISTL.
 */
template <class TypeTag>
class ParallelAmgBackend : public ParallelBaseBackend<TypeTag>
{
    using ParentType = ParallelBaseBackend<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LinearSolverScalar = GetPropType<TypeTag, Properties::LinearSolverScalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Overlap = GetPropType<TypeTag, Properties::Overlap>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

    using ParallelOperator = typename ParentType::ParallelOperator;
    using OverlappingVector = typename ParentType::OverlappingVector;
    using ParallelPreconditioner = typename ParentType::ParallelPreconditioner;
    using ParallelScalarProduct = typename ParentType::ParallelScalarProduct;

    static constexpr int numEq = getPropValue<TypeTag, Properties::NumEq>();
    using VectorBlock = Dune::FieldVector<LinearSolverScalar, numEq>;
    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;
    using IstlMatrix = typename SparseMatrixAdapter::IstlMatrix;

    using Vector = Dune::BlockVector<VectorBlock>;

    // define the smoother used for the AMG and specify its
    // arguments
    using SequentialSmoother = Dune::SeqSOR<IstlMatrix, Vector, Vector>;
// using SequentialSmoother = Dune::SeqSSOR<IstlMatrix,Vector,Vector>;
// using SequentialSmoother = Dune::SeqJac<IstlMatrix,Vector,Vector>;
// using SequentialSmoother = Dune::SeqILU<IstlMatrix,Vector,Vector>;

#if HAVE_MPI
    using OwnerOverlapCopyCommunication = Dune::OwnerOverlapCopyCommunication<Opm::Linear::Index>;
    using FineOperator = Dune::OverlappingSchwarzOperator<IstlMatrix,
                                                          Vector,
                                                          Vector,
                                                          OwnerOverlapCopyCommunication>;
    using FineScalarProduct = Dune::OverlappingSchwarzScalarProduct<Vector,
                                                                    OwnerOverlapCopyCommunication>;
    using ParallelSmoother = Dune::BlockPreconditioner<Vector,
                                                       Vector,
                                                       OwnerOverlapCopyCommunication,
                                                       SequentialSmoother>;
    using AMG = Dune::Amg::AMG<FineOperator,
                               Vector,
                               ParallelSmoother,
                               OwnerOverlapCopyCommunication>;
#else
    using FineOperator = Dune::MatrixAdapter<IstlMatrix, Vector, Vector>;
    using FineScalarProduct = Dune::SeqScalarProduct<Vector>;
    using ParallelSmoother = SequentialSmoother;
    using AMG = Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother>;
#endif

    using RawLinearSolver = BiCGStabSolver<ParallelOperator,
                                           OverlappingVector,
                                           AMG> ;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock> >::value,
                  "The ParallelAmgBackend linear solver backend requires the IstlSparseMatrixAdapter");

public:
    ParallelAmgBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    static void registerParameters()
    {
        ParentType::registerParameters();

        Parameters::Register<Parameters::LinearSolverMaxError<Scalar>>
            ("The maximum residual error which the linear solver tolerates "
             "without giving up");
        Parameters::Register<Parameters::AmgCoarsenTarget>
            ("The coarsening target for the agglomerations of "
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
        using CCC = CombinedCriterion<OverlappingVector, decltype(gridView.comm())>;

        Scalar linearSolverTolerance = Parameters::Get<Parameters::LinearSolverTolerance<Scalar>>();
        Scalar linearSolverAbsTolerance = Parameters::Get<Parameters::LinearSolverAbsTolerance<Scalar>>();
        if (linearSolverAbsTolerance < 0.0) {
            linearSolverAbsTolerance = this->simulator_.model().newtonMethod().tolerance()/100.0;
        }

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

#if HAVE_MPI
    template <class ParallelIndexSet>
    void setupAmgIndexSet_(const Overlap& overlap, ParallelIndexSet& istlIndices)
    {
        using GridAttributes = Dune::OwnerOverlapCopyAttributeSet;
        using GridAttributeSet = Dune::OwnerOverlapCopyAttributeSet::AttributeSet;

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

    // trailing return type with decltype used for detecting existence of setUseFixedOrder member function by overloading the setUseFixedOrder function
    template <typename C>
    auto setUseFixedOrder(C criterion, bool booleanValue) -> decltype(criterion.setUseFixedOrder(booleanValue))
    {
        return criterion.setUseFixedOrder(booleanValue); // Set flag to ensure that the matrices in the AMG hierarchy are constructed with deterministic indices.
    }
    template <typename C>
    void setUseFixedOrder(C, ...)
    {
        // do nothing, since the function setUseFixedOrder does not exist yet
    }

    void setupAmg_()
    {
        if (amg_)
            amg_.reset();

        int verbosity = 0;
        if (this->simulator_.vanguard().gridView().comm().rank() == 0)
            verbosity = Parameters::Get<Parameters::LinearSolverVerbosity>();

        using SmootherArgs = typename Dune::Amg::SmootherTraits<ParallelSmoother>::Arguments;

        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1.0;

        // specify the coarsen criterion:
        //
        // using CoarsenCriterion =
        // Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<IstlMatrix,
        //                             Dune::Amg::FirstDiagonal>>
        using CoarsenCriterion = Dune::Amg::
            CoarsenCriterion<Dune::Amg::SymmetricCriterion<IstlMatrix, Dune::Amg::FrobeniusNorm> >;
        int coarsenTarget = Parameters::Get<Parameters::AmgCoarsenTarget>();
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
        setUseFixedOrder(coarsenCriterion, true); // If possible, set flag to ensure that the matrices in the AMG hierarchy are constructed with deterministic indices.

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

} // namespace Opm::Linear

#endif
