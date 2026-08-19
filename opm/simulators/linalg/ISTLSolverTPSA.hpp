// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef ISTL_SOLVER_TPSA_HPP
#define ISTL_SOLVER_TPSA_HPP

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solver.hh>

#include <opm/common/CriticalError.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/linalg/AbstractISTLSolver.hpp>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>
#include <opm/simulators/linalg/tpsa/TpsaPreconditionerFactory.hpp>
#include <opm/simulators/linalg/tpsa/TpsaTypes.hpp>

#include <algorithm>
#include <any>
#include <cctype>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \brief Linear solver for the field-split TPSA system.
 *
 * The Jacobian handed in by TpsaLinearizer is a Linear::TpsaMatrix, which
 * already stores the system as the 19 sub-matrices of the five-field split
 * (see TpsaTypes.hpp), and the residual is a Linear::TpsaVector wrapping the
 * matching Dune::MultiTypeBlockVector.
 */
template <class TypeTag>
class ISTLSolverTPSA :
    public AbstractISTLSolver<GetPropType<TypeTag, Properties::SparseMatrixAdapterTPSA>,
                              GetPropType<TypeTag, Properties::GlobalEqVectorTPSA> >
{
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapterTPSA>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVectorTPSA>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using MultiVector = Linear::TpsaMultiVector<Scalar>;

    using SolverType = Dune::InverseOperator<MultiVector, MultiVector>;
    using PrecondType = Dune::PreconditionerWithUpdate<MultiVector, MultiVector>;

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    //! \brief No pressure equation to single out in the mechanics system.
    static constexpr std::size_t pressureIndex = 0;

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Simulator object
    */
    explicit ISTLSolverTPSA(const Simulator& simulator)
        : simulator_(simulator)
          , solveCount_(0)
          , iterations_(0)
          , matrix_(nullptr)
          , rhs_(nullptr)
    {
        // Init parameters
        parameters_.init();

        // Initialize linear solver
        initialize();
    }

    /*!
    * \brief Register runtime/default parameters for linear solver
    */
    static void registerParameters()
    {
        TpsaLinearSolverParameters::registerParameters();
    }

    /*!
    * \brief Setup linear solver object based on runtime/default parameters
    */
    void initialize()
    {
        // Setup property tree for FlexibleSolver
        prm_ = setupPropertyTree(parameters_, false, false, /*tpsaSetup=*/true);

        // Reset comm_ pointer
#if HAVE_MPI
        comm_.reset(new CommunicationType(simulator_.vanguard().grid().comm()));
#endif

        // Extract and copy parallel grid information
        extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
#if HAVE_MPI
        if (isParallel()) {
            const std::size_t size = simulator_.vanguard().grid().leafGridView().size(0);
            detail::copyParValues(parallelInformation_, size, *comm_);
        }
#endif

        // Get info on overlapping rows
        ElementMapper elemMapper(simulator_.vanguard().gridView(), Dune::mcmgElementLayout());
        std::vector<int> dummyInteriorRows;
        detail::findOverlapAndInterior(simulator_.vanguard().grid(),
                                       elemMapper,
                                       overlapRows_,
                                       dummyInteriorRows);

        detail::printLinearSolverParameters(parameters_, prm_, simulator_.gridView().comm());
    }

    /*!
    * \brief Prepare matix and rhs vector for linear solve
    *
    * \param M System matrix
    * \param b Right-hand side vector
    */
    void initPrepare(const SparseMatrixAdapter& M, Vector& b)
    {
        if (matrix_ == nullptr) {
            // The model reuses the same matrix object for the whole run; keep a
            // pointer to it. The const_cast mirrors ISTLSolver: the solver has
            // to be able to modify the assembled system (overlap rows).
            matrix_ = const_cast<SparseMatrixAdapter*>(&M);
        } else if (&M != matrix_) {
            OPM_THROW(std::logic_error,
                      "TPSA: Matrix objects are expected to be reused when reassembling!");
        }

        // Set right-hand side vector
        rhs_ = &b;

        // Zero out the overlapping cells (not for the overlapping ILU variant,which handles them
        // itself).
        auto type = prm_.get<std::string>("preconditioner.type", "paroverilu0");
        std::ranges::transform(type, type.begin(), ::tolower);
        if (isParallel() && type != "paroverilu0") {
            matrix_->makeOverlapRowsInvalid(overlapRows_);
        }
    }

    /*!
    * \copydoc AbstractISTLSolver::prepare
    */
    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        OPM_TIMEBLOCK(istlSolverPrepare);
        try {
            initPrepare(M, b);
            prepareSolver();
        }
        OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR
            ("TPSA: Failure likely due to a faulty linear solver JSON specification. "
             "Check for errors related to missing nodes.");
    }

    /*!
    * \copydoc AbstractISTLSolver::prepare
    */
    void prepare(const Matrix& /*M*/, Vector& /*b*/) override
    {
        // The five-field view is not enough to prepare the solve: the overlap
        // handling needs the owning TpsaMatrix
        OPM_THROW(std::logic_error,
                  "TPSA: prepare() needs the TpsaMatrix itself, not just its field-split view");
    }

    /*!
    * \copydoc AbstractISTLSolver::solve
    */
    bool solve(Vector& x) override
    {
        OPM_TIMEBLOCK(istlSolverSolve);
        ++solveCount_;

        x = 0.0;

        // Solve linear system
        Dune::InverseOperatorResult result;
        assert(solver_);
        solver_->apply(x.istlVector(), rhs_->istlVector(), result);

        // Store no. linear iterations
        iterations_ = result.iterations;

        // Return result for convergence check (boolean)
        return checkConvergence(result);
    }

    /*!
    * \brief Reset number of solver calls to zero
    */
    void resetSolveCount()
    {
        solveCount_ = 0;
    }

    // ////
    // Unused AbstractISTLSolver functions
    // ////
    /*!
    * \copydoc AbstractISTLSolver::eraseMatrix
    */
    void eraseMatrix() override
    {
    }

    /*!
    * \copydoc AbstractISTLSolver::setActiveSolver
    */
    void setActiveSolver(int /*num*/) override
    {
    }

    /*!
    * \copydoc AbstractISTLSolver::numAvailableSolvers
    */
    int numAvailableSolvers() const override
    {
        return 1;
    }

    /*!
    * \copydoc AbstractISTLSolver::setResidual
    */
    void setResidual(Vector& /*b*/) override
    {
    }

    /*!
    * \copydoc AbstractISTLSolver::setMatrix
    */
    void setMatrix(const SparseMatrixAdapter& /*M*/) override
    {
    }

    // ///
    // Public get functions
    // ///
    /*!
    * \copydoc AbstractISTLSolver::getResidual
    */
    void getResidual(Vector& b) const override
    {
        b = *rhs_;
    }

    /*!
    * \copydoc AbstractISTLSolver::getSolveCount
    */
    int getSolveCount() const override
    {
        return solveCount_;
    }

    /*!
    * \copydoc AbstractISTLSolver::iterations
    */
    int iterations() const override
    {
        return iterations_;
    }

    /*!
    * \copydoc AbstractISTLSolver::comm
    */
    const CommunicationType* comm() const override
    {
        return comm_.get();
    }

protected:
    /*!
    * \brief Check for parallel session
    *
    * \returns Bool indicating if this is a parallel session
    *
    * \warning comm_ must be set before using this function
    */
    bool isParallel() const
    {
#if HAVE_MPI
        return comm_->communicator().size() > 1;
#else
        return false;
#endif
    }

    /*!
    * \brief Check for linear solver convergence
    *
    * \param result Linear solver result container
    * \returns Bool indicating convergence
    */
    bool checkConvergence(const Dune::InverseOperatorResult& result) const
    {
        return AbstractISTLSolver<SparseMatrixAdapter, Vector>::checkConvergence(result, parameters_);
    }

    // ///
    // Protected get functions
    // ///
    /*!
     * \brief Create the solver on the first call, refresh the preconditioner on
     *        every later one.
     */
    void prepareSolver()
    {
        OPM_TIMEBLOCK(flexibleSolverPrepare);

        if (solver_ == nullptr) {
            OPM_TIMEBLOCK(flexibleSolverCreate);
            createSolver();
        } else {
            OPM_TIMEBLOCK(flexibleSolverUpdate);
            precond_->update();
        }
    }


    /*!
     * \brief Build the linear operator and the flexible solver for the TPSA system
     *
     * \warning matrix_ and comm_ must be set before calling this function. The
     *          matrix is handed to the operator by reference, so it must outlive
     *          the solver.
     */
    void createSolver()
    {
        // The mechanics system has no equation weights.
        std::function<MultiVector()> weightCalculator;

        Matrix& matrixView = matrix_->istlMatrix();

        if (isParallel()) {
#if HAVE_MPI
            // All five fields live on the same dof partition, so they share the same
            // communicator.
            multiComm_ = std::make_unique<TpsaComm>(*comm_, *comm_, *comm_, *comm_, *comm_);

            parOperator_ = std::make_unique<TpsaParOperator<Scalar> >(matrixView, *multiComm_);
            parSolver_ = std::make_unique<
                Dune::FlexibleSolver<TpsaParOperator<Scalar> > >(*parOperator_,
                                                                 *multiComm_,
                                                                 prm_,
                                                                 weightCalculator,
                                                                 pressureIndex);

            solver_ = parSolver_.get();
            precond_ = &parSolver_->preconditioner();
#endif
        } else {
            seqOperator_ = std::make_unique<TpsaSeqOperator<Scalar> >(matrixView);
            seqSolver_ = std::make_unique<
                Dune::FlexibleSolver<TpsaSeqOperator<Scalar> > >(*seqOperator_,
                                                                 prm_,
                                                                 weightCalculator,
                                                                 pressureIndex);

            solver_ = seqSolver_.get();
            precond_ = &seqSolver_->preconditioner();
        }
    }

    const Simulator& simulator_;
    std::any parallelInformation_;
    int solveCount_;
    int iterations_;

    TpsaLinearSolverParameters parameters_;
    PropertyTree prm_;

    SparseMatrixAdapter* matrix_;
    Vector* rhs_;

    std::shared_ptr<CommunicationType> comm_;
    std::vector<int> overlapRows_;

    // Serial solver components
    std::unique_ptr<TpsaSeqOperator<Scalar> > seqOperator_;
    std::unique_ptr<Dune::FlexibleSolver<TpsaSeqOperator<Scalar> > > seqSolver_;

    // Parallel solver components
#if HAVE_MPI
    std::unique_ptr<TpsaComm> multiComm_;
    std::unique_ptr<TpsaParOperator<Scalar> > parOperator_;
    std::unique_ptr<Dune::FlexibleSolver<TpsaParOperator<Scalar> > > parSolver_;
#endif

    SolverType* solver_ = nullptr;
    PrecondType* precond_ = nullptr;
};

}  // namespace Opm

#endif // ISTL_SOLVER_TPSA_HPP
