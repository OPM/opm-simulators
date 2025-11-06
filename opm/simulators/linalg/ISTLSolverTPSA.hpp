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

#include <opm/common/CriticalError.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTreeTPSA.hpp>
#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>

#include <algorithm>
#include <any>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>


namespace Opm {

/*!
* \brief Class for setting up ISTL linear solvers for TPSA
*
* \note Most of the code is copied from ISTLSolver class
*/
template <class TypeTag>
class ISTLSolverTPSA
{
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapterTPSA>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVectorTPSA>;

    using Matrix = typename SparseMatrixAdapter::IstlMatrix;


#if HAVE_MPI
        using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
        using CommunicationType = Dune::Communication<int>;
#endif

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Simulator object
    */
    ISTLSolverTPSA(const Simulator& simulator)
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
        prm_ = setupPropertyTreeTPSA(parameters_);

        // Reset comm_ pointer
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
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
        detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, dummyInteriorRows);

        // Set number of interior cells in FlexibleSolverInfo
        flexibleSolver_.interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), true);

        // Print parameters to PRT/DBG logs.
        detail::printLinearSolverParameters(parameters_, prm_, simulator_.gridView().comm());
    }

    /*!
    * \brief Prepare matix and rhs vector for linear solve
    *
    * \param M System matrix
    * \param b Right-hand side vector
    */
    void initPrepare(const Matrix& M, Vector& b)
    {
        // Update matrix entries if it has not been initialized yet
        const bool firstcall = (matrix_ == nullptr);
        if (firstcall) {
            // Model will not change the matrix object. Hence simply store a pointer to the original one with a deleter
            // that does nothing.
            // OBS: We need to be able to scale the linear system, hence const_cast
            matrix_ = const_cast<Matrix*>(&M);
        }
        else {
            // Pointers should not change; throw if the case
            if (&M != matrix_) {
                OPM_THROW(std::logic_error, "TPSA: Matrix objects are expected to be reused when reassembling!");
            }

        }

        // Set right-hand side vector
        rhs_ = &b;

        // Zero out the overlapping cells in matrix (not in ilu0 case)
        std::string type = prm_.template get<std::string>("preconditioner.type", "paroverilu0");
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);
        if (isParallel() && type != "paroverilu0") {
            detail::makeOverlapRowsInvalid(getMatrix(), overlapRows_);
        }
    }

    /*!
    * \brief Prepare matix and rhs vector for linear solve
    *
    * \param M System matrix
    * \param b Right-hand side vector
    *
    * \note No setResidual() or setMatrix() functions in this class. Must be handled here!
    */
    void prepare(const Matrix& M, Vector& b)
    {
        try {
            initPrepare(M, b);

            prepareFlexibleSolver();
        }
        OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR
            ("TPSA: Failure likely due to a faulty linear solver JSON specification. "
             "Check for errors related to missing nodes.");
    }

    /*!
    * \brief Prepare matix and rhs vector for linear solve
    *
    * \param M System matrix
    * \param b Right-hand side vector
    *
    * \note No setResidual() or setMatrix() functions in this class. Must be handled here!
    */
    void prepare(const SparseMatrixAdapter& M, Vector& b)
    {
        prepare(M.istlMatrix(), b);
    }

    /*!
    * \brief Prepare linear solver
    */
    void prepareFlexibleSolver()
    {
        // Create solver or just update preconditioner
        if (!flexibleSolver_.solver_) {
            // Dummy weights calculator
            // TODO: TPSA specific weights calculation for AMG preconditioner
            std::function<Vector()> weightCalculator;

            // Create FlexibleSolver
            flexibleSolver_.create(getMatrix(),
                                   isParallel(),
                                   prm_,
                                   /*pressureIndex=*/0,
                                   weightCalculator,
                                   /*forceSerial_=*/false,
                                   comm_.get());
        }
        else {
            // Update preconditioner
            flexibleSolver_.pre_->update();
        }
    }

    /*!
    * \brief Solve the linear system and store result in the input Vector x
    *
    * \param x Linear system solution will be stored here
    */
    bool solve(Vector& x)
    {
        // Increase solver count
        ++solveCount_;

        // Write system matrix if verbosity level is high
        const int verbosity = prm_.get("verbosity", 0);
        if (verbosity > 10) {
            // simulator_ is only used to get names
            Helper::writeSystem(simulator_,
                                getMatrix(),
                                *rhs_,
                                comm_.get());
        }

        // Solve linear system
        Dune::InverseOperatorResult result;
        assert(flexibleSolver_.solver_);
        flexibleSolver_.solver_->apply(x, *rhs_, result);

        // Store no. linear iterations
        iterations_ = result.iterations;

        // Return result for convergence check (boolean)
        return checkConvergence(result);
    }

    /*!
    * \brief Reset number of solver calls to zero
    */
    void resetSolveCount() {
        solveCount_ = 0;
    }

    // ///
    // Public get functions
    // ///
    /*!
    * \brief Copy right-hand side vector (rhs_) to incomming vector (b)
    *
    * \param b Vector to copy rhs_ to
    */
    void getResidual(Vector& b) const
    {
        b = *rhs_;
    }

    /*!
    * \brief Get number of solver calls
    */
    int getSolveCount() const
    {
        return solveCount_;
    }

    /*!
    * \brief Get number of linear solver iterations
    */
    int iterations () const
    {
        return iterations_;
    }

protected:
    /*!
    * \brief Check for parallel session
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
    */
    bool checkConvergence(const Dune::InverseOperatorResult& result) const
    {
        // Check relaxed linear solver tolerance
        if (!result.converged && result.reduction < parameters_.relaxed_linear_solver_reduction_) {
            std::string msg = fmt::format("Full linear solver tolerance not achieved. The reduction is {} "
                                            "after {} iterations.");
            OpmLog::warning(msg);
            return true;
        }

        // If we have failed and we don't ignore failures, throw error
        if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
            const std::string msg("Convergence failure for linear solver.");
            OPM_THROW_NOLOG(NumericalProblem, msg);
        }

        // Return convergence bool from linear solver result
        return result.converged;
    }

    // ///
    // Protected get functions
    // ///
    /*!
    * \brief Get reference to system matrix object
    */
    Matrix& getMatrix()
    {
        if (!matrix_) {
            OPM_THROW(std::runtime_error, "TPSA: System matrix \"M\" not defined!");
        }
        return *matrix_;
    }

    /*!
    * \brief Get reference to system matrix object
    */
    const Matrix& getMatrix() const
    {
        if (!matrix_) {
            OPM_THROW(std::runtime_error, "TPSA: System matrix \"M\" not defined!");
        }
        return *matrix_;
    }

    const Simulator& simulator_;
    std::any parallelInformation_;
    int solveCount_;
    int iterations_;

    detail::FlexibleSolverInfo<Matrix, Vector, CommunicationType> flexibleSolver_;
    TpsaLinearSolverParameters parameters_;
    PropertyTree prm_;

    Matrix* matrix_;
    Vector* rhs_;

    std::shared_ptr<CommunicationType> comm_;
    std::vector<int> overlapRows_;
};

}  // namespace Opm

#endif