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

#ifndef OPM_ISTLSOLVERGPUISTL_HEADER_INCLUDED
#define OPM_ISTLSOLVERGPUISTL_HEADER_INCLUDED

#include <dune/istl/operators.hh>
#include <memory>
#include <opm/simulators/linalg/AbstractISTLSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>

#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#endif

#include <opm/simulators/linalg/gpuistl/detail/FlexibleSolverWrapper.hpp>
#include <opm/simulators/linalg/printlinearsolverparameter.hpp>

namespace Opm::gpuistl
{

/**
* \brief ISTL solver for GPU using the GPU ISTL backend.
*
* This class implements the AbstractISTLSolver interface and provides
* methods to prepare the solver, set and get residuals, solve the system,
* and manage communication.
*
* \tparam TypeTag The type tag for the properties used in this solver.
*
* \note This solver takes CPU matrices and vectors, but uses GPU
*       matrices and vectors internally for computations.
*/
template <class TypeTag>
class ISTLSolverGPUISTL : public AbstractISTLSolver<TypeTag>
{
public:
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;

    using real_type = typename Vector::field_type;

    using GPUMatrix = Opm::gpuistl::GpuSparseMatrix<real_type>;
    using GPUVector = Opm::gpuistl::GpuVector<real_type>;

    constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;


#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    using SolverType = Opm::gpuistl::detail::FlexibleSolverWrapper<GPUMatrix, GPUVector, CommunicationType>;

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    /// \param[in] parameters  Explicit parameters for solver setup, do not
    ///                        read them from command line parameters.
    /// \param[in] forceSerial If true, will set up a serial linear solver only,
    ///                        local to the current rank, instead of creating a
    ///                        parallel (MPI distributed) linear solver.
    ISTLSolverGPUISTL(const Simulator& simulator,
                      const FlowLinearSolverParameters& parameters,
                      bool forceSerial = false)
        : m_parameters(parameters)
        , m_forceSerial(forceSerial)
    {
        m_parameters.init(simulator.vanguard().eclState().getSimulationConfig().useCPR());
        m_propertyTree = setupPropertyTree(m_parameters,
                                           Parameters::IsSet<Parameters::LinearSolverMaxIter>(),
                                           Parameters::IsSet<Parameters::LinearSolverReduction>());

        Opm::detail::printLinearSolverParameters(m_parameters, m_propertyTree, simulator.gridView().comm());
    }

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    explicit ISTLSolverGPUISTL(const Simulator& simulator)
        : ISTLSolverGPUISTL(simulator, FlowLinearSolverParameters(), false)
    {
    }

    /**
     * \copydoc AbstractISTLSolver::eraseMatrix
     * 
     * \note This method will not do anything.
     */
    void eraseMatrix() override
    {
        // Nothing, this is the same as the ISTLSolver
    }

    /**
     * \copydoc AbstractISTLSolver::setActiveSolver
     * 
     * \note This method will throw an exception if the solver number is not 0,
     *       as only one solver is available for the GPU backend.
     */
    void setActiveSolver(int num) override
    {
        if (num != 0) {
            OPM_THROW(std::logic_error, "Only one solver available for the GPU backend.");
        }
    }

    /**
     * \copydoc AbstractISTLSolver::numAvailableSolvers
     * 
     * \note This method always returns 1, as only one solver is available for the GPU backend.
     */
    int numAvailableSolvers() const override
    {
        return 1;
    }

    /**
     * \brief Prepare the solver with the given matrix and right-hand side vector.
     *
     * This method initializes the solver with the provided matrix and right-hand side vector.
     * It updates the internal GPU matrix and right-hand side vector.
     *
     * \param M The matrix to be used in the solver.
     * \param b The right-hand side vector.
     */
    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        prepare(M.istlMatrix(), b);
    }

    /**
     * \brief Prepare the solver with the given matrix and right-hand side vector.
     *
     * This method initializes the solver with the provided matrix and right-hand side vector.
     * It updates the internal GPU matrix and right-hand side vector.
     *
     * \param M The matrix to be used in the solver.
     * \param b The right-hand side vector.
     */
    void prepare(const Matrix& M, Vector& b) override
    {
        try {
            updateMatrix(M);
            updateRhs(b);
        }
        OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR("This is likely due to a faulty linear solver JSON specification. "
                                                "Check for errors related to missing nodes.");
    }

    /** 
     * \copydoc AbstractISTLSolver::setResidual
     *
     * \note Unused in this implementation.
     */
    void setResidual(Vector&) override
    {
        // Should be handled in prepare() instead.
    }

    /**
     * \brief Get the residual vector.
     *
     * This method retrieves the current residual vector from the solver.
     * It copies the data from the internal GPU vector to the provided vector.
     *
     * \param b The vector to store the residual.
     */
    void getResidual(Vector& b) const override
    {
        if (!m_rhs) {
            OPM_THROW(std::runtime_error, "m_rhs not initialized, prepare(matrix, rhs); needs to be called");
        }
        m_rhs->copyToHost(b);
    }

    /**
     * \copydoc AbstractISTLSolver::setMatrix
     * 
     * \note This method does not set the matrix directly, as it should be handled in prepare().
     */
    void setMatrix(const SparseMatrixAdapter&) override
    {
        // Should be handled in prepare() instead.
    }

    /**
     * \brief Solve the system of linear equations Ax = b.
     *
     * This method solves the linear system represented by the matrix A and the right-hand side vector b,
     * storing the solution in vector x.
     *
     * \param x The vector to store the solution.
     * \return true if the solver converged, false otherwise.
     *
     * Before this function is called, prepare() should have been called with a valid matrix and right-hand side vector.
     */
    bool solve(Vector& x) override
    {
        // TODO: Write matrix to disk if needed
        Dune::InverseOperatorResult result;
        if (!m_matrix) {
            OPM_THROW(std::runtime_error, "m_matrix not initialized, prepare(matrix, rhs); needs to be called");
        }
        if (!m_rhs) {
            OPM_THROW(std::runtime_error, "m_rhs not initialized, prepare(matrix, rhs); needs to be called");
        }
        if (!m_gpuSolver) {
            OPM_THROW(std::runtime_error,
                      "m_gpuFlexibleSolver not initialized, prepare(matrix, rhs); needs to be called");
        }

        if (!m_x) {
            m_x = std::make_unique<GpuVector<real_type>>(x);
        } else {
            m_x->copyFromHost(x);
        }
        m_gpuSolver->apply(*m_x, *m_rhs, result);

        m_x->copyToHost(x);

        ++m_solveCount;

        m_lastSeenIterations = result.iterations;
        return checkConvergence(result);
    }

    /**
     * \copydoc AbstractISTLSolver::post
     *
     * Returns the actual number of iterations used in the last solve.
     */
    int iterations() const override
    {
        return m_lastSeenIterations;
    }

    /**
     * \copydoc AbstractISTLSolver::comm
     */
    const CommunicationType* comm() const override
    {
        // TODO: Implement this if needed
        OPM_THROW(std::logic_error, "Communication not implemented for GPU ISTL solver.");
    }

    /**
     * \copydoc AbstractISTLSolver::getSolveCount
     */
    int getSolveCount() const override
    {
        return m_solveCount;
    }

private:
    bool checkConvergence(const Dune::InverseOperatorResult& result) const
    {
        return AbstractISTLSolver<TypeTag>::checkConvergence(result, m_parameters);
    }

    void updateMatrix(const Matrix& M)
    {
        if (!m_matrix) {

            m_matrix.reset(new auto(GPUMatrix::fromMatrix(M)));
            std::function<GPUVector()> weightsCalculator = {};
            const bool parallel = false;
            m_gpuSolver = std::make_unique<SolverType>(
                *m_matrix, parallel, m_propertyTree, pressureIndex, weightsCalculator, m_forceSerial, nullptr);
        } else {
            m_matrix->updateNonzeroValues(M);
            m_gpuSolver->update();
        }
    }

    void updateRhs(const Vector& b)
    {
        if (!m_rhs) {
            m_rhs = std::make_unique<GPUVector>(b);
        } else {
            m_rhs->copyFromHost(b);
        }
    }

    FlowLinearSolverParameters m_parameters;
    PropertyTree m_propertyTree;

    int m_lastSeenIterations = 0;
    int m_solveCount = 0;

    std::unique_ptr<GPUMatrix> m_matrix;

    std::unique_ptr<SolverType> m_gpuSolver;

    std::unique_ptr<GPUVector> m_rhs;
    std::unique_ptr<GPUVector> m_x;

    const bool m_forceSerial;
};
} // namespace Opm::gpuistl

#endif
