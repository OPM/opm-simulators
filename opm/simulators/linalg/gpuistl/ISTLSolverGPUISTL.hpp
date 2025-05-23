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
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>

#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrix.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#endif

#include <opm/simulators/linalg/gpuistl/detail/FlexibleSolverWrapper.hpp>

namespace Opm::gpuistl
{


template <class TypeTag>
class ISTLSolverGPUISTL : public AbstractISTLSolver<TypeTag>
{
public:
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using real_type = typename Vector::field_type;

    using GPUMatrix = Opm::gpuistl::GpuSparseMatrix<real_type>;
    using GPUVector = Opm::gpuistl::GpuVector<real_type>;
    using GPUVectorInt = Opm::gpuistl::GpuVector<int>;

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
        , m_simulator(simulator)
    {
        m_parameters.init(simulator.vanguard().eclState().getSimulationConfig().useCPR());
        m_propertyTree = setupPropertyTree(m_parameters,
                                           Parameters::IsSet<Parameters::LinearSolverMaxIter>(),
                                           Parameters::IsSet<Parameters::LinearSolverReduction>());

        useWellConn_ = Parameters::Get<Parameters::MatrixAddWellContributions>();

        if (!useWellConn_) {
            OPM_THROW(std::logic_error, "Well operators are currently not supported for the GPU backend. "
            "Use --matrix-add-well-contributions=true to add well contributions to the matrix instead.");
        }
    }

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    explicit ISTLSolverGPUISTL(const Simulator& simulator)
        : ISTLSolverGPUISTL(simulator, FlowLinearSolverParameters(), false)
    {
    }



    void eraseMatrix() override
    {
        // Nothing, this is the same as the ISTLSolver
    }

    void setActiveSolver(int num) override
    {
        if (num != 0) {
            OPM_THROW(std::logic_error, "Only one solver available for the GPU backend.");
        }
    }

    int numAvailableSolvers() const override
    {
        return 1;
    }

    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        prepare(M.istlMatrix(), b);
    }

    void prepare(const Matrix& M, Vector& b) override
    {
        try {
            updateMatrix(M);
            updateRhs(b);
        }
        OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR("This is likely due to a faulty linear solver JSON specification. "
                                                "Check for errors related to missing nodes.");
    }

    void setResidual(Vector&) override
    {
        // Should be handled in prepare() instead.
    }

    void getResidual(Vector& b) const override
    {
        if (!m_rhs) {
            OPM_THROW(std::runtime_error, "m_rhs not initialized, prepare(matrix, rhs); needs to be called");
        }
        m_rhs->copyToHost(b);
    }

    void setMatrix(const SparseMatrixAdapter&) override
    {
        // Should be handled in prepare() instead.
    }

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
            OPM_THROW(std::runtime_error, "m_gpuFlexibleSolver not initialized, prepare(matrix, rhs); needs to be called");
        }
    
        if (!m_x) {
            m_x = std::make_unique<GPUVector>(x);
        } else {
            m_x->copyFromHost(x);
        }
        m_gpuSolver->apply(*m_x, *m_rhs, result);
    
        m_x->copyToHost(x);

        ++m_solveCount;


        m_lastSeenIterations = result.iterations;
        if (!result.converged) {
            if (result.reduction < m_parameters.relaxed_linear_solver_reduction_) {
                std::stringstream ss;
                ss << "Full linear solver tolerance not achieved. The reduction is:" << result.reduction << " after "
                   << result.iterations << " iterations ";
                OpmLog::warning(ss.str());
                return true;
            }
        }
        // Check for failure of linear solver.
        if (!m_parameters.ignoreConvergenceFailure_ && !result.converged) {
            const std::string msg("Convergence failure for linear solver.");
            OPM_THROW_NOLOG(NumericalProblem, msg);
        }

        return result.converged;
    }

    int iterations() const override
    {
        return m_lastSeenIterations;
    }

    const CommunicationType* comm() const override
    {
        // TODO: Implement this if needed
        return nullptr;
    }

    int getSolveCount() const override
    {
        return m_solveCount;
    }

private:
    // Weights to make approximate pressure equations.
    std::function<GPUVector&()> getWeightsCalculator()
    {
        std::function<GPUVector&()> weightsCalculator;

        using namespace std::string_literals;

        auto preconditionerType = m_propertyTree.get("preconditioner.type"s, "cpr"s);
        // Make the preconditioner type lowercase for internal canonical representation
        std::transform(preconditionerType.begin(), preconditionerType.end(), preconditionerType.begin(), ::tolower);
        if (preconditionerType == "cpr" || preconditionerType == "cprt"
            || preconditionerType == "cprw" || preconditionerType == "cprwt") {
            const bool transpose = preconditionerType == "cprt" || preconditionerType == "cprwt";
            const auto weightsType = m_propertyTree.get("preconditioner.weight_type"s, "quasiimpes"s);
            if (weightsType == "quasiimpes") {
                m_weights = std::make_unique<GPUVector>(m_matrix->N() * m_matrix->blockSize());
                // Pre-compute diagonal indices once when setting up the calculator
                auto diagonalIndices = Opm::Amg::precomputeDiagonalIndices(*m_matrix);
                m_diagonalIndices = std::make_unique<GPUVectorInt>(diagonalIndices);

                weightsCalculator = [this, transpose]() -> GPUVector& {
                    // GPU implementation for quasiimpes weights
                    Amg::getQuasiImpesWeights(*m_matrix, pressureIndex, transpose, *m_weights, *m_diagonalIndices);
                    return *m_weights;
                };
            } else if (weightsType == "trueimpes") {
                // Create CPU vector for the weights and initialize GPU vector
                m_cpuWeights.resize(m_matrix->N());
                m_weights = std::make_unique<GPUVector>(m_cpuWeights);

                // CPU implementation wrapped for GPU
                weightsCalculator = [this]() -> GPUVector& {

                    // Use the CPU implementation to calculate the weights
                    ElementContext elemCtx(m_simulator);
                    Amg::getTrueImpesWeights(pressureIndex, m_cpuWeights,
                                             m_simulator.vanguard().gridView(),
                                             elemCtx, m_simulator.model(),
                                             ThreadManager::threadId());
                    // Copy CPU vector to GPU vector
                    m_weights->copyFromHost(m_cpuWeights);
                    return *m_weights;
                };
            } else if (weightsType == "trueimpesanalytic") {
                // Create CPU vector for the weights and initialize GPU vector
                m_cpuWeights.resize(m_matrix->N());
                m_weights = std::make_unique<GPUVector>(m_cpuWeights);

                // CPU implementation wrapped for GPU
                weightsCalculator = [this]() -> GPUVector& {

                    // Use the CPU implementation to calculate the weights
                    ElementContext elemCtx(m_simulator);
                    Amg::getTrueImpesWeightsAnalytic(pressureIndex, m_cpuWeights,
                                                     m_simulator.vanguard().gridView(),
                                                     elemCtx, m_simulator.model(),
                                                     ThreadManager::threadId());
                    // Copy CPU vector to GPU vector
                    m_weights->copyFromHost(m_cpuWeights);
                    return *m_weights;
                };
            } else {
                OPM_THROW(std::invalid_argument,
                          "Weights type " + weightsType +
                          " not implemented for cpr."
                          " Please use quasiimpes, trueimpes or trueimpesanalytic.");
            }
        }
        return weightsCalculator;
    }

    void updateMatrix(const Matrix& M)
    {
        if (!m_matrix) {
            m_matrix.reset(new auto(GPUMatrix::fromMatrix(M)));
            std::function<GPUVector&()> weightsCalculator = getWeightsCalculator();
            const bool parallel = false;
            m_gpuSolver = std::make_unique<SolverType>(
                *m_matrix, parallel, m_propertyTree, pressureIndex, weightsCalculator, m_forceSerial, nullptr);
        } else {
            m_matrix->updateNonzeroValues(M);
        }
    
        m_gpuSolver->update();
    }

    void updateRhs(const Vector& b) {
        if (!m_rhs) {
            m_rhs = std::make_unique<GPUVector>(b);
        } else {
            m_rhs->copyFromHost(b);
        }
    }

    FlowLinearSolverParameters m_parameters;
    const bool m_forceSerial;
    const Simulator& m_simulator;
    PropertyTree m_propertyTree;
    bool useWellConn_;

    int m_lastSeenIterations = 0;
    int m_solveCount = 0;

    std::unique_ptr<GPUMatrix> m_matrix;

    std::unique_ptr<SolverType> m_gpuSolver;

    std::unique_ptr<GPUVector> m_rhs;
    std::unique_ptr<GPUVector> m_x;

    Vector m_cpuWeights;
    std::unique_ptr<GPUVector> m_weights;
    std::unique_ptr<GPUVectorInt> m_diagonalIndices;

};
} // namespace Opm::gpuistl

#endif
