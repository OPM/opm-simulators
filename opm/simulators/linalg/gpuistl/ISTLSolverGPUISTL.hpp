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

#include <opm/simulators/linalg/gpuistl/detail/ISTLSolverGPUISTLImplementation.hpp>

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
    constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;


#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif




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
    {
        m_parameters.init(simulator.vanguard().eclState().getSimulationConfig().useCPR());
        m_propertyTree = setupPropertyTree(m_parameters,
                                           Parameters::IsSet<Parameters::LinearSolverMaxIter>(),
                                           Parameters::IsSet<Parameters::LinearSolverReduction>());
        m_implementation = std::make_unique<detail::ISTLSolverGPUISTLImplementation<Matrix, Vector, CommunicationType>>(
            m_parameters, m_propertyTree, forceSerial, pressureIndex);
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
            m_implementation->prepare(M, b);
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
        m_implementation->getResidual(b);
    }

    void setMatrix(const SparseMatrixAdapter&) override
    {
        // Should be handled in prepare() instead.
    }

    bool solve(Vector& x) override
    {
        // TODO: Write matrix to disk if needed
        Dune::InverseOperatorResult result;
        m_implementation->apply(x, result);
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
    FlowLinearSolverParameters m_parameters;
    PropertyTree m_propertyTree;

    std::unique_ptr<detail::ISTLSolverGPUISTLImplementation<Matrix, Vector, CommunicationType>> m_implementation;
    int m_lastSeenIterations = 0;
    int m_solveCount = 0;
};
} // namespace Opm::gpuistl

#endif
