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

#ifndef OPM_ISTLSOLVERRUNTIMEOPTIONPROXY_HEADER_INCLUDED
#define OPM_ISTLSOLVERRUNTIMEOPTIONPROXY_HEADER_INCLUDED

#include "opm/simulators/linalg/FlowLinearSolverParameters.hpp"
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/AbstractISTLSolver.hpp>
#include <opm/simulators/linalg/ISTLSolver.hpp>
#if COMPILE_GPU_BRIDGE
#include <opm/simulators/linalg/ISTLSolverGpuBridge.hpp>
#endif

#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/ISTLSolverGPUISTL.hpp>
#endif

namespace Opm
{

/**
 * \brief ISTLSolverRuntimeOptionProxy selects the appropriate ISTLSolver runtime based on CLI options
 */
template <class TypeTag>
class ISTLSolverRuntimeOptionProxy : public AbstractISTLSolver<TypeTag>
{
public:
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif


    static void registerParameters()
    {
        FlowLinearSolverParameters::registerParameters();
    }

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    /// \param[in] parameters  Explicit parameters for solver setup, do not
    ///                        read them from command line parameters.
    /// \param[in] forceSerial If true, will set up a serial linear solver only,
    ///                        local to the current rank, instead of creating a
    ///                        parallel (MPI distributed) linear solver.
    ISTLSolverRuntimeOptionProxy(const Simulator& simulator,
                                 const FlowLinearSolverParameters& parameters,
                                 bool forceSerial = false)
    {
        createSolver(simulator, parameters, forceSerial);
    }

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    explicit ISTLSolverRuntimeOptionProxy(const Simulator& simulator)
    {
        createSolver(simulator);
    }


    void eraseMatrix() override
    {
        istlSolver_->eraseMatrix();
    }

    void setActiveSolver(int num) override
    {
        istlSolver_->setActiveSolver(num);
    }

    int numAvailableSolvers() const override
    {
        return istlSolver_->numAvailableSolvers();
    }

    void prepare(const SparseMatrixAdapter& M, Vector& b) override
    {
        istlSolver_->prepare(M, b);
    }

    void prepare(const Matrix& M, Vector& b) override
    {
        istlSolver_->prepare(M, b);
    }

    void setResidual(Vector& b) override
    {
        istlSolver_->setResidual(b);
    }

    void getResidual(Vector& b) const override
    {
        istlSolver_->getResidual(b);
    }

    void setMatrix(const SparseMatrixAdapter& M) override
    {
        istlSolver_->setMatrix(M);
    }

    bool solve(Vector& x) override
    {
        return istlSolver_->solve(x);
    }

    int iterations() const override
    {
        return istlSolver_->iterations();
    }

    const CommunicationType* comm() const override
    {
        return istlSolver_->comm();
    }

    int getSolveCount() const override
    {
        return istlSolver_->getSolveCount();
    }

private:
    std::unique_ptr<AbstractISTLSolver<TypeTag>> istlSolver_;


    template <class... Args>
    void createSolver(const Simulator& simulator, Args&&... args)
    {
        const auto backend = Parameters::linearSolverAcceleratorTypeFromCLI();
        if (backend == Parameters::LinearSolverAcceleratorType::CPU) {
        // Note that for now we keep the old behavior of using the bridge solver if it is available.
#if COMPILE_GPU_BRIDGE
            istlSolver_ = std::make_unique<ISTLSolverGpuBridge<TypeTag>>(simulator, std::forward<Args>(args)...);
#else
            istlSolver_ = std::make_unique<ISTLSolver<TypeTag>>(simulator, std::forward<Args>(args)...);
#endif
        }
#if HAVE_CUDA
        else if (backend == Parameters::LinearSolverAcceleratorType::GPU) {
            istlSolver_ = std::make_unique<gpuistl::ISTLSolverGPUISTL<TypeTag>>(simulator, std::forward<Args>(args)...);
        }
#endif
        else {
            // If we reach here, it means the backend is not supported. This could be because we have added a third backend
            // that we need to handle. A user error would be handled in the linearSolverAcceleratorTypeFromString function called above.
            OPM_THROW(std::invalid_argument, fmt::format("Unknown backend: {}", Parameters::toString(backend)));
        }
    }
};
} // namespace Opm

#endif
