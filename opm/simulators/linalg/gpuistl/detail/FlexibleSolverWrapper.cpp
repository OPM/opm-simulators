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

#include "config.h"

#include <functional>
#include <opm/simulators/linalg/gpuistl/detail/FlexibleSolverWrapper.hpp>

#include <dune/common/parallel/communication.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuOwnerOverlapCopy.hpp>

namespace Opm::gpuistl::detail
{

namespace
{
    template <class Matrix, class Vector, class Comm>
    std::tuple<typename FlexibleSolverWrapper<Matrix, Vector, Comm>::AbstractOperatorPtrType,
               typename FlexibleSolverWrapper<Matrix, Vector, Comm>::AbstractSolverPtrType,
               std::reference_wrapper<typename FlexibleSolverWrapper<Matrix, Vector, Comm>::AbstractPreconditionerType>,
               std::shared_ptr<typename FlexibleSolverWrapper<Matrix, Vector, Comm>::GpuCommunicationType>>
    createOperatorAndSolver(const Matrix& matrix,
                            [[maybe_unused]] bool parallel,
                            const PropertyTree& prm,
                            std::size_t pressureIndex,
                            const std::function<Vector()>& weightCalculator,
                            [[maybe_unused]] bool forceSerial,
                            [[maybe_unused]] const Comm* comm)
    {
        using Wrapper = FlexibleSolverWrapper<Matrix, Vector, Comm>;
        if (!parallel || forceSerial) {
            // For now only matrix adapter is supported
            using OperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
            using SolverType = Dune::FlexibleSolver<OperatorType>;

            auto operatorPtr = std::make_unique<OperatorType>(matrix);

            auto solverPtr = std::make_unique<SolverType>(*operatorPtr, prm, weightCalculator, pressureIndex);
            auto preconditioner = std::ref(solverPtr->preconditioner());

            return std::make_tuple(std::move(operatorPtr), std::move(solverPtr), preconditioner, 
                                   std::shared_ptr<typename Wrapper::GpuCommunicationType>{});
        }
#if HAVE_MPI
        else {
            using real_type = typename Matrix::field_type;
            using return_type = std::tuple<typename Wrapper::AbstractOperatorPtrType,
            typename Wrapper::AbstractSolverPtrType,
            std::reference_wrapper<typename Wrapper::AbstractPreconditionerType>,
            std::shared_ptr<typename Wrapper::GpuCommunicationType>>;

            // We need the block size at compile time to instantiate the correct types
            // hence we need to dispatch on the block size
            return matrix.dispatchOnBlocksize([&](auto blockSizeVal) -> return_type {
                // Get the block size from the decltype of the blockSizeVal, 
                // making it a compile time constant
                constexpr int block_size = decltype(blockSizeVal)::value;

                // We are running in parallel (MPI) so we need a communication type, this
                // time for GPU
                using CudaCommunication = GpuOwnerOverlapCopy<real_type, Comm>;
                using SchwarzOperator
                    = Dune::OverlappingSchwarzOperator<GpuSparseMatrix<real_type>, Vector, Vector, CudaCommunication>;
  
                using SolverType = Dune::FlexibleSolver<SchwarzOperator>;

                // Create the communication object that will handle the GPU and MPI communication
                auto cudaCommunication = makeGpuOwnerOverlapCopy<real_type, block_size, Comm>(*comm);
                // Create the operator that will (through the communication object) handle the
                // GPU and MPI communication
                auto operatorPtr = std::make_unique<SchwarzOperator>(matrix, *cudaCommunication);
                auto solverPtr = std::make_unique<SolverType>(*operatorPtr, *cudaCommunication, prm, weightCalculator, pressureIndex);
                auto preconditioner = std::ref(solverPtr->preconditioner());

                return std::make_tuple(
                    std::move(operatorPtr), std::move(solverPtr), preconditioner,
                    std::static_pointer_cast<typename Wrapper::GpuCommunicationType>(cudaCommunication));
            });
        }
#else
        else {
            OPM_THROW(std::logic_error, "Parallel GPU ISTL solver requires MPI support.");
        }
#endif
    }

} // namespace

template <class Matrix, class Vector, class Comm>
FlexibleSolverWrapper<Matrix, Vector, Comm>::FlexibleSolverWrapper(const Matrix& matrix,
                                                                   bool parallel,
                                                                   const PropertyTree& prm,
                                                                   std::size_t pressureIndex,
                                                                   const std::function<Vector()>& weightCalculator,
                                                                   bool forceSerial,
                                                                   const Comm* comm)
    : FlexibleSolverWrapper(
          createOperatorAndSolver(matrix, parallel, prm, pressureIndex, weightCalculator, forceSerial, comm))
{
}

template <class Matrix, class Vector, class Comm>
FlexibleSolverWrapper<Matrix, Vector, Comm>::FlexibleSolverWrapper(
    std::tuple<AbstractOperatorPtrType, AbstractSolverPtrType, std::reference_wrapper<AbstractPreconditionerType>, std::shared_ptr<GpuCommunicationType>>&&
        solverTuple)
    : m_operator(std::move(std::get<0>(solverTuple)))
    , m_solver(std::move(std::get<1>(solverTuple)))
    , m_preconditioner(std::get<2>(solverTuple))
    , m_gpuCommunication(std::get<3>(solverTuple))
{
}

template <class Matrix, class Vector, class Comm>
void
FlexibleSolverWrapper<Matrix, Vector, Comm>::apply(Vector& x, Vector& y, Dune::InverseOperatorResult& result)
{
    m_solver->apply(x, y, result);
}

template <class Matrix, class Vector, class Comm>
void
FlexibleSolverWrapper<Matrix, Vector, Comm>::update()
{
    m_preconditioner.update();
}


} // namespace Opm::gpuistl::detail

#if HAVE_MPI
using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
using CommunicationType = Dune::Communication<int>;
#endif

#define INSTANTIATE_FLEXIBLE_SOLVER_WRAPPER(real_type)                                                                 \
    template class ::Opm::gpuistl::detail::FlexibleSolverWrapper<::Opm::gpuistl::GpuSparseMatrix<real_type>,           \
                                                                 ::Opm::gpuistl::GpuVector<real_type>,                 \
                                                                 CommunicationType>

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_FLEXIBLE_SOLVER_WRAPPER(float);
#endif
INSTANTIATE_FLEXIBLE_SOLVER_WRAPPER(double);
