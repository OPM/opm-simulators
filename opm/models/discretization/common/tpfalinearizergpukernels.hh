/*
  Copyright 2026 Equinor ASA
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

#ifndef TPFA_LINEARIZER_GPUKERNELS_HH
#define TPFA_LINEARIZER_GPUKERNELS_HH

/*
    This file contains the GPU kernels that handles the GPU parallelization of the linearization.
    This is effectively just the GPU alternative to what is a for-loop on the CPU
*/

namespace Opm
{

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
template <class TpfaLinearizerType,
          class ModelClass,
          class LocalResidualT,
          class DiagPtrType,
          class ScalarType,
          class DomainType,
          class NeighborSparseTable,
          class ResidualType,
          class ProblemT>
__global__ __launch_bounds__(256) void kernel_linearize(const unsigned int numCells,
                                                        const DomainType domain,
                                                        const NeighborSparseTable neighborInfo,
                                                        DiagPtrType diagMatAddress,
                                                        ResidualType residual,
                                                        ModelClass model,
                                                        ScalarType dt,
                                                        ProblemT problem)
{
    const unsigned int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if (ii < numCells) {
        // velocityInfo is unused on GPU (guarded by if constexpr (!useGPU)),
        // but the parameter is T& so we need an lvalue to bind to.
        std::nullptr_t dummyVelocityInfo = nullptr;
        TpfaLinearizerType::template linearize_cell<true, LocalResidualT>(ii,
                                                                          domain,
                                                                          neighborInfo,
                                                                          diagMatAddress,
                                                                          residual,
                                                                          model,
                                                                          dummyVelocityInfo,
                                                                          dt,
                                                                          problem);
    }
}

template <class TpfaLinearizerType,
          class IntensiveQuantities,
          class ModelT,
          class LocalResidualT,
          class DiagPtrType,
          class ResidualT,
          class BoundaryInfo,
          class ProblemT>
__global__ void
linearize_bc_threadsafe(DiagPtrType diagMatAddress,
                        ResidualT residual,
                        const BoundaryInfo boundaryInfo,
                        ModelT model,
                        ProblemT gpuProblem)
{
    const unsigned int ii = blockIdx.x * blockDim.x + threadIdx.x;
    if (ii < boundaryInfo.size()) {
        TpfaLinearizerType::template linearize_bc_threadsafe_single_cell<IntensiveQuantities,
                                                                         LocalResidualT>(
            diagMatAddress, residual, boundaryInfo[ii], model, gpuProblem);
    }
}

#endif

} // namespace Opm

#endif // TPFA_LINEARIZER_GPUKERNELS_HH
