/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
#define OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <fmt/format.h>

#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/models/parallel/threadmanager.hpp>
#include <opm/simulators/linalg/TrueImpesWeights.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/detail/cpr_amg_operations.hpp>
#else
#include <opm/simulators/linalg/gpuistl/detail/cpr_amg_operations.hpp>
#endif
#endif


namespace Opm
{

namespace Details
{
    template <class DenseMatrix>
    DenseMatrix transposeDenseMatrix(const DenseMatrix& M)
    {
        DenseMatrix tmp;
        for (int i = 0; i < M.rows; ++i)
            for (int j = 0; j < M.cols; ++j)
                tmp[j][i] = M[i][j];

        return tmp;
    }
} // namespace Details

namespace Amg
{
    template <class Matrix, class Vector>
    void getQuasiImpesWeights(const Matrix& matrix,
                              const int pressureVarIndex,
                              const bool transpose,
                              Vector& weights,
                              [[maybe_unused]] bool enable_thread_parallel)
    {
        using VectorBlockType = typename Vector::block_type;
        using MatrixBlockType = typename Matrix::block_type;
        const Matrix& A = matrix;

        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;

        // Declare variables outside the loop to avoid repetitive allocation
        MatrixBlockType diag_block;
        VectorBlockType bweights;
        MatrixBlockType diag_block_transpose;

        // Use OpenMP to parallelize over matrix rows (runtime controlled via if clause)
#ifdef _OPENMP
#pragma omp parallel for private(diag_block, bweights, diag_block_transpose) if(enable_thread_parallel)
#endif
        for (int row_idx = 0; row_idx < static_cast<int>(A.N()); ++row_idx) {
            diag_block = MatrixBlockType(0.0);
            // Find diagonal block for this row
            const auto row_it = A.begin() + row_idx;
            const auto endj = (*row_it).end();
            for (auto j = (*row_it).begin(); j != endj; ++j) {
                if (row_it.index() == j.index()) {
                    diag_block = (*j);
                    break;
                }
            }
            if (transpose) {
                diag_block.solve(bweights, rhs);
            } else {
                diag_block_transpose = Details::transposeDenseMatrix(diag_block);
                diag_block_transpose.solve(bweights, rhs);
            }

            const double abs_max =
                *std::ranges::max_element(bweights,
                                          [](double a, double b)
                                          { return std::fabs(a) < std::fabs(b); });
            bweights /= std::fabs(abs_max);
            weights[row_idx] = bweights;
        }
    }

    template <class Matrix, class Vector>
    Vector getQuasiImpesWeights(const Matrix& matrix,
                                const int pressureVarIndex,
                                const bool transpose,
                                bool enable_thread_parallel)
    {
        Vector weights(matrix.N());
        getQuasiImpesWeights(matrix, pressureVarIndex, transpose, weights, enable_thread_parallel);
        return weights;
    }

#if HAVE_CUDA
    template <typename T>
    std::vector<int> precomputeDiagonalIndices(const gpuistl::GpuSparseMatrixWrapper<T>& matrix) {
        std::vector<int> diagonalIndices(matrix.N(), -1);
        const auto rowIndices = matrix.getRowIndices().asStdVector();
        const auto colIndices = matrix.getColumnIndices().asStdVector();

        for (auto row = 0; row < Opm::gpuistl::detail::to_int(matrix.N()); ++row) {
            for (auto i = rowIndices[row]; i < rowIndices[row+1]; ++i) {
                if (colIndices[i] == row) {
                    diagonalIndices[row] = i;
                    break;
                }
            }
        }
        return diagonalIndices;
    }

    // GPU version that delegates to the GPU implementation
    template <typename T, bool transpose>
    void getQuasiImpesWeights(const gpuistl::GpuSparseMatrixWrapper<T>& matrix,
                             const int pressureVarIndex,
                             gpuistl::GpuVector<T>& weights,
                             const gpuistl::GpuVector<int>& diagonalIndices)
    {
        gpuistl::detail::getQuasiImpesWeights<T, transpose>(matrix, pressureVarIndex, weights, diagonalIndices);
    }

    template <typename T, bool transpose>
    gpuistl::GpuVector<T> getQuasiImpesWeights(const gpuistl::GpuSparseMatrixWrapper<T>& matrix,
                                              const int pressureVarIndex,
                                              const gpuistl::GpuVector<int>& diagonalIndices)
    {
        gpuistl::GpuVector<T> weights(matrix.N() * matrix.blockSize());
        getQuasiImpesWeights<T, transpose>(matrix, pressureVarIndex, weights, diagonalIndices);
        return weights;
    }
#endif

    template<class Vector, class ElementContext, class Model, class ElementChunksType>
    void getTrueImpesWeights(int pressureVarIndex, Vector& weights,
                             const ElementContext& elemCtx,
                             const Model& model,
                             const ElementChunksType& element_chunks,
                             [[maybe_unused]] bool enable_thread_parallel)
    {
        using VectorBlockType = typename Vector::block_type;
        using Matrix = typename std::decay_t<decltype(model.linearizer().jacobian())>;
        using MatrixBlockType = typename Matrix::MatrixBlock;
        constexpr int numEq = VectorBlockType::size();
        using Evaluation = typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>
            ::block_type;

        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;

        // Declare variables outside the loop to avoid repetitive allocation
        MatrixBlockType block;
        VectorBlockType bweights;
        MatrixBlockType block_transpose;
        Dune::FieldVector<Evaluation, numEq> storage;

        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(block, bweights, block_transpose, storage) if(enable_thread_parallel)
#endif
        for (std::size_t ci = 0; ci < element_chunks.size(); ++ci) {
            const auto& chunk = element_chunks[ci];
            const std::size_t thread_id = ThreadManager::threadId();
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                model.localLinearizer(thread_id).localResidual().computeStorage(storage, localElemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);

                auto extrusionFactor = localElemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
                auto scvVolume = localElemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;
                auto storage_scale = scvVolume / localElemCtx.simulator().timeStepSize();
                detail::formTrueImpesStorageJacobian<numEq>(storage, storage_scale,
                                                            pressureVarIndex, block_transpose);
                try {
                    block_transpose.solve(bweights, rhs);
                }
                catch (const Dune::FMatrixError&) {
                    // Rank-deficient storage derivatives: no combination of the
                    // mass balances has a storage term that depends on pressure
                    // alone, so there is no weight vector to compute here.
                    //
                    // Deliberately no fallback value.  bweights is indexed by
                    // equation while rhs is indexed by primary variable, so
                    // substituting rhs would put unit weight on whichever
                    // equation happens to share the pressure variable's index -
                    // an arbitrary choice that goes on to make the CPR pressure
                    // system itself singular.  Name the cell instead, so the
                    // cause can be found.
                    const auto globalDofIdx =
                        localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);

                    throw std::runtime_error {
                        fmt::format("Singular storage matrix when forming the CPR "
                                    "pressure weights for cell {} (Cartesian index "
                                    "{}).  The storage derivatives of that cell are "
                                    "rank deficient, so the pressure equation cannot "
                                    "be formed there.",
                                    globalDofIdx,
                                    localElemCtx.simulator().vanguard()
                                        .cartesianIndex(globalDofIdx))
                    };
                }

                const double abs_max =
                    *std::ranges::max_element(bweights,
                                              [](double a, double b)
                                              { return std::fabs(a) < std::fabs(b); });
                // probably a scaling which could give approximately total compressibility would be better
                bweights /=  std::fabs(abs_max); // given normal densities this scales weights to about 1.

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

    template <class Vector, class ElementContext, class Model, class ElementChunksType>
    void getTrueImpesWeightsAnalytic(int /*pressureVarIndex*/,
                                     Vector& weights,
                                     const ElementContext& elemCtx,
                                     const Model& model,
                                     const ElementChunksType& element_chunks,
                                     [[maybe_unused]] bool enable_thread_parallel)
    {
        using VectorBlockType = typename Vector::block_type;
        const auto& solution = model.solution(/*timeIdx*/ 0);

        // Use OpenMP to parallelize over element chunks (runtime controlled via if clause)
        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for if(enable_thread_parallel)
#endif
        for (std::size_t ci = 0; ci < element_chunks.size(); ++ci) {
            const auto& chunk = element_chunks[ci];

            // Each thread gets a unique copy of elemCtx
            ElementContext localElemCtx(elemCtx.simulator());
            // Extra equations (e.g. energy) have no analytic component weight.
            VectorBlockType bweights(0.0);

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = localElemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                detail::analyticWeights(intQuants.getFluidSystem(), fs, solution[index], bweights);
                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesAnalyticWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }
} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
