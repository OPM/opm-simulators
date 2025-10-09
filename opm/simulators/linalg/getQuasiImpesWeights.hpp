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

#include <dune/common/fvector.hh>

#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/models/parallel/threadmanager.hpp>
#include <algorithm>
#include <cmath>

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
    void getQuasiImpesWeights(const Matrix& matrix, const int pressureVarIndex, const bool transpose, Vector& weights, bool enable_thread_parallel)
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

            double abs_max = *std::max_element(
                bweights.begin(), bweights.end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); });
            bweights /= std::fabs(abs_max);
            weights[row_idx] = bweights;
        }
    }

    template <class Matrix, class Vector>
    Vector getQuasiImpesWeights(const Matrix& matrix, const int pressureVarIndex, const bool transpose, bool enable_thread_parallel)
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
                             bool enable_thread_parallel)
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
        for (const auto& chunk : element_chunks) {
            const std::size_t thread_id = ThreadManager::threadId();
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                model.localLinearizer(thread_id).localResidual().computeStorage(storage, localElemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);

                auto extrusionFactor = localElemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
                auto scvVolume = localElemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;
                auto storage_scale = scvVolume / localElemCtx.simulator().timeStepSize();
                const double pressure_scale = 50e5;

                // Build the transposed matrix directly to avoid separate transpose step
                for (int ii = 0; ii < numEq; ++ii) {
                    for (int jj = 0; jj < numEq; ++jj) {
                        block_transpose[jj][ii] = storage[ii].derivative(jj)/storage_scale;
                        if (jj == pressureVarIndex) {
                            block_transpose[jj][ii] *= pressure_scale;
                        }
                    }
                }
                block_transpose.solve(bweights, rhs);

                double abs_max = *std::max_element(
                    bweights.begin(), bweights.end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); });
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
                                     bool enable_thread_parallel)
    {
        // The sequential residual is a linear combination of the
        // mass balance residuals, with coefficients equal to (for
        // water, oil, gas):
        //    1/bw,
        //    (1/bo - rs/bg)/(1-rs*rv)
        //    (1/bg - rv/bo)/(1-rs*rv)
        // These coefficients must be applied for both the residual and
        // Jacobian.
        using FluidSystem = typename Model::FluidSystem;
        using LhsEval = double;

        using PrimaryVariables = typename Model::PrimaryVariables;
        using VectorBlockType = typename Vector::block_type;
        using Evaluation =
            typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>::block_type;
        using Toolbox = MathToolbox<Evaluation>;

        const auto& solution = model.solution(/*timeIdx*/ 0);
        VectorBlockType bweights;

        // Use OpenMP to parallelize over element chunks (runtime controlled via if clause)
        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(bweights) if(enable_thread_parallel)
#endif
        for (const auto& chunk : element_chunks) {

            // Each thread gets a unique copy of elemCtx
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = localElemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::waterPhaseIdx));
                    bweights[activeCompIdx]
                        = Toolbox::template decay<LhsEval>(1 / fs.invB(FluidSystem::waterPhaseIdx));
                }

                double denominator = 1.0;
                double rs = Toolbox::template decay<double>(fs.Rs());
                double rv = Toolbox::template decay<double>(fs.Rv());
                const auto& priVars = solution[index];
                if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                    rs = 0.0;
                }
                if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
                    rv = 0.0;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                    && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    denominator = Toolbox::template decay<LhsEval>(1 - rs * rv);
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::oilPhaseIdx));
                    bweights[activeCompIdx] = Toolbox::template decay<LhsEval>(
                        (1 / fs.invB(FluidSystem::oilPhaseIdx) - rs / fs.invB(FluidSystem::gasPhaseIdx))
                        / denominator);
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                    bweights[activeCompIdx] = Toolbox::template decay<LhsEval>(
                        (1 / fs.invB(FluidSystem::gasPhaseIdx) - rv / fs.invB(FluidSystem::oilPhaseIdx))
                        / denominator);
                }

                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesAnalyticWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }
} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
