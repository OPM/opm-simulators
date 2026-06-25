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

    /// \brief Compute IMPES pressure-equation weights via AD storage Jacobian (derivative-based).
    ///
    /// For each cell, the storage Jacobian dM/dx (where M_i is the accumulation term
    /// for equation i and x_j is primary variable j) is assembled from the automatic
    /// differentiation (AD) derivatives of the storage residual.  The transposed system
    ///
    ///   (dM/dx)^T  w = e_pressure
    ///
    /// is then solved for the weight vector w.  This is the "numerical" or
    /// derivative-based IMPES approach: it works for any blackoil extension
    /// (dissolved gas in water, vaporised water, thermal, solvents, etc.)
    /// without requiring a closed-form analytic formula.
    ///
    /// This method is more expensive than the analytic alternatives
    /// (getTrueImpesWeightsAnalytic, getCoatsWeightsBlackoil) because it
    /// requires an element-context update and a small dense solve per cell.
    ///
    /// References:
    ///   - Wallis, J.R. (1983). "Incomplete Gaussian Elimination as a
    ///     Preconditioning for Generalized Conjugate Gradient Acceleration."
    ///     SPE-12265-MS. https://doi.org/10.2118/12265-MS
    ///   - Cao, H., Tchelepi, H.A., Wallis, J.R., Yardumian, H. (2005).
    ///     "Parallel Scalable Unstructured CPR-Type Linear Solver for
    ///     Reservoir Simulation." SPE-96809-MS.
    ///     https://doi.org/10.2118/96809-MS
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

    /// \brief Compute analytic IMPES/Coats pressure-equation weights for the blackoil model.
    ///
    /// Equivalent to getCoatsWeightsBlackoil() for the standard blackoil model
    /// (no dissolved gas in water, no vaporised water).  The weights are:
    ///
    ///   w_water = B_w
    ///   w_oil   = (B_o - R_s B_g) / (1 - R_s R_v)
    ///   w_gas   = (B_g - R_v B_o) / (1 - R_s R_v)
    ///
    /// where B_alpha = 1/invB_alpha is the formation-volume factor.
    /// These weights satisfy the volumetric-balance linear system
    ///
    ///   | invBw   0       0       | | w_w |   | 1 |
    ///   | 0       invBo   Rs*invBo| | w_o | = | 1 |
    ///   | 0       Rv*invBg invBg  | | w_g |   | 1 |
    ///
    /// so that the weighted sum of the component mass balances equals the
    /// total reservoir pore-volume balance  phi*(S_w + S_o + S_g).
    ///
    /// Undersaturation is handled by switching off the appropriate ratio:
    ///   - gas dissolved only (GasMeaning::Rs): set R_v = 0
    ///   - oil vaporised only (GasMeaning::Rv): set R_s = 0
    ///
    /// For the extended blackoil model with dissolved gas in water (R_sw)
    /// and vaporised water in gas (R_vw), use getCoatsWeightsBlackoil()
    /// instead, which adds the full cross-terms.
    ///
    /// References:
    ///   - Wallis, J.R. (1983). "Incomplete Gaussian Elimination as a
    ///     Preconditioning for Generalized Conjugate Gradient Acceleration."
    ///     SPE-12265-MS. https://doi.org/10.2118/12265-MS
    ///   - Coats, K.H. (1980). "An Equation of State Compositional Model."
    ///     SPE Journal, 20(5), 363-376. SPE-8284-PA.
    ///     https://doi.org/10.2118/8284-PA
    template <class Vector, class ElementContext, class Model, class ElementChunksType>
    void getTrueImpesWeightsAnalytic(int /*pressureVarIndex*/,
                                     Vector& weights,
                                     const ElementContext& elemCtx,
                                     const Model& model,
                                     const ElementChunksType& element_chunks,
                                     [[maybe_unused]] bool enable_thread_parallel)
    {
        // The sequential residual is a linear combination of the
        // mass balance residuals, with coefficients equal to (for
        // water, oil, gas) -- here bw means B_w (formation-volume factor):
        //    B_w                         = 1/invBw
        //    (B_o - Rs*B_g)/(1-Rs*Rv)   = (1/invBo - Rs/invBg)/(1-Rs*Rv)
        //    (B_g - Rv*B_o)/(1-Rs*Rv)   = (1/invBg - Rv/invBo)/(1-Rs*Rv)
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
                    rs = 0.0;// then oil saturation is zero an rs undetermined, in simulator it is saturated 
                }
                if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
                    rv = 0.0;// then gas saturation is zero but rv is undetermined, in simulator it is saturated
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
                // in the understaturated case it is possible to modify the weights between saturated and undersaturated equation.
                // a natural since bouth only depend on one saturation.


                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesAnalyticWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

    /// \brief Compute Coats pressure-equation weights for the blackoil model,
    ///        including the extended Rsw/Rvw cross terms.
    ///
    /// Based on Coats (1980), the pressure equation is formed by taking a linear
    /// combination of the component mass-balance equations with coefficients that
    /// convert them into a total reservoir-volume balance:
    ///
    ///   sum_i  w_i * (component i mass balance) = d(total pore volume)/dt + flux terms
    ///
    /// The weights ensure that the weighted combination of storage terms satisfies
    ///
    ///   w_w * invBw * S_w  +  w_o * invBo * S_o  +  w_g * invBg * S_g  =  S_w + S_o + S_g
    ///
    /// For a 3-phase blackoil system the weights are derived by solving:
    ///
    ///   | invBw      0        Rsw*invBw | | ww |   | 1 |
    ///   | 0          invBo    Rs*invBo  | | wo | = | 1 |
    ///   | Rvw*invBg  Rv*invBg invBg     | | wg |   | 1 |
    ///
    /// giving (with D = 1 - Rs*Rv - Rvw*Rsw):
    ///   ww = [ Bw*(1 - Rs*Rv) - Rsw*(Bg - Rv*Bo) ] / D
    ///   wo = [ Bo*(1 - Rvw*Rsw) - Rs*(Bg - Rvw*Bw) ] / D
    ///   wg = ( Bg - Rvw*Bw - Rv*Bo ) / D
    ///
    /// where B_alpha = 1/invB_alpha.  The function handles all active-phase
    /// combinations and the optional dissolved-gas-in-water (Rsw) and
    /// vaporized-water-in-gas (Rvw) extensions.
    ///
    /// Without Rsw and Rvw this formula reduces to getTrueImpesWeightsAnalytic().
    ///
    /// Undersaturation is handled by switching off the appropriate ratio:
    ///   - oil undersaturated (GasMeaning::Rv encoding): set R_s = 0
    ///   - gas undersaturated (GasMeaning::Rs encoding): set R_v = 0
    ///
    /// References:
    ///   - Coats, K.H. (1980). "An Equation of State Compositional Model."
    ///     SPE Journal, 20(5), 363-376. SPE-8284-PA.
    ///     https://doi.org/10.2118/8284-PA
    ///   - Coats, K.H. (1999). "A Note on IMPES and Some IMPES-Based
    ///     Simulation Models." SPE Journal, 5(3), 245-251. SPE-65998-PA.
    ///     https://doi.org/10.2118/65998-PA
    template<class Vector, class ElementContext, class Model, class ElementChunksType>
    void getCoatsWeightsBlackoil(int /*pressureVarIndex*/,
                                 Vector& weights,
                                 const ElementContext& elemCtx,
                                 const Model& model,
                                 const ElementChunksType& element_chunks,
                                 [[maybe_unused]] bool enable_thread_parallel)
    {
        using FluidSystem = typename Model::FluidSystem;
        using LhsEval = double;
        using PrimaryVariables = typename Model::PrimaryVariables;
        using VectorBlockType = typename Vector::block_type;
        using Evaluation =
            typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>::block_type;
        using Toolbox = MathToolbox<Evaluation>;

        const auto& solution = model.solution(/*timeIdx=*/0);
        VectorBlockType bweights;

        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(bweights) if(enable_thread_parallel)
#endif
        for (const auto& chunk : element_chunks) {
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = localElemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();
                const auto& priVars = solution[index];

                bweights = 0.0;

                const bool waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
                const bool oilActive   = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
                const bool gasActive   = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);

                // Formation volume factors  B_alpha = 1 / invB_alpha
                const double Bw = waterActive
                    ? Toolbox::template decay<LhsEval>(1.0 / fs.invB(FluidSystem::waterPhaseIdx))
                    : 0.0;
                const double Bo = oilActive
                    ? Toolbox::template decay<LhsEval>(1.0 / fs.invB(FluidSystem::oilPhaseIdx))
                    : 0.0;
                const double Bg = gasActive
                    ? Toolbox::template decay<LhsEval>(1.0 / fs.invB(FluidSystem::gasPhaseIdx))
                    : 0.0;

                // Solution-gas/vaporized-oil ratios.  Zero out if the primary
                // variable encodes the saturated value instead (switching variables).
                double rs = 0.0, rv = 0.0, rsw = 0.0, rvw = 0.0;

                if (oilActive && gasActive) {
                    if (FluidSystem::enableDissolvedGas()) {
                        rs = Toolbox::template decay<double>(fs.Rs());
                        // When the gas primary variable represents Rv (oil vaporized in gas),
                        // the oil phase is undersaturated, so Rs should be treated as zero.
                        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                            rs = 0.0;
                        }
                    }
                    if (FluidSystem::enableVaporizedOil()) {
                        rv = Toolbox::template decay<double>(fs.Rv());
                        // When the gas primary variable represents Rs (gas dissolved in oil),
                        // the gas phase is undersaturated, so Rv should be treated as zero.
                        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
                            rv = 0.0;
                        }
                    }
                }
                if (waterActive && gasActive) {
                    if (FluidSystem::enableDissolvedGasInWater()) {
                        rsw = Toolbox::template decay<double>(fs.Rsw());
                    }
                    if (FluidSystem::enableVaporizedWater()) {
                        rvw = Toolbox::template decay<double>(fs.Rvw());
                    }
                }

                // Denominator D = 1 - Rs*Rv - Rvw*Rsw
                const double D = 1.0 - rs * rv - rvw * rsw;

                // Compute weights by solving the linear system derived from the
                // requirement that sum_i w_i * M_i = total reservoir saturation.
                //
                // 3-phase (w, o, g) with full cross-terms:
                //   ww = [ Bw*(1 - Rs*Rv) - Rsw*(Bg - Rv*Bo) ] / D
                //   wo = [ Bo*(1 - Rvw*Rsw) - Rs*(Bg - Rvw*Bw) ] / D
                //   wg = ( Bg - Rvw*Bw - Rv*Bo ) / D
                //
                // Degenerate (fewer active phases) cases collapse correctly:
                //   water-only:    ww = Bw
                //   gas+oil:       wo = (Bo - Rs*Bg)/(1-Rs*Rv),  wg = (Bg - Rv*Bo)/(1-Rs*Rv)
                //   water+gas:     ww = (Bw - Rsw*Bg)/(1-Rsw*Rvw), wg = (Bg - Rvw*Bw)/(1-Rsw*Rvw)
                //   water+oil:     ww = Bw, wo = Bo  (no cross-dissolution)

                if (waterActive) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::waterPhaseIdx));
                    bweights[activeCompIdx] = static_cast<LhsEval>(
                        (Bw * (1.0 - rs * rv) - rsw * (Bg - rv * Bo)) / D);
                }
                if (oilActive) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::oilPhaseIdx));
                    bweights[activeCompIdx] = static_cast<LhsEval>(
                        (Bo * (1.0 - rvw * rsw) - rs * (Bg - rvw * Bw)) / D);
                }
                if (gasActive) {
                    const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                        FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                    bweights[activeCompIdx] = static_cast<LhsEval>(
                        (Bg - rvw * Bw - rv * Bo) / D);
                }

                // Normalize so that the maximum absolute weight equals one.
                const double abs_max =
                    *std::ranges::max_element(bweights,
                                              [](double a, double b)
                                              { return std::fabs(a) < std::fabs(b); });
                if (std::fabs(abs_max) > 0.0) {
                    bweights /= std::fabs(abs_max);
                }

                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getCoatsWeightsBlackoil() failed: ",
                                    elemCtx.simulator().vanguard().grid().comm());
    }

    /// \brief Compute Coats pressure-equation weights for the compositional (ptflash) model.
    ///
    /// Based on Coats (1980), the pressure equation is obtained by forming a
    /// volumetric balance from the component mass balances.  Unlike the blackoil
    /// formulation (getCoatsWeightsBlackoil), the compositional model uses
    /// saturation-weighted specific volumes derived from the phase compositions
    /// and densities provided by the ptflash equilibrium calculation.
    ///
    /// References:
    ///   - Coats, K.H. (1980). "An Equation of State Compositional Model."
    ///     SPE Journal, 20(5), 363-376. SPE-8284-PA.
    ///     https://doi.org/10.2118/8284-PA
    ///
    /// For each hydrocarbon
    /// component i with mass-based storage term
    ///
    ///   M_i = phi * sum_{alpha in HC} omega_{alpha,i} * rho_alpha * S_alpha
    ///
    /// the Coats weight is chosen so that
    ///
    ///   sum_i  w_i * M_i = phi * (S_o + S_g)   (total HC pore volume)
    ///
    /// which is satisfied by taking
    ///
    ///   w_i = (sum_{alpha in HC} omega_{alpha,i} * S_alpha)
    ///         / (sum_{alpha in HC} omega_{alpha,i} * rho_alpha * S_alpha)
    ///       = saturation-weighted mass fraction / saturation-weighted mass density
    ///
    /// For the water equation (if present):  w_water = 1 / rho_w
    ///
    /// The energy equation (if present) is left with zero weight because
    /// temperature preconditioning is handled separately.
    template<class Vector, class ElementContext, class Model, class ElementChunksType>
    void getCoatsWeightsCompositional(int /*pressureVarIndex*/,
                                      Vector& weights,
                                      const ElementContext& elemCtx,
                                      const Model& model,
                                      const ElementChunksType& element_chunks,
                                      [[maybe_unused]] bool enable_thread_parallel)
    {
        using FluidSystem = typename Model::FluidSystem;
        using LhsEval = double;
        using VectorBlockType = typename Vector::block_type;
        using Evaluation =
            typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>::block_type;
        using Toolbox = MathToolbox<Evaluation>;
        using Indices = typename Model::Indices;

        constexpr int numComponents = FluidSystem::numComponents;
        constexpr int numPhases     = FluidSystem::numPhases;
        constexpr bool waterEnabled = Indices::waterEnabled;

        VectorBlockType bweights;

        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(bweights) if(enable_thread_parallel)
#endif
        for (const auto& chunk : element_chunks) {
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = localElemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                bweights = 0.0;

                // Iterate over HC components (equations conti0EqIdx + 0 .. numComponents-1).
                // For each component i compute:
                //   num_i = sum_{alpha in HC} massFraction(alpha,i) * saturation(alpha)
                //   den_i = sum_{alpha in HC} massFraction(alpha,i) * density(alpha) * saturation(alpha)
                //   w_i   = num_i / den_i   (= saturation-weighted specific volume)
                //
                // When den_i is zero (component absent from all HC phases), set w_i = 0.

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    double num_i = 0.0;
                    double den_i = 0.0;

                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        // Skip the immiscible water phase for HC component weights
                        if (waterEnabled &&
                            phaseIdx == static_cast<int>(FluidSystem::waterPhaseIdx))
                        {
                            continue;
                        }
                        const double S_alpha =
                            Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx));
                        const double omega_ai =
                            Toolbox::template decay<LhsEval>(fs.massFraction(phaseIdx, compIdx));
                        const double rho_alpha =
                            Toolbox::template decay<LhsEval>(fs.density(phaseIdx));

                        num_i += omega_ai * S_alpha;
                        den_i += omega_ai * rho_alpha * S_alpha;
                    }

                    const int eqIdx = Indices::conti0EqIdx + compIdx;
                    bweights[eqIdx] = (den_i > 0.0) ? static_cast<LhsEval>(num_i / den_i)
                                                    : static_cast<LhsEval>(0.0);
                }

                // Water equation weight: w_water = 1 / rho_w  (ensures phi*S_w contribution)
                if constexpr (waterEnabled) {
                    const int waterEqIdx = Indices::conti0EqIdx + numComponents;
                    const double rho_w =
                        Toolbox::template decay<LhsEval>(
                            fs.density(FluidSystem::waterPhaseIdx));
                    bweights[waterEqIdx] = (rho_w > 0.0)
                        ? static_cast<LhsEval>(1.0 / rho_w)
                        : static_cast<LhsEval>(0.0);
                }

                // Energy equation (if present, numEq > numComponents + waterEnabled):
                // left at zero -- temperature preconditioning is a separate task.

                // Normalize so that the maximum absolute weight equals one.
                const double abs_max =
                    *std::ranges::max_element(bweights,
                                              [](double a, double b)
                                              { return std::fabs(a) < std::fabs(b); });
                if (std::fabs(abs_max) > 0.0) {
                    bweights /= std::fabs(abs_max);
                }

                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getCoatsWeightsCompositional() failed: ",
                                    elemCtx.simulator().vanguard().grid().comm());
    }

} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
