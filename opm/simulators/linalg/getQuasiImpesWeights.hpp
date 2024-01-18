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

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <algorithm>
#include <cmath>

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
    void getQuasiImpesWeights(const Matrix& matrix, const int pressureVarIndex, const bool transpose, Vector& weights)
    {
        using VectorBlockType = typename Vector::block_type;
        using MatrixBlockType = typename Matrix::block_type;
        const Matrix& A = matrix;
        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;
        const auto endi = A.end();
        for (auto i = A.begin(); i != endi; ++i) {
            const auto endj = (*i).end();
            MatrixBlockType diag_block(0.0);
            for (auto j = (*i).begin(); j != endj; ++j) {
                if (i.index() == j.index()) {
                    diag_block = (*j);
                    break;
                }
            }
            VectorBlockType bweights;
            if (transpose) {
                diag_block.solve(bweights, rhs);
            } else {
                auto diag_block_transpose = Details::transposeDenseMatrix(diag_block);
                diag_block_transpose.solve(bweights, rhs);
            }
            double abs_max = *std::max_element(
                bweights.begin(), bweights.end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); });
            bweights /= std::fabs(abs_max);
            weights[i.index()] = bweights;
        }
        // return weights;
    }

    template <class Matrix, class Vector>
    Vector getQuasiImpesWeights(const Matrix& matrix, const int pressureVarIndex, const bool transpose)
    {
        Vector weights(matrix.N());
        getQuasiImpesWeights(matrix, pressureVarIndex, transpose, weights);
        return weights;
    }

    template<class Vector, class GridView, class ElementContext, class Model>
    void getTrueImpesWeights(int pressureVarIndex, Vector& weights, const GridView& gridView,
                             ElementContext& elemCtx, const Model& model, std::size_t threadId)
    {
        using VectorBlockType = typename Vector::block_type;
        using Matrix = typename std::decay_t<decltype(model.linearizer().jacobian())>;
        using MatrixBlockType = typename Matrix::MatrixBlock;
        constexpr int numEq = VectorBlockType::size();
        using Evaluation = typename std::decay_t<decltype(model.localLinearizer(threadId).localResidual().residual(0))>
            ::block_type;
        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;
        int index = 0;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            Dune::FieldVector<Evaluation, numEq> storage;
            model.localLinearizer(threadId).localResidual().computeStorage(storage,elemCtx,/*spaceIdx=*/0, /*timeIdx=*/0);
            auto extrusionFactor = elemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
            auto scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;
            auto storage_scale = scvVolume / elemCtx.simulator().timeStepSize();
            MatrixBlockType block;
            double pressure_scale = 50e5;
            for (int ii = 0; ii < numEq; ++ii) {
                for (int jj = 0; jj < numEq; ++jj) {
                    block[ii][jj] = storage[ii].derivative(jj)/storage_scale;
                    if (jj == pressureVarIndex) {
                        block[ii][jj] *= pressure_scale;
                    }
                }
            }
            VectorBlockType bweights;
            MatrixBlockType block_transpose = Details::transposeDenseMatrix(block);
            block_transpose.solve(bweights, rhs);
            double abs_max = *std::max_element(
                bweights.begin(), bweights.end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); });
            // probably a scaling which could give approximately total compressibility would be better
            bweights /=  std::fabs(abs_max); // given normal densities this scales weights to about 1.

            weights[index] = bweights;
            ++index;
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

    template <class Vector, class GridView, class ElementContext, class Model>
    void getTrueImpesWeightsAnalytic(int pressureVarIndex,
                                     Vector& weights,
                                     const GridView& gridView,
                                     ElementContext& elemCtx,
                                     const Model& model,
                                     std::size_t threadId)
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
        using Indices = typename Model::Indices;

        using PrimaryVariables = typename Model::PrimaryVariables;
        using VectorBlockType = typename Vector::block_type;
        constexpr int numEq = VectorBlockType::size();
        using Evaluation =
            typename std::decay_t<decltype(model.localLinearizer(threadId).localResidual().residual(0))>::block_type;
        using Toolbox = MathToolbox<Evaluation>;
        VectorBlockType rhs(0.0);
        const auto& solution = model.solution(/*timeIdx*/ 0);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& index = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            VectorBlockType bweights;

            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(
                    FluidSystem::solventComponentIndex(FluidSystem::waterPhaseIdx));
                bweights[FluidSystem::waterPhaseIdx]
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
                unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(
                    FluidSystem::solventComponentIndex(FluidSystem::oilPhaseIdx));
                bweights[activeCompIdx] = Toolbox::template decay<LhsEval>(
                    (1 / fs.invB(FluidSystem::oilPhaseIdx) - rs / fs.invB(FluidSystem::gasPhaseIdx))
                    / denominator);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(
                    FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                bweights[activeCompIdx] = Toolbox::template decay<LhsEval>(
                    (1 / fs.invB(FluidSystem::gasPhaseIdx) - rv / fs.invB(FluidSystem::oilPhaseIdx))
                    / denominator);
            }
            weights[index] = bweights;
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesAnalyticWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }
} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
