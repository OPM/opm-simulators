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

    template<class Evaluation, class Vector, class GridView, class ElementContext, class Model>
    void getTrueImpesWeights(int pressureVarIndex, Vector& weights, const GridView& gridView,
                             ElementContext& elemCtx, const Model& model, std::size_t threadId)
    {
        using VectorBlockType = typename Vector::block_type;
        using Matrix = typename std::decay_t<decltype(model.linearizer().jacobian())>;
        using MatrixBlockType = typename Matrix::MatrixBlock;
        constexpr int numEq = VectorBlockType::size();
//        using Evaluation = typename std::decay_t<decltype(model.localLinearizer(threadId).localResidual().residual(0))>
//            ::block_type;
        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;
        int index = 0;
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updatePrimaryStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            Dune::FieldVector<Evaluation, numEq> storage;
            model.localLinearizer(threadId).localResidual().computeStorage(storage,elemCtx,/*spaceIdx=*/0, /*timeIdx=*/0);
            //auto extrusionFactor = elemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
            auto scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume();// * extrusionFactor;
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
            bweights /= 1000.0; // given normal densities this scales weights to about 1.
            weights[index] = bweights;
            ++index;
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

    template<class Evaluation, class Vector, class Model>
    void getTrueImpesWeights(int pressureVarIndex, Vector& weights, const Model& model)
    {
        using VectorBlockType = typename Vector::block_type;
        using Matrix = typename std::decay_t<decltype(model.linearizer().jacobian())>;
        using MatrixBlockType = typename Matrix::MatrixBlock;
        constexpr int numEq = VectorBlockType::size();
        unsigned numCells = model.numTotalDof();
        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;
        //NB !!OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP        
#pragma omp parallel for
#endif        
        for(unsigned globI = 0; globI < numCells; globI++){
            Dune::FieldVector<Evaluation, numEq> storage;
            const auto* intQuantsInP = model.cachedIntensiveQuantities(globI, /*timeIdx*/0);
            assert(intQuantsInP);
            const auto& intQuantsIn = *intQuantsInP;
            Model::LocalResidual::computeStorage(storage,intQuantsIn, 0);           
            double scvVolume = model.dofTotalVolume(globI);
            double dt = 3600*24;
            auto storage_scale = scvVolume / dt;
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
            bweights /= 1000.0; // given normal densities this scales weights to about 1.
            weights[globI] = bweights;
        }
        //NB!! OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
