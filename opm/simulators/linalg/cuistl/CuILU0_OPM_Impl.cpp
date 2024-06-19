/*
  Copyright 2024 SINTEF AS

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
#include <cuda.h>
#include <cuda_runtime.h>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/coloringAndReorderingUtils.hpp>
#include <opm/simulators/linalg/cuistl/CuILU0_OPM_Impl.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/preconditionerKernels/ILU0Kernels.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <vector>

#include <config.h>
#include <chrono>
#include <tuple>
namespace Opm::cuistl
{

template <class M, class X, class Y, int l>
CuILU0_OPM_Impl<M, X, Y, l>::CuILU0_OPM_Impl(const M& A, bool splitMatrix, bool tuneKernels)
    : m_cpuMatrix(A)
    , m_levelSets(Opm::getMatrixRowColoring(m_cpuMatrix, Opm::ColoringType::LOWER))
    , m_reorderedToNatural(detail::createReorderedToNatural(m_levelSets))
    , m_naturalToReordered(detail::createNaturalToReordered(m_levelSets))
    , m_gpuMatrix(CuSparseMatrix<field_type>::fromMatrix(m_cpuMatrix, true))
    , m_gpuMatrixReordered(nullptr)
    , m_gpuMatrixReorderedLower(nullptr)
    , m_gpuMatrixReorderedUpper(nullptr)
    , m_gpuNaturalToReorder(m_naturalToReordered)
    , m_gpuReorderToNatural(m_reorderedToNatural)
    , m_gpuDInv(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize())
    , m_splitMatrix(splitMatrix)
    , m_tuneThreadBlockSizes(tuneKernels)
{
    // TODO: Should in some way verify that this matrix is symmetric, only do it debug mode?
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_gpuMatrix.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_gpuMatrix.N(), A.N()));
    OPM_ERROR_IF(A[0][0].N() != m_gpuMatrix.blockSize(),
                 fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.blockSize(),
                             A[0][0].N()));
    OPM_ERROR_IF(A.N() * A[0][0].N() != m_gpuMatrix.dim(),
                 fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.dim(),
                             A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != m_gpuMatrix.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_gpuMatrix.nonzeroes(),
                             A.nonzeroes()));
    if (m_splitMatrix) {
        m_gpuMatrixReorderedDiag.emplace(CuVector<field_type>(blocksize_ * blocksize_ * m_cpuMatrix.N()));
        detail::extractLowerAndUpperMatrices<M, field_type, CuSparseMatrix<field_type>>(
            m_cpuMatrix, m_reorderedToNatural, m_gpuMatrixReorderedLower, m_gpuMatrixReorderedUpper);
    } else {
        detail::createReorderedMatrix<M, field_type, CuSparseMatrix<field_type>>(
            m_cpuMatrix, m_reorderedToNatural, m_gpuReorderedLU);
    }
    computeDiagAndMoveReorderedData();

#ifdef USE_HIP
    if (m_tuneThreadBlockSizes){
        tuneThreadBlockSizes();
    }
#endif
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::apply(X& v, const Y& d)
{
    // ScopeTimer timer("Apply");
    OPM_TIMEBLOCK(prec_apply);
    {
        int levelStartIdx = 0;
        for (int level = 0; level < m_levelSets.size(); ++level) {
            const int numOfRowsInLevel = m_levelSets[level].size();
            if (m_splitMatrix) {
                detail::ILU0::ILULowerSolveLevelSetSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLower->getRowIndices().data(),
                    m_gpuMatrixReorderedLower->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuMatrixReorderedDiag.value().data(),
                    d.data(),
                    v.data(),
                    m_applyThreadBlockSize);
            } else {
                detail::ILU0::ILULowerSolveLevelSet<field_type, blocksize_>(
                    m_gpuReorderedLU->getNonZeroValues().data(),
                    m_gpuReorderedLU->getRowIndices().data(),
                    m_gpuReorderedLU->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    d.data(),
                    v.data(),
                    m_applyThreadBlockSize);
            }
            levelStartIdx += numOfRowsInLevel;
        }

        levelStartIdx = m_cpuMatrix.N();
        //  upper triangular solve: (D + U_A) v = Dy
        for (int level = m_levelSets.size() - 1; level >= 0; --level) {
            const int numOfRowsInLevel = m_levelSets[level].size();
            levelStartIdx -= numOfRowsInLevel;
            if (m_splitMatrix) {
                detail::ILU0::ILUUpperSolveLevelSetSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpper->getRowIndices().data(),
                    m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuMatrixReorderedDiag.value().data(),
                    v.data(),
                    m_applyThreadBlockSize);
            } else {
                detail::ILU0::ILUUpperSolveLevelSet<field_type, blocksize_>(
                    m_gpuReorderedLU->getNonZeroValues().data(),
                    m_gpuReorderedLU->getRowIndices().data(),
                    m_gpuReorderedLU->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    v.data(),
                    m_applyThreadBlockSize);
            }
        }
    }
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
CuILU0_OPM_Impl<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::update()
{
    OPM_TIMEBLOCK(prec_update);
    {
        m_gpuMatrix.updateNonzeroValues(m_cpuMatrix, true); // send updated matrix to the gpu
        computeDiagAndMoveReorderedData();
    }
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::computeDiagAndMoveReorderedData()
{
    OPM_TIMEBLOCK(prec_update);
    {
        if (m_splitMatrix) {
            detail::copyMatDataToReorderedSplit<field_type, blocksize_>(
                m_gpuMatrix.getNonZeroValues().data(),
                m_gpuMatrix.getRowIndices().data(),
                m_gpuMatrix.getColumnIndices().data(),
                m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                m_gpuMatrixReorderedLower->getRowIndices().data(),
                m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                m_gpuMatrixReorderedUpper->getRowIndices().data(),
                m_gpuMatrixReorderedDiag.value().data(),
                m_gpuNaturalToReorder.data(),
                m_gpuMatrixReorderedLower->N(),
                m_updateThreadBlockSize);
        } else {
            detail::copyMatDataToReordered<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                                    m_gpuMatrix.getRowIndices().data(),
                                                                    m_gpuReorderedLU->getNonZeroValues().data(),
                                                                    m_gpuReorderedLU->getRowIndices().data(),
                                                                    m_gpuNaturalToReorder.data(),
                                                                    m_gpuReorderedLU->N(),
                                                                    m_updateThreadBlockSize);
        }
        int levelStartIdx = 0;
        for (int level = 0; level < m_levelSets.size(); ++level) {
            const int numOfRowsInLevel = m_levelSets[level].size();

            if (m_splitMatrix) {
                detail::ILU0::LUFactorizationSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLower->getRowIndices().data(),
                    m_gpuMatrixReorderedLower->getColumnIndices().data(),
                    m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpper->getRowIndices().data(),
                    m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                    m_gpuMatrixReorderedDiag.value().data(),
                    m_gpuReorderToNatural.data(),
                    m_gpuNaturalToReorder.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_updateThreadBlockSize);

            } else {
                detail::ILU0::LUFactorization<field_type, blocksize_>(m_gpuReorderedLU->getNonZeroValues().data(),
                                                                        m_gpuReorderedLU->getRowIndices().data(),
                                                                        m_gpuReorderedLU->getColumnIndices().data(),
                                                                        m_gpuNaturalToReorder.data(),
                                                                        m_gpuReorderToNatural.data(),
                                                                        numOfRowsInLevel,
                                                                        levelStartIdx,
                                                                        m_updateThreadBlockSize);
            }
            levelStartIdx += numOfRowsInLevel;
        }
    }
}

template <class M, class X, class Y, int l>
void
CuILU0_OPM_Impl<M, X, Y, l>::tuneThreadBlockSizes()
{
    auto start = std::chrono::high_resolution_clock::now();
    long long bestApplyTime = __LONG_LONG_MAX__;
    long long bestUpdateTime = __LONG_LONG_MAX__;
    int bestApplyBlockSize = -1;
    int bestUpdateBlockSize = -1;
    int interval = 64;

    //temporary buffers for the apply
    CuVector<field_type> tmpV(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    CuVector<field_type> tmpD(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    tmpD = 1;

    for (int thrBlockSize = interval; thrBlockSize <= 1024; thrBlockSize += interval){
        // sometimes the first kernel launch kan be slower, so take the time twice
        for (int i = 0; i < 2; ++i){

            auto beforeUpdate = std::chrono::high_resolution_clock::now();
            m_updateThreadBlockSize = thrBlockSize;
            update();
            std::ignore = cudaDeviceSynchronize();
            auto afterUpdate = std::chrono::high_resolution_clock::now();
            if (cudaSuccess == cudaGetLastError()){ // kernel launch was valid
                long long durationInMicroSec = std::chrono::duration_cast<std::chrono::microseconds>(afterUpdate - beforeUpdate).count();
                if (durationInMicroSec < bestUpdateTime){
                    bestUpdateTime = durationInMicroSec;
                    bestUpdateBlockSize = thrBlockSize;
                }
            }

            auto beforeApply = std::chrono::high_resolution_clock::now();
            m_applyThreadBlockSize = thrBlockSize;
            apply(tmpV, tmpD);
            std::ignore = cudaDeviceSynchronize();
            auto afterApply = std::chrono::high_resolution_clock::now();
            if (cudaSuccess == cudaGetLastError()){ // kernel launch was valid
                long long durationInMicroSec = std::chrono::duration_cast<std::chrono::microseconds>(afterApply - beforeApply).count();
                if (durationInMicroSec < bestApplyTime){
                    bestApplyTime = durationInMicroSec;
                    bestApplyBlockSize = thrBlockSize;
                }
            }
        }
    }

    m_applyThreadBlockSize = bestApplyBlockSize;
    m_updateThreadBlockSize = bestUpdateBlockSize;
    auto end = std::chrono::high_resolution_clock::now();
    long long durationInMicroSec = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    printf("calibration time %lld, results: {apply: %d, update: %d}\n", durationInMicroSec, m_applyThreadBlockSize, m_updateThreadBlockSize);
}

} // namespace Opm::cuistl
#define INSTANTIATE_CUDILU_DUNE(realtype, blockdim)                                                                    \
    template class ::Opm::cuistl::CuILU0_OPM_Impl<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,            \
                                         ::Opm::cuistl::CuVector<realtype>,                                            \
                                         ::Opm::cuistl::CuVector<realtype>>;                                           \
    template class ::Opm::cuistl::CuILU0_OPM_Impl<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,             \
                                         ::Opm::cuistl::CuVector<realtype>,                                            \
                                         ::Opm::cuistl::CuVector<realtype>>

INSTANTIATE_CUDILU_DUNE(double, 1);
INSTANTIATE_CUDILU_DUNE(double, 2);
INSTANTIATE_CUDILU_DUNE(double, 3);
INSTANTIATE_CUDILU_DUNE(double, 4);
INSTANTIATE_CUDILU_DUNE(double, 5);
INSTANTIATE_CUDILU_DUNE(double, 6);

INSTANTIATE_CUDILU_DUNE(float, 1);
INSTANTIATE_CUDILU_DUNE(float, 2);
INSTANTIATE_CUDILU_DUNE(float, 3);
INSTANTIATE_CUDILU_DUNE(float, 4);
INSTANTIATE_CUDILU_DUNE(float, 5);
INSTANTIATE_CUDILU_DUNE(float, 6);
