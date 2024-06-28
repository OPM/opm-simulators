/*
  Copyright 2022-2023 SINTEF AS

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
#include <opm/simulators/linalg/cuistl/CuDILU.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <vector>
#include <config.h>
#include <chrono>
#include <tuple>


namespace
{
std::vector<int>
createReorderedToNatural(Opm::SparseTable<size_t>& levelSets)
{
    auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[globCnt++] = static_cast<int>(col);
        }
    }
    return res;
}

std::vector<int>
createNaturalToReordered(Opm::SparseTable<size_t>& levelSets)
{
    auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[col] = globCnt++;
        }
    }
    return res;
}

template <class M, class field_type, class GPUM>
void
createReorderedMatrix(const M& naturalMatrix,
                      std::vector<int>& reorderedToNatural,
                      std::unique_ptr<GPUM>& reorderedGpuMat)
{
    M reorderedMatrix(naturalMatrix.N(), naturalMatrix.N(), naturalMatrix.nonzeroes(), M::row_wise);
    for (auto dstRowIt = reorderedMatrix.createbegin(); dstRowIt != reorderedMatrix.createend(); ++dstRowIt) {
        auto srcRow = naturalMatrix.begin() + reorderedToNatural[dstRowIt.index()];
        // For elements in A
        for (auto elem = srcRow->begin(); elem != srcRow->end(); elem++) {
            dstRowIt.insert(elem.index());
        }
    }

    reorderedGpuMat.reset(new auto (GPUM::fromMatrix(reorderedMatrix, true)));
}

template <class M, class field_type, class GPUM>
void
extractLowerAndUpperMatrices(const M& naturalMatrix,
                             std::vector<int>& reorderedToNatural,
                             std::unique_ptr<GPUM>& lower,
                             std::unique_ptr<GPUM>& upper)
{
    const size_t new_nnz = (naturalMatrix.nonzeroes() - naturalMatrix.N()) / 2;

    M reorderedLower(naturalMatrix.N(), naturalMatrix.N(), new_nnz, M::row_wise);
    M reorderedUpper(naturalMatrix.N(), naturalMatrix.N(), new_nnz, M::row_wise);

    for (auto lowerIt = reorderedLower.createbegin(), upperIt = reorderedUpper.createbegin();
         lowerIt != reorderedLower.createend();
         ++lowerIt, ++upperIt) {

        auto srcRow = naturalMatrix.begin() + reorderedToNatural[lowerIt.index()];

        for (auto elem = srcRow->begin(); elem != srcRow->end(); ++elem) {
            if (elem.index() < srcRow.index()) { // add index to lower matrix if under the diagonal
                lowerIt.insert(elem.index());
            } else if (elem.index() > srcRow.index()) { // add element to upper matrix if above the diagonal
                upperIt.insert(elem.index());
            }
        }
    }

    lower.reset(new auto (GPUM::fromMatrix(reorderedLower, true)));
    upper.reset(new auto (GPUM::fromMatrix(reorderedUpper, true)));
    return;
}

} // NAMESPACE

namespace Opm::cuistl
{

template <class M, class X, class Y, int l>
CuDILU<M, X, Y, l>::CuDILU(const M& A, bool splitMatrix, bool tuneKernels)
    : m_cpuMatrix(A)
    , m_levelSets(Opm::getMatrixRowColoring(m_cpuMatrix, Opm::ColoringType::LOWER))
    , m_reorderedToNatural(createReorderedToNatural(m_levelSets))
    , m_naturalToReordered(createNaturalToReordered(m_levelSets))
    , m_gpuMatrix(CuSparseMatrix<field_type>::fromMatrix(m_cpuMatrix, true))
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
        m_gpuMatrixReorderedDiag.reset(new auto(CuVector<field_type>(blocksize_ * blocksize_ * m_cpuMatrix.N())));
        extractLowerAndUpperMatrices<M, field_type, CuSparseMatrix<field_type>>(
            m_cpuMatrix, m_reorderedToNatural, m_gpuMatrixReorderedLower, m_gpuMatrixReorderedUpper);
    } else {
        createReorderedMatrix<M, field_type, CuSparseMatrix<field_type>>(
            m_cpuMatrix, m_reorderedToNatural, m_gpuMatrixReordered);
    }
    computeDiagAndMoveReorderedData();

    // HIP does currently not support automtically picking thread block sizes as well as CUDA
    // So only when tuning and using hip should we do our own manual tuning
#ifdef USE_HIP
    if (m_tuneThreadBlockSizes){
        tuneThreadBlockSizes();
    }
#endif
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::apply(X& v, const Y& d)
{
    OPM_TIMEBLOCK(prec_apply);
    {
        int levelStartIdx = 0;
        for (int level = 0; level < m_levelSets.size(); ++level) {
            const int numOfRowsInLevel = m_levelSets[level].size();
            if (m_splitMatrix) {
                detail::computeLowerSolveLevelSetSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLower->getRowIndices().data(),
                    m_gpuMatrixReorderedLower->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    d.data(),
                    v.data(),
                    m_applyThreadBlockSize);
            } else {
                detail::computeLowerSolveLevelSet<field_type, blocksize_>(
                    m_gpuMatrixReordered->getNonZeroValues().data(),
                    m_gpuMatrixReordered->getRowIndices().data(),
                    m_gpuMatrixReordered->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
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
                detail::computeUpperSolveLevelSetSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpper->getRowIndices().data(),
                    m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    v.data(),
                    m_applyThreadBlockSize);
            } else {
                detail::computeUpperSolveLevelSet<field_type, blocksize_>(
                    m_gpuMatrixReordered->getNonZeroValues().data(),
                    m_gpuMatrixReordered->getRowIndices().data(),
                    m_gpuMatrixReordered->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    v.data(),
                    m_applyThreadBlockSize);
            }
        }
    }
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
CuDILU<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::update()
{
    OPM_TIMEBLOCK(prec_update);
    {
        m_gpuMatrix.updateNonzeroValues(m_cpuMatrix, true); // send updated matrix to the gpu
        computeDiagAndMoveReorderedData();
    }
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::computeDiagAndMoveReorderedData()
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
                m_gpuMatrixReorderedDiag->data(),
                m_gpuNaturalToReorder.data(),
                m_gpuMatrixReorderedLower->N(),
                m_updateThreadBlockSize);
        } else {
            detail::copyMatDataToReordered<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                                   m_gpuMatrix.getRowIndices().data(),
                                                                   m_gpuMatrixReordered->getNonZeroValues().data(),
                                                                   m_gpuMatrixReordered->getRowIndices().data(),
                                                                   m_gpuNaturalToReorder.data(),
                                                                   m_gpuMatrixReordered->N(),
                                                                   m_updateThreadBlockSize);
        }

        int levelStartIdx = 0;
        for (int level = 0; level < m_levelSets.size(); ++level) {
            const int numOfRowsInLevel = m_levelSets[level].size();
            if (m_splitMatrix) {
                detail::computeDiluDiagonalSplit<field_type, blocksize_>(
                    m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLower->getRowIndices().data(),
                    m_gpuMatrixReorderedLower->getColumnIndices().data(),
                    m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpper->getRowIndices().data(),
                    m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                    m_gpuMatrixReorderedDiag->data(),
                    m_gpuReorderToNatural.data(),
                    m_gpuNaturalToReorder.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    m_updateThreadBlockSize);
            } else {
                detail::computeDiluDiagonal<field_type, blocksize_>(m_gpuMatrixReordered->getNonZeroValues().data(),
                                                                    m_gpuMatrixReordered->getRowIndices().data(),
                                                                    m_gpuMatrixReordered->getColumnIndices().data(),
                                                                    m_gpuReorderToNatural.data(),
                                                                    m_gpuNaturalToReorder.data(),
                                                                    levelStartIdx,
                                                                    numOfRowsInLevel,
                                                                    m_gpuDInv.data(),
                                                                    m_updateThreadBlockSize);
            }
            levelStartIdx += numOfRowsInLevel;
        }
    }
}

template <class M, class X, class Y, int l>
void
CuDILU<M, X, Y, l>::tuneThreadBlockSizes()
{
    // TODO: generalize this code and put it somewhere outside of this class
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
}

} // namespace Opm::cuistl
#define INSTANTIATE_CUDILU_DUNE(realtype, blockdim)                                                                    \
    template class ::Opm::cuistl::CuDILU<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,            \
                                         ::Opm::cuistl::CuVector<realtype>,                                            \
                                         ::Opm::cuistl::CuVector<realtype>>;                                           \
    template class ::Opm::cuistl::CuDILU<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,             \
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
