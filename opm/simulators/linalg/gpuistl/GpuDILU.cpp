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
#include <chrono>
#include <config.h>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <fmt/core.h>
#include <limits>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/GraphColoring.hpp>
#include <opm/simulators/linalg/gpuistl/detail/autotuner.hpp>
#include <opm/simulators/linalg/gpuistl/GpuDILU.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/coloringAndReorderingUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/DILUKernels.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <tuple>
#include <functional>
#include <utility>
#include <string>
namespace Opm::gpuistl
{

template <class M, class X, class Y, int l>
GpuDILU<M, X, Y, l>::GpuDILU(const M& A, bool splitMatrix, bool tuneKernels, bool storeFactorizationAsFloat)
    : m_cpuMatrix(A)
    , m_levelSets(Opm::getMatrixRowColoring(m_cpuMatrix, Opm::ColoringType::LOWER))
    , m_reorderedToNatural(detail::createReorderedToNatural(m_levelSets))
    , m_naturalToReordered(detail::createNaturalToReordered(m_levelSets))
    , m_gpuMatrix(GpuSparseMatrix<field_type>::fromMatrix(m_cpuMatrix, true))
    , m_gpuNaturalToReorder(m_naturalToReordered)
    , m_gpuReorderToNatural(m_reorderedToNatural)
    , m_gpuDInv(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize())
    , m_splitMatrix(splitMatrix)
    , m_tuneThreadBlockSizes(tuneKernels)
    , m_storeFactorizationAsFloat(storeFactorizationAsFloat)

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
        m_gpuMatrixReorderedDiag = std::make_unique<GpuVector<field_type>>(blocksize_ * blocksize_ * m_cpuMatrix.N());
        std::tie(m_gpuMatrixReorderedLower, m_gpuMatrixReorderedUpper)
            = detail::extractLowerAndUpperMatrices<M, field_type, GpuSparseMatrix<field_type>>(m_cpuMatrix,
                                                                                              m_reorderedToNatural);
    }
    else {
        m_gpuMatrixReordered = detail::createReorderedMatrix<M, field_type, GpuSparseMatrix<field_type>>(
            m_cpuMatrix, m_reorderedToNatural);
    }

    if (m_storeFactorizationAsFloat) {
        OPM_THROW(std::runtime_error, "Matrix must be split when storing as float.");
        m_gpuMatrixReorderedLowerFloat = std::make_unique<FloatMat>(m_gpuMatrixReorderedLower->getRowIndices(), m_gpuMatrixReorderedLower->getColumnIndices(), blocksize_);
        m_gpuMatrixReorderedUpperFloat = std::make_unique<FloatMat>(m_gpuMatrixReorderedUpper->getRowIndices(), m_gpuMatrixReorderedUpper->getColumnIndices(), blocksize_);
        m_gpuMatrixReorderedDiagFloat = std::make_unique<FloatVec>(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize());
    }

    computeDiagAndMoveReorderedData(m_moveThreadBlockSize, m_DILUFactorizationThreadBlockSize);

    if (m_tuneThreadBlockSizes) {
        tuneThreadBlockSizes();
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::apply(X& v, const Y& d)
{
    OPM_TIMEBLOCK(prec_apply);
    {
        apply(v, d, m_lowerSolveThreadBlockSize, m_upperSolveThreadBlockSize);
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::apply(X& v, const Y& d, int lowerSolveThreadBlockSize, int upperSolveThreadBlockSize)
{
    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        if (m_splitMatrix) {
            detail::DILU::solveLowerLevelSetSplit<field_type, blocksize_>(
                m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                m_gpuMatrixReorderedLower->getRowIndices().data(),
                m_gpuMatrixReorderedLower->getColumnIndices().data(),
                m_gpuReorderToNatural.data(),
                levelStartIdx,
                numOfRowsInLevel,
                m_gpuDInv.data(),
                d.data(),
                v.data(),
                lowerSolveThreadBlockSize);
        } else {
            detail::DILU::solveLowerLevelSet<field_type, blocksize_>(
                m_gpuMatrixReordered->getNonZeroValues().data(),
                m_gpuMatrixReordered->getRowIndices().data(),
                m_gpuMatrixReordered->getColumnIndices().data(),
                m_gpuReorderToNatural.data(),
                levelStartIdx,
                numOfRowsInLevel,
                m_gpuDInv.data(),
                d.data(),
                v.data(),
                lowerSolveThreadBlockSize);
        }
        levelStartIdx += numOfRowsInLevel;
    }

    levelStartIdx = m_cpuMatrix.N();
    //  upper triangular solve: (D + U_A) v = Dy
    for (int level = m_levelSets.size() - 1; level >= 0; --level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        levelStartIdx -= numOfRowsInLevel;
        if (m_splitMatrix) {
            detail::DILU::solveUpperLevelSetSplit<field_type, blocksize_>(
                m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                m_gpuMatrixReorderedUpper->getRowIndices().data(),
                m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                m_gpuReorderToNatural.data(),
                levelStartIdx,
                numOfRowsInLevel,
                m_gpuDInv.data(),
                v.data(),
                upperSolveThreadBlockSize);
        } else {
            detail::DILU::solveUpperLevelSet<field_type, blocksize_>(
                m_gpuMatrixReordered->getNonZeroValues().data(),
                m_gpuMatrixReordered->getRowIndices().data(),
                m_gpuMatrixReordered->getColumnIndices().data(),
                m_gpuReorderToNatural.data(),
                levelStartIdx,
                numOfRowsInLevel,
                m_gpuDInv.data(),
                v.data(),
                upperSolveThreadBlockSize);
        }
    }
}


template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
GpuDILU<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::update()
{
    OPM_TIMEBLOCK(prec_update);
    {
        update(m_moveThreadBlockSize, m_DILUFactorizationThreadBlockSize);
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::update(int moveThreadBlockSize, int factorizationBlockSize)
{
    m_gpuMatrix.updateNonzeroValues(m_cpuMatrix, true); // send updated matrix to the gpu
    computeDiagAndMoveReorderedData(moveThreadBlockSize, factorizationBlockSize);
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::computeDiagAndMoveReorderedData(int moveThreadBlockSize, int factorizationBlockSize)
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
            moveThreadBlockSize);
    } else {
        detail::copyMatDataToReordered<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                                m_gpuMatrix.getRowIndices().data(),
                                                                m_gpuMatrixReordered->getNonZeroValues().data(),
                                                                m_gpuMatrixReordered->getRowIndices().data(),
                                                                m_gpuNaturalToReorder.data(),
                                                                m_gpuMatrixReordered->N(),
                                                                moveThreadBlockSize);
    }

    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        if (m_splitMatrix) {
            detail::DILU::computeDiluDiagonalSplit<field_type, blocksize_>(
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
                factorizationBlockSize);
        } else {
            detail::DILU::computeDiluDiagonal<field_type, blocksize_>(
                m_gpuMatrixReordered->getNonZeroValues().data(),
                m_gpuMatrixReordered->getRowIndices().data(),
                m_gpuMatrixReordered->getColumnIndices().data(),
                m_gpuReorderToNatural.data(),
                m_gpuNaturalToReorder.data(),
                levelStartIdx,
                numOfRowsInLevel,
                m_gpuDInv.data(),
                factorizationBlockSize);
        }
        levelStartIdx += numOfRowsInLevel;
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::tuneThreadBlockSizes()
{
    // tune the thread-block size of the update function
    auto tuneMoveThreadBlockSizeInUpdate = [this](int moveThreadBlockSize){
        this->update(moveThreadBlockSize, m_DILUFactorizationThreadBlockSize);
    };
    m_moveThreadBlockSize = detail::tuneThreadBlockSize(tuneMoveThreadBlockSizeInUpdate, "Kernel moving data to reordered matrix");

    auto tuneFactorizationThreadBlockSizeInUpdate = [this](int factorizationThreadBlockSize){
        this->update(m_moveThreadBlockSize, factorizationThreadBlockSize);
    };
    m_DILUFactorizationThreadBlockSize = detail::tuneThreadBlockSize(tuneFactorizationThreadBlockSizeInUpdate, "Kernel computing DILU factorization");

    // tune the thread-block size of the apply
    GpuVector<field_type> tmpV(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    GpuVector<field_type> tmpD(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    tmpD = 1;

    auto tuneLowerSolveThreadBlockSizeInApply = [this, &tmpV, &tmpD](int lowerSolveThreadBlockSize){
        this->apply(tmpV, tmpD, lowerSolveThreadBlockSize, m_DILUFactorizationThreadBlockSize);
    };
    m_lowerSolveThreadBlockSize = detail::tuneThreadBlockSize(tuneLowerSolveThreadBlockSizeInApply, "Kernel computing a lower triangular solve for a level set");

    auto tuneUpperSolveThreadBlockSizeInApply = [this, &tmpV, &tmpD](int upperSolveThreadBlockSize){
        this->apply(tmpV, tmpD, m_moveThreadBlockSize, upperSolveThreadBlockSize);
    };
    m_upperSolveThreadBlockSize = detail::tuneThreadBlockSize(tuneUpperSolveThreadBlockSizeInApply, "Kernel computing an upper triangular solve for a level set");
}

} // namespace Opm::gpuistl
#define INSTANTIATE_CUDILU_DUNE(realtype, blockdim)                                                                    \
    template class ::Opm::gpuistl::GpuDILU<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,            \
                                         ::Opm::gpuistl::GpuVector<realtype>,                                            \
                                         ::Opm::gpuistl::GpuVector<realtype>>;                                           \
    template class ::Opm::gpuistl::GpuDILU<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,             \
                                         ::Opm::gpuistl::GpuVector<realtype>,                                            \
                                         ::Opm::gpuistl::GpuVector<realtype>>

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
