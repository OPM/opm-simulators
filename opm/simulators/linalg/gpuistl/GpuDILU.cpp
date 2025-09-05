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
#include <functional>
#include <limits>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/GraphColoring.hpp>
#include <opm/simulators/linalg/gpuistl/GpuDILU.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/autotuner.hpp>
#include <opm/simulators/linalg/gpuistl/detail/coloringAndReorderingUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/DILUKernels.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <string>
#include <tuple>
#include <utility>

namespace Opm::gpuistl
{

template <class M, class X, class Y, int l>
GpuDILU<M, X, Y, l>::GpuDILU(const typename GpuDILU<M, X, Y, l>::GPUMatrix& gpuMatrix,
                             const M& cpuMatrix,
                             bool splitMatrix,
                             bool tuneKernels,
                             int mixedPrecisionScheme,
                             bool reorder)
    : m_levelSets(Opm::getMatrixRowColoring(cpuMatrix, Opm::ColoringType::LOWER))
    , m_reorderedToNatural(detail::createReorderedToNatural(m_levelSets))
    , m_naturalToReordered(detail::createNaturalToReordered(m_levelSets))
    , m_gpuMatrix(gpuMatrix)
    , m_gpuNaturalToReorder(m_naturalToReordered)
    , m_gpuReorderToNatural(m_reorderedToNatural)
    , m_gpuLevelSets(m_levelSets.dataPtr(), m_levelSets.dataSize())
    , m_gpuDInv(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize())
    , m_splitMatrix(splitMatrix)
    , m_tuneThreadBlockSizes(tuneKernels)
    , m_mixedPrecisionScheme(makeMatrixStorageMPScheme(mixedPrecisionScheme))
    , m_reorder(reorder)
    , m_diagIdxs(m_gpuMatrix.N())
{
    // TODO: Should in some way verify that this matrix is symmetric, only do it debug mode?
    // Some sanity check
    OPM_ERROR_IF(
        cpuMatrix.N() != m_gpuMatrix.N(),
        fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_gpuMatrix.N(), cpuMatrix.N()));
    OPM_ERROR_IF(cpuMatrix[0][0].N() != m_gpuMatrix.blockSize(),
                 fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.blockSize(),
                             cpuMatrix[0][0].N()));
    OPM_ERROR_IF(cpuMatrix.N() * cpuMatrix[0][0].N() != m_gpuMatrix.dim(),
                 fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.dim(),
                             cpuMatrix.N() * cpuMatrix[0][0].N()));
    OPM_ERROR_IF(cpuMatrix.nonzeroes() != m_gpuMatrix.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_gpuMatrix.nonzeroes(),
                             cpuMatrix.nonzeroes()));
    if (!m_reorder && (m_splitMatrix || m_mixedPrecisionScheme != MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG)) {
        OPM_THROW(std::runtime_error,
                  "Reordering is only supported for full precision matrices without matrix splitting.");
    }
    if (m_reorder) {
        if (m_splitMatrix) {
            m_gpuMatrixReorderedDiag = std::make_unique<GpuVector<field_type>>(blocksize_ * blocksize_ * cpuMatrix.N());
            std::tie(m_gpuMatrixReorderedLower, m_gpuMatrixReorderedUpper)
                = detail::extractLowerAndUpperMatrices<M, field_type,GpuSparseMatrixWrapper<field_type>>(cpuMatrix,
                                                                                                m_reorderedToNatural);
        } else {
            m_gpuMatrixReordered = detail::createReorderedMatrix<M, field_type,GpuSparseMatrixWrapper<field_type>>(
                cpuMatrix, m_reorderedToNatural);
        }

        if (m_mixedPrecisionScheme != MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG) {
            if (!m_splitMatrix) {
                OPM_THROW(std::runtime_error, "Matrix must be split when storing as float.");
            }
            m_gpuMatrixReorderedLowerFloat = std::make_unique<FloatMat>(
                m_gpuMatrixReorderedLower->getRowIndices(), m_gpuMatrixReorderedLower->getColumnIndices(), blocksize_);
            m_gpuMatrixReorderedUpperFloat = std::make_unique<FloatMat>(
                m_gpuMatrixReorderedUpper->getRowIndices(), m_gpuMatrixReorderedUpper->getColumnIndices(), blocksize_);

            if (m_mixedPrecisionScheme == MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG) {
                m_gpuDInvFloat
                    = std::make_unique<FloatVec>(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize());
            }
        }
    }

    // Precompute diagaonal block element indices and handle the actual reordering if that is enabled.
    if (m_reorder) {
        if (!m_splitMatrix) {
            detail::computeDiagIndices<field_type>(m_gpuMatrixReordered->getRowIndices().data(),
                                                m_gpuMatrixReordered->getColumnIndices().data(),
                                                m_gpuReorderToNatural.data(),
                                                m_gpuMatrix.N(),
                                                m_diagIdxs.data());
        }
        reorderAndSplitMatrix(m_moveThreadBlockSize);
    } else {
        detail::computeDiagIndicesNoReorder<field_type>(m_gpuMatrix.getRowIndices().data(),
                                                        m_gpuMatrix.getColumnIndices().data(),
                                                        m_gpuLevelSets.data(),
                                                        m_gpuMatrix.N(),
                                                        m_diagIdxs.data());
    }
    computeDiagonal(m_DILUFactorizationThreadBlockSize);

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
    // ensure that this stream only starts doing work when main stream is completed up to this point
    OPM_GPU_SAFE_CALL(cudaEventRecord(m_before.get(), 0));
    OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(m_stream.get(), m_before.get(), 0));

    OPM_TIMEBLOCK(prec_apply);
    {
        const auto ptrs = std::make_pair(v.data(), d.data());

        auto it = m_apply_graphs.find(ptrs);

        if (it == m_apply_graphs.end()) {
            OPM_GPU_SAFE_CALL(cudaStreamBeginCapture(m_stream.get(), cudaStreamCaptureModeGlobal));

            // The apply functions contains lots of small function calls which call a kernel each
            apply(v, d, m_lowerSolveThreadBlockSize, m_upperSolveThreadBlockSize);

            OPM_GPU_SAFE_CALL(cudaStreamEndCapture(m_stream.get(), &m_apply_graphs[ptrs].get()));
            OPM_GPU_SAFE_CALL(
                cudaGraphInstantiate(&m_executableGraphs[ptrs].get(), m_apply_graphs[ptrs].get(), nullptr, nullptr, 0));
        }
        OPM_GPU_SAFE_CALL(cudaGraphLaunch(m_executableGraphs[ptrs].get(), 0));
    }

    // ensure that main stream only continues after this stream is completed
    OPM_GPU_SAFE_CALL(cudaEventRecord(m_after.get(), m_stream.get()));
    OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(0, m_after.get(), 0));
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::apply(X& v, const Y& d, int lowerSolveThreadBlockSize, int upperSolveThreadBlockSize)
{
    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        if (m_splitMatrix) {
            if (m_mixedPrecisionScheme == MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::solveLowerLevelSetSplit<blocksize_, field_type, float, float>(
                    m_gpuMatrixReorderedLowerFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLowerFloat->getRowIndices().data(),
                    m_gpuMatrixReorderedLowerFloat->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInvFloat->data(),
                    d.data(),
                    v.data(),
                    lowerSolveThreadBlockSize,
                    m_stream.get());
            } else if (m_mixedPrecisionScheme == MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::solveLowerLevelSetSplit<blocksize_, field_type, float, field_type>(
                    m_gpuMatrixReorderedLowerFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLowerFloat->getRowIndices().data(),
                    m_gpuMatrixReorderedLowerFloat->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    d.data(),
                    v.data(),
                    lowerSolveThreadBlockSize,
                    m_stream.get());
            } else {
                detail::DILU::solveLowerLevelSetSplit<blocksize_, field_type, field_type, field_type>(
                    m_gpuMatrixReorderedLower->getNonZeroValues().data(),
                    m_gpuMatrixReorderedLower->getRowIndices().data(),
                    m_gpuMatrixReorderedLower->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    d.data(),
                    v.data(),
                    lowerSolveThreadBlockSize,
                    m_stream.get());
            }
        } else {
            if (m_reorder) {
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
                    lowerSolveThreadBlockSize,
                    m_stream.get());
            } else {
                detail::DILU::solveLowerLevelSetNoReorder<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                                                  m_gpuMatrix.getRowIndices().data(),
                                                                                  m_gpuMatrix.getColumnIndices().data(),
                                                                                  m_gpuLevelSets.data(),
                                                                                  levelStartIdx,
                                                                                  numOfRowsInLevel,
                                                                                  m_gpuDInv.data(),
                                                                                  d.data(),
                                                                                  v.data(),
                                                                                  lowerSolveThreadBlockSize,
                                                                                  m_stream.get());
            }
        }
        levelStartIdx += numOfRowsInLevel;
    }

    levelStartIdx = m_gpuMatrix.N();
    //  upper triangular solve: (D + U_A) v = Dy
    for (int level = m_levelSets.size() - 1; level >= 0; --level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        levelStartIdx -= numOfRowsInLevel;
        if (m_splitMatrix) {
            if (m_mixedPrecisionScheme == MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::solveUpperLevelSetSplit<blocksize_, field_type, float>(
                    m_gpuMatrixReorderedUpperFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpperFloat->getRowIndices().data(),
                    m_gpuMatrixReorderedUpperFloat->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInvFloat->data(),
                    v.data(),
                    upperSolveThreadBlockSize,
                    m_stream.get());
            } else if (m_mixedPrecisionScheme == MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::solveUpperLevelSetSplit<blocksize_, field_type, float>(
                    m_gpuMatrixReorderedUpperFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpperFloat->getRowIndices().data(),
                    m_gpuMatrixReorderedUpperFloat->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    v.data(),
                    upperSolveThreadBlockSize,
                    m_stream.get());
            } else {
                detail::DILU::solveUpperLevelSetSplit<blocksize_, field_type, field_type>(
                    m_gpuMatrixReorderedUpper->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpper->getRowIndices().data(),
                    m_gpuMatrixReorderedUpper->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    v.data(),
                    upperSolveThreadBlockSize,
                    m_stream.get());
            }
        } else {
            if (m_reorder) {
                detail::DILU::solveUpperLevelSet<field_type, blocksize_>(
                    m_gpuMatrixReordered->getNonZeroValues().data(),
                    m_gpuMatrixReordered->getRowIndices().data(),
                    m_gpuMatrixReordered->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    v.data(),
                    upperSolveThreadBlockSize,
                    m_stream.get());
            } else {
                detail::DILU::solveUpperLevelSetNoReorder<field_type, blocksize_>(m_gpuMatrix.getNonZeroValues().data(),
                                                                                  m_gpuMatrix.getRowIndices().data(),
                                                                                  m_gpuMatrix.getColumnIndices().data(),
                                                                                  m_gpuLevelSets.data(),
                                                                                  levelStartIdx,
                                                                                  numOfRowsInLevel,
                                                                                  m_gpuDInv.data(),
                                                                                  v.data(),
                                                                                  upperSolveThreadBlockSize,
                                                                                  m_stream.get());
            }
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
        if (m_reorder) {
            reorderAndSplitMatrix(m_moveThreadBlockSize);
        }
        computeDiagonal(m_DILUFactorizationThreadBlockSize);
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::update(int moveThreadBlockSize, int factorizationBlockSize)
{
    // ensure that this stream only starts doing work when main stream is completed up to this point
    OPM_GPU_SAFE_CALL(cudaEventRecord(m_before.get(), 0));
    OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(m_stream.get(), m_before.get(), 0));

    if (m_reorder) {
        reorderAndSplitMatrix(moveThreadBlockSize);
    }
    computeDiagonal(factorizationBlockSize);

    // ensure that main stream only continues after this stream is completed
    OPM_GPU_SAFE_CALL(cudaEventRecord(m_after.get(), m_stream.get()));
    OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(0, m_after.get(), 0));
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::reorderAndSplitMatrix(int moveThreadBlockSize)
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
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::computeDiagonal(int factorizationBlockSize)
{
    int levelStartIdx = 0;
    for (int level = 0; level < m_levelSets.size(); ++level) {
        const int numOfRowsInLevel = m_levelSets[level].size();
        if (m_splitMatrix) {
            if (m_mixedPrecisionScheme == MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::computeDiluDiagonalSplit<blocksize_,
                                                       field_type,
                                                       float,
                                                       MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG>(
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
                    m_gpuDInvFloat->data(),
                    m_gpuMatrixReorderedLowerFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpperFloat->getNonZeroValues().data(),
                    factorizationBlockSize);
            } else if (m_mixedPrecisionScheme == MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG) {
                detail::DILU::computeDiluDiagonalSplit<blocksize_,
                                                       field_type,
                                                       float,
                                                       MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG>(
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
                    nullptr,
                    m_gpuMatrixReorderedLowerFloat->getNonZeroValues().data(),
                    m_gpuMatrixReorderedUpperFloat->getNonZeroValues().data(),
                    factorizationBlockSize);
            } else {
                // TODO: should this be field type twice or field type then float in the template?
                detail::DILU::computeDiluDiagonalSplit<blocksize_,
                                                       field_type,
                                                       float,
                                                       MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG>(
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
                    nullptr,
                    nullptr,
                    nullptr,
                    factorizationBlockSize);
            }
        } else {
            if (m_reorder) {
                detail::DILU::computeDiluDiagonal<field_type, blocksize_>(
                    m_gpuMatrixReordered->getNonZeroValues().data(),
                    m_gpuMatrixReordered->getRowIndices().data(),
                    m_gpuMatrixReordered->getColumnIndices().data(),
                    m_gpuReorderToNatural.data(),
                    m_gpuNaturalToReorder.data(),
                    m_diagIdxs.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    factorizationBlockSize);
            } else {
                detail::DILU::computeDiluDiagonalNoReorder<field_type, blocksize_>(
                    m_gpuMatrix.getNonZeroValues().data(),
                    m_gpuMatrix.getRowIndices().data(),
                    m_gpuMatrix.getColumnIndices().data(),
                    m_gpuLevelSets.data(),
                    m_diagIdxs.data(),
                    levelStartIdx,
                    numOfRowsInLevel,
                    m_gpuDInv.data(),
                    factorizationBlockSize);
            }
        }
        levelStartIdx += numOfRowsInLevel;
    }
}

template <class M, class X, class Y, int l>
void
GpuDILU<M, X, Y, l>::tuneThreadBlockSizes()
{
    if (m_reorder) {
        // tune the thread-block size of the update function
        auto tuneMoveThreadBlockSizeInUpdate = [this](int moveThreadBlockSize) {
            this->update(moveThreadBlockSize, m_DILUFactorizationThreadBlockSize);
        };
        m_moveThreadBlockSize = detail::tuneThreadBlockSize(tuneMoveThreadBlockSizeInUpdate,
                                                            "(in DILU update) Move data to reordered matrix");
    }

    auto tuneFactorizationThreadBlockSizeInUpdate = [this](int factorizationThreadBlockSize) {
        this->update(m_moveThreadBlockSize, factorizationThreadBlockSize);
    };
    m_DILUFactorizationThreadBlockSize
        = detail::tuneThreadBlockSize(tuneFactorizationThreadBlockSizeInUpdate, "(in DILU update) DILU factorization");

    // tune the thread-block size of the apply
    GpuVector<field_type> tmpV(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    GpuVector<field_type> tmpD(m_gpuMatrix.N() * m_gpuMatrix.blockSize());
    tmpD = 1;

    auto tuneLowerSolveThreadBlockSizeInApply = [this, &tmpV, &tmpD](int lowerSolveThreadBlockSize) {
        this->apply(tmpV, tmpD, lowerSolveThreadBlockSize, m_DILUFactorizationThreadBlockSize);
    };
    m_lowerSolveThreadBlockSize
        = detail::tuneThreadBlockSize(tuneLowerSolveThreadBlockSizeInApply, "(in DILU apply) Triangular lower solve");

    auto tuneUpperSolveThreadBlockSizeInApply = [this, &tmpV, &tmpD](int upperSolveThreadBlockSize) {
        this->apply(tmpV, tmpD, m_lowerSolveThreadBlockSize, upperSolveThreadBlockSize);
    };
    m_upperSolveThreadBlockSize
        = detail::tuneThreadBlockSize(tuneUpperSolveThreadBlockSizeInApply, "(in DILU apply) Triangular upper solve");
}

} // namespace Opm::gpuistl
#define INSTANTIATE_CUDILU_DUNE(realtype, blockdim)                                                                    \
    template class ::Opm::gpuistl::GpuDILU<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,          \
                                           ::Opm::gpuistl::GpuVector<realtype>,                                        \
                                           ::Opm::gpuistl::GpuVector<realtype>>;                                       \
    template class ::Opm::gpuistl::GpuDILU<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,           \
                                           ::Opm::gpuistl::GpuVector<realtype>,                                        \
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
