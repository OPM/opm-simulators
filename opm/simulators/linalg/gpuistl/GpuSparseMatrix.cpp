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
#include <config.h>

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixGeneric.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_constants.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_wrapper.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_utilities.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <cuda.h>
#include <fmt/core.h>
namespace Opm::gpuistl
{

template <class T>
GpuSparseMatrix<T>::GpuSparseMatrix(const T* nonZeroElements,
                                  const int* rowIndices,
                                  const int* columnIndices,
                                  size_t numberOfNonzeroBlocks,
                                  size_t blockSize,
                                  size_t numberOfRows)
    : m_nonZeroElements(nonZeroElements, numberOfNonzeroBlocks * blockSize * blockSize)
    , m_columnIndices(columnIndices, numberOfNonzeroBlocks)
    , m_rowIndices(rowIndices, numberOfRows + 1)
    , m_numberOfNonzeroBlocks(detail::to_int(numberOfNonzeroBlocks))
    , m_numberOfRows(detail::to_int(numberOfRows))
    , m_blockSize(detail::to_int(blockSize))
    , m_matrixDescription(detail::createMatrixDescription())
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    if (detail::to_size_t(rowIndices[numberOfRows]) != numberOfNonzeroBlocks) {
        OPM_THROW(std::invalid_argument, "Wrong sparsity format. Needs to be CSR compliant. ");
    }

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (blockSize == 1) {
        m_genericMatrixForBlockSize1 = std::make_unique<GpuSparseMatrixGeneric<T>>(
            nonZeroElements, rowIndices, columnIndices, numberOfNonzeroBlocks, blockSize, numberOfRows);
    }
}

template <class T>
GpuSparseMatrix<T>::GpuSparseMatrix(const GpuVector<int>& rowIndices,
                                  const GpuVector<int>& columnIndices,
                                  size_t blockSize)
    : m_nonZeroElements(columnIndices.dim() * blockSize * blockSize)
    , m_columnIndices(columnIndices)
    , m_rowIndices(rowIndices)
    , m_numberOfNonzeroBlocks(detail::to_int(columnIndices.dim()))
    , m_numberOfRows(detail::to_int(rowIndices.dim()-1))
    , m_blockSize(detail::to_int(blockSize))
    , m_matrixDescription(detail::createMatrixDescription())
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (blockSize == 1) {
        m_genericMatrixForBlockSize1 = std::make_unique<GpuSparseMatrixGeneric<T>>(
            rowIndices, columnIndices, blockSize);
    }
}

template<class T>
GpuSparseMatrix<T>::GpuSparseMatrix(const GpuSparseMatrix<T>& other)
    : m_nonZeroElements(other.m_nonZeroElements)
    , m_columnIndices(other.m_columnIndices)
    , m_rowIndices(other.m_rowIndices)
    , m_numberOfNonzeroBlocks(other.m_numberOfNonzeroBlocks)
    , m_numberOfRows(other.m_numberOfRows)
    , m_blockSize(other.m_blockSize)
    , m_matrixDescription(other.m_matrixDescription)
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (other.blockSize() == 1) {
        OPM_ERROR_IF(!other.m_genericMatrixForBlockSize1,
            "Internal error, other.blockSize() is 1 but other.m_genericMatrixForBlockSize1 has not been set.");
        m_genericMatrixForBlockSize1 = std::make_unique<GpuSparseMatrixGeneric<T>>(*other.m_genericMatrixForBlockSize1);
    }
}


template <class T>
GpuSparseMatrix<T>::~GpuSparseMatrix()
{
    // empty
}

template <typename T>
template <typename MatrixType>
GpuSparseMatrix<T>
GpuSparseMatrix<T>::fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    // TODO: Do we need this intermediate storage? Or this shuffling of data?
    auto [columnIndices, rowIndices] = detail::extractSparsityPattern(matrix);

    constexpr size_t blockSize = MatrixType::block_type::rows;
    const size_t numberOfRows = matrix.N();
    const size_t numberOfNonzeroBlocks = matrix.nonzeroes();

    detail::validateMatrixConversion(matrix, columnIndices, rowIndices);

    if (copyNonZeroElementsDirectly) {
        // Find pointer to first element for direct copying
        T* nonZeroElementsTmp = detail::findFirstElementPointer<T>(matrix);
        OPM_ERROR_IF(!nonZeroElementsTmp, "error converting DUNE matrix to CuSparse matrix");

        const T* nonZeroElements = nonZeroElementsTmp;
        return GpuSparseMatrix<T>(
            nonZeroElements, rowIndices.data(), columnIndices.data(), numberOfNonzeroBlocks, blockSize, numberOfRows);
    } else {
        auto nonZeroElementData = detail::extractNonzeroValues<T>(matrix);
        return GpuSparseMatrix<T>(nonZeroElementData.data(),
                                 rowIndices.data(),
                                 columnIndices.data(),
                                 numberOfNonzeroBlocks,
                                 blockSize,
                                 numberOfRows);
    }
}

template <class T>
template <class MatrixType>
void
GpuSparseMatrix<T>::updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    detail::validateMatrixCompatibility(nonzeroes(), blockSize(), N(),
                                       matrix.nonzeroes(), matrix[0][0].N(), matrix.N());

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (m_genericMatrixForBlockSize1) {
        m_genericMatrixForBlockSize1->updateNonzeroValues(matrix, copyNonZeroElementsDirectly);
        return;
    }

    if (!copyNonZeroElementsDirectly) {
        auto nonZeroElementsData = detail::extractNonzeroValues<T>(matrix);
        m_nonZeroElements.copyFromHost(nonZeroElementsData.data(), nonzeroes() * blockSize() * blockSize());

    } else {
        const T* newNonZeroElements = static_cast<const T*>(&((matrix[0][0][0][0])));
        // copy from host to device using default stream and asynchronous transfer
        m_nonZeroElements.copyFromHostAsync(newNonZeroElements, nonzeroes() * blockSize() * blockSize());
    }
}

template <class T>
void
GpuSparseMatrix<T>::updateNonzeroValues(const GpuSparseMatrix<T>& matrix)
{
    detail::validateMatrixCompatibility(nonzeroes(), blockSize(), N(),
                                       matrix.nonzeroes(), matrix.blockSize(), matrix.N());

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (m_genericMatrixForBlockSize1) {
        OPM_ERROR_IF(!matrix.m_genericMatrixForBlockSize1,
            "Internal error, matrix.blockSize() is 1 but matrix.m_genericMatrixForBlockSize1 has not been set.");
        m_genericMatrixForBlockSize1->updateNonzeroValues(*matrix.m_genericMatrixForBlockSize1);
        return;
    }

    m_nonZeroElements.copyFromDeviceToDevice(matrix.getNonZeroValues());
}

template <typename T>
void
GpuSparseMatrix<T>::setUpperTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(m_matrixDescription->get(), CUSPARSE_FILL_MODE_UPPER));
}

template <typename T>
void
GpuSparseMatrix<T>::setLowerTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(m_matrixDescription->get(), CUSPARSE_FILL_MODE_LOWER));
}

template <typename T>
void
GpuSparseMatrix<T>::setUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(m_matrixDescription->get(), CUSPARSE_DIAG_TYPE_UNIT));
}

template <typename T>
void
GpuSparseMatrix<T>::setNonUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(m_matrixDescription->get(), CUSPARSE_DIAG_TYPE_NON_UNIT));
}

template <typename T>
void
GpuSparseMatrix<T>::mv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (m_genericMatrixForBlockSize1) {
        m_genericMatrixForBlockSize1->mv(x, y);
        return;
    }
    const auto nonzeroValues = getNonZeroValues().data();

    auto rowIndices = getRowIndices().data();
    auto columnIndices = getColumnIndices().data();
    T alpha = 1.0;
    T beta = 0.0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrmv(m_cusparseHandle.get(),
                                                 detail::CUSPARSE_MATRIX_ORDER,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 m_numberOfRows,
                                                 m_numberOfRows,
                                                 m_numberOfNonzeroBlocks,
                                                 &alpha,
                                                 m_matrixDescription->get(),
                                                 nonzeroValues,
                                                 rowIndices,
                                                 columnIndices,
                                                 blockSize(),
                                                 x.data(),
                                                 &beta,
                                                 y.data()));
}

template <typename T>
void
GpuSparseMatrix<T>::umv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (m_genericMatrixForBlockSize1) {
        m_genericMatrixForBlockSize1->umv(x, y);
        return;
    }

    const auto nonzeroValues = getNonZeroValues().data();

    auto rowIndices = getRowIndices().data();
    auto columnIndices = getColumnIndices().data();
    T alpha = 1.0;
    T beta = 1.0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrmv(m_cusparseHandle.get(),
                                                 detail::CUSPARSE_MATRIX_ORDER,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 m_numberOfRows,
                                                 m_numberOfRows,
                                                 m_numberOfNonzeroBlocks,
                                                 &alpha,
                                                 m_matrixDescription->get(),
                                                 nonzeroValues,
                                                 rowIndices,
                                                 columnIndices,
                                                 m_blockSize,
                                                 x.data(),
                                                 &beta,
                                                 y.data()));
}

template <typename T>
void
GpuSparseMatrix<T>::usmv(T alpha, const GpuVector<T>& x, GpuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);

    // For blockSize == 1, use GpuSparseMatrixGeneric
    if (m_genericMatrixForBlockSize1) {
        m_genericMatrixForBlockSize1->usmv(alpha, x, y);
        return;
    }
    const auto numberOfRows = N();
    const auto numberOfNonzeroBlocks = nonzeroes();
    const auto nonzeroValues = getNonZeroValues().data();

    auto rowIndices = getRowIndices().data();
    auto columnIndices = getColumnIndices().data();

    T beta = 1.0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrmv(m_cusparseHandle.get(),
                                                 detail::CUSPARSE_MATRIX_ORDER,
                                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 numberOfRows,
                                                 numberOfRows,
                                                 numberOfNonzeroBlocks,
                                                 &alpha,
                                                 m_matrixDescription->get(),
                                                 nonzeroValues,
                                                 rowIndices,
                                                 columnIndices,
                                                 blockSize(),
                                                 x.data(),
                                                 &beta,
                                                 y.data()));
}

template <class T>
template <class VectorType>
void
GpuSparseMatrix<T>::assertSameSize(const VectorType& x) const
{
    // Assume square matrices: numberOfColumns == numberOfRows
    detail::validateVectorMatrixSizes(x.dim(), blockSize(), N());
}



template class GpuSparseMatrix<float>;
template class GpuSparseMatrix<double>;

INSTANTIATE_FOR_TYPE_AND_CLASS(GpuSparseMatrix, float);
INSTANTIATE_FOR_TYPE_AND_CLASS(GpuSparseMatrix, double);

} // namespace Opm::gpuistl
