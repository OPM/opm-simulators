/*
  Copyright 2025 Equinor ASA

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

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixGeneric.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_utilities.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <cuda.h>
#include <fmt/core.h>

#include <cstddef>
#include <memory>

namespace Opm::gpuistl
{

template <class T>
GpuSparseMatrixGeneric<T>::GpuSparseMatrixGeneric(const T* nonZeroElements,
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
    , m_matrixDescriptor(detail::makeSafeMatrixDescriptor())
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    if (detail::to_size_t(rowIndices[numberOfRows]) != numberOfNonzeroBlocks) {
        OPM_THROW(std::invalid_argument, "Wrong sparsity format. Needs to be CSR compliant.");
    }

    initializeMatrixDescriptor();
    preprocessSpMV();
}

template <class T>
GpuSparseMatrixGeneric<T>::GpuSparseMatrixGeneric(const GpuVector<int>& rowIndices,
                                                  const GpuVector<int>& columnIndices,
                                                  size_t blockSize)
    : m_nonZeroElements(columnIndices.dim() * blockSize * blockSize)
    , m_columnIndices(columnIndices)
    , m_rowIndices(rowIndices)
    , m_numberOfNonzeroBlocks(detail::to_int(columnIndices.dim()))
    , m_numberOfRows(detail::to_int(rowIndices.dim() - 1))
    , m_blockSize(detail::to_int(blockSize))
    , m_matrixDescriptor(detail::makeSafeMatrixDescriptor())
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    initializeMatrixDescriptor();
    preprocessSpMV();
}

template <class T>
GpuSparseMatrixGeneric<T>::GpuSparseMatrixGeneric(const GpuSparseMatrixGeneric<T>& other)
    : m_nonZeroElements(other.m_nonZeroElements)
    , m_columnIndices(other.m_columnIndices)
    , m_rowIndices(other.m_rowIndices)
    , m_numberOfNonzeroBlocks(other.m_numberOfNonzeroBlocks)
    , m_numberOfRows(other.m_numberOfRows)
    , m_blockSize(other.m_blockSize)
    , m_matrixDescriptor(detail::makeSafeMatrixDescriptor())
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
{
    initializeMatrixDescriptor();
    preprocessSpMV();
}

template <class T>
void
GpuSparseMatrixGeneric<T>::initializeMatrixDescriptor()
{
    // Create matrix descriptor based on blockSize
    if (m_blockSize > 1) {
#if !USE_HIP && CUDA_VERSION >= 12030
        // Use BSR format for blocked matrices (requires CUDA 12.3+)
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateBsr(m_matrixDescriptor.get(),
                                                 m_numberOfRows,
                                                 m_numberOfRows,
                                                 m_numberOfNonzeroBlocks,
                                                 m_blockSize,
                                                 m_blockSize,
                                                 m_rowIndices.data(),
                                                 m_columnIndices.data(),
                                                 m_nonZeroElements.data(),
                                                 CUSPARSE_INDEX_32I,
                                                 CUSPARSE_INDEX_32I,
                                                 CUSPARSE_INDEX_BASE_ZERO,
                                                 getDataType(),
                                                 CUSPARSE_ORDER_ROW));
#else
        OPM_THROW(std::invalid_argument, "BSR format not supported for HIP or CUDA < 12.3 with Generic API");
#endif
    } else {
        // Use CSR format for scalar matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateCsr(m_matrixDescriptor.get(),
                                                 m_numberOfRows,
                                                 m_numberOfRows,
                                                 m_numberOfNonzeroBlocks,
                                                 m_rowIndices.data(),
                                                 m_columnIndices.data(),
                                                 m_nonZeroElements.data(),
                                                 CUSPARSE_INDEX_32I,
                                                 CUSPARSE_INDEX_32I,
                                                 CUSPARSE_INDEX_BASE_ZERO,
                                                 getDataType()));
    }
}

template <class T>
void
GpuSparseMatrixGeneric<T>::preprocessSpMV()
{
    // Initialize buffer for SpMV operations and preprocess
    size_t vecSize = m_numberOfRows * m_blockSize;

    // Create temporary vectors for preprocessing
    GpuVector<T> tempX(vecSize);
    GpuVector<T> tempY(vecSize);

    // Create vector descriptors with RAII cleanup
    auto tempVecX = detail::makeSafeVectorDescriptor();
    auto tempVecY = detail::makeSafeVectorDescriptor();

    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(tempVecX.get(), vecSize, tempX.data(), getDataType()));
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(tempVecY.get(), vecSize, tempY.data(), getDataType()));

    // Determine buffer size for SpMV
    T alpha = 1.0;
    T beta = 0.0;
    size_t bufferSize = 0;
    // TODO: Check buffersize for different alpha, beta values
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV_bufferSize(m_cusparseHandle.get(),
                                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &alpha,
                                                   *m_matrixDescriptor,
                                                   *tempVecX,
                                                   &beta,
                                                   *tempVecY,
                                                   getDataType(),
                                                   CUSPARSE_SPMV_ALG_DEFAULT,
                                                   &bufferSize));

    m_buffer.resize(bufferSize);

#if CUDA_VERSION >= 12040
    // Preprocess SpMV operation to optimize for this sparsity pattern (requires CUDA 12.4+)
    // TODO: Does the value of the alpha and beta matter in this preprocessing?
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV_preprocess(m_cusparseHandle.get(),
                                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &alpha,
                                                   *m_matrixDescriptor,
                                                   *tempVecX,
                                                   &beta,
                                                   *tempVecY,
                                                   getDataType(),
                                                   CUSPARSE_SPMV_ALG_DEFAULT,
                                                   m_buffer.data()));
#endif

    // Descriptors automatically cleaned up by unique_ptr destructors
}

template <class T>
GpuSparseMatrixGeneric<T>::~GpuSparseMatrixGeneric()
{
    // Matrix descriptor is automatically cleaned up by unique_ptr deleter
}

template <typename T>
template <typename MatrixType>
GpuSparseMatrixGeneric<T>
GpuSparseMatrixGeneric<T>::fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    // TODO: Do we need this intermediate storage? Or this shuffling of data?
    auto [columnIndices, rowIndices] = detail::extractSparsityPattern(matrix);

    constexpr size_t blockSize = MatrixType::block_type::rows;
    const size_t numberOfRows = matrix.N();
    const size_t numberOfNonzeroBlocks = matrix.nonzeroes();

    detail::validateMatrixConversion(matrix, columnIndices, rowIndices);

    const T* nonZeroElements;
    std::vector<T> nonZeroElementsData;

    if (copyNonZeroElementsDirectly) {
        nonZeroElements = detail::findFirstElementPointer<T>(matrix);
        OPM_ERROR_IF(nonZeroElements == nullptr, "error converting DUNE matrix to CuSparse matrix");
    } else {
        nonZeroElementsData = detail::extractNonzeroValues<T>(matrix);
        nonZeroElements = nonZeroElementsData.data();
    }

    return GpuSparseMatrixGeneric<T>(
        nonZeroElements, rowIndices.data(), columnIndices.data(), numberOfNonzeroBlocks, blockSize, numberOfRows);
}

template <class T>
template <class MatrixType>
void
GpuSparseMatrixGeneric<T>::updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    // Validate matrix compatibility
    detail::validateMatrixCompatibility(
        nonzeroes(), blockSize(), N(), matrix.nonzeroes(), matrix[0][0].N(), matrix.N());

    if (copyNonZeroElementsDirectly) {
        const T* newNonZeroElements = static_cast<const T*>(&((matrix[0][0][0][0])));
        m_nonZeroElements.copyFromHost(newNonZeroElements, nonzeroes() * blockSize() * blockSize());
    } else {
        auto nonZeroElementsData = detail::extractNonzeroValues<T>(matrix);
        m_nonZeroElements.copyFromHost(nonZeroElementsData.data(), nonzeroes() * blockSize() * blockSize());
    }
}

template <class T>
void
GpuSparseMatrixGeneric<T>::updateNonzeroValues(const GpuSparseMatrixGeneric<T>& matrix)
{
    detail::validateMatrixCompatibility(
        nonzeroes(), blockSize(), N(), matrix.nonzeroes(), matrix.blockSize(), matrix.N());

    m_nonZeroElements.copyFromDeviceToDevice(matrix.getNonZeroValues());
}



template <typename T>
void
GpuSparseMatrixGeneric<T>::spMV(T alpha, const GpuVector<T>& x, T beta, GpuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);

    // Create vector descriptors with RAII cleanup
    auto vecX = detail::makeSafeVectorDescriptor();
    auto vecY = detail::makeSafeVectorDescriptor();

    // Create dense vector descriptors
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(vecX.get(), x.dim(), const_cast<T*>(x.data()), getDataType()));
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(vecY.get(), y.dim(), y.data(), getDataType()));

    // Execute SpMV operation - buffer already prepared and preprocessed in constructor
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV(m_cusparseHandle.get(),
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &alpha,
                                        *m_matrixDescriptor,
                                        *vecX,
                                        &beta,
                                        *vecY,
                                        getDataType(),
                                        CUSPARSE_SPMV_ALG_DEFAULT,
                                        m_buffer.data()));

    // Descriptors automatically cleaned up by unique_ptr destructors
}

template <typename T>
void
GpuSparseMatrixGeneric<T>::mv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    T alpha = 1.0;
    T beta = 0.0;
    spMV(alpha, x, beta, y);
}

template <typename T>
void
GpuSparseMatrixGeneric<T>::umv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    T alpha = 1.0;
    T beta = 1.0;
    spMV(alpha, x, beta, y);
}

template <typename T>
void
GpuSparseMatrixGeneric<T>::usmv(T alpha, const GpuVector<T>& x, GpuVector<T>& y) const
{
    T beta = 1.0;
    spMV(alpha, x, beta, y);
}

template <class T>
template <class VectorType>
void
GpuSparseMatrixGeneric<T>::assertSameSize(const VectorType& vector) const
{
    // Assume square matrices: numberOfColumns == numberOfRows
    detail::validateVectorMatrixSizes(vector.dim(), blockSize(), N());
}

template class GpuSparseMatrixGeneric<float>;
template class GpuSparseMatrixGeneric<double>;

INSTANTIATE_FOR_TYPE_AND_CLASS(GpuSparseMatrixGeneric, float);
INSTANTIATE_FOR_TYPE_AND_CLASS(GpuSparseMatrixGeneric, double);

} // namespace Opm::gpuistl
