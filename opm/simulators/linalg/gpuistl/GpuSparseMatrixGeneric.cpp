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
#include "opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp"
#include <cuda.h>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <fmt/core.h>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixGeneric.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm::gpuistl
{

namespace
{
    template <class T, class M>
    std::vector<T> extractNonzeroValues(const M& matrix)
    {
        const size_t blockSize = matrix[0][0].N();
        const size_t numberOfNonzeroBlocks = matrix.nonzeroes();
        const size_t numberOfNonzeroElements = blockSize * blockSize * numberOfNonzeroBlocks;

        std::vector<T> nonZeroElementsData;
        nonZeroElementsData.reserve(numberOfNonzeroElements);
        for (auto& row : matrix) {
            for (auto columnIterator = row.begin(); columnIterator != row.end(); ++columnIterator) {
                for (size_t c = 0; c < blockSize; ++c) {
                    for (size_t d = 0; d < blockSize; ++d) {
                        nonZeroElementsData.push_back((*columnIterator)[c][d]);
                    }
                }
            }
        }

        return nonZeroElementsData;
    }
} // namespace

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
    , m_matrixDescriptor(nullptr)
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
    , m_buffer(nullptr)
    , m_bufferSize(0)
{
    if (detail::to_size_t(rowIndices[numberOfRows]) != numberOfNonzeroBlocks) {
        OPM_THROW(std::invalid_argument, "Wrong sparsity format. Needs to be CSR compliant. ");
    }

    // Create matrix descriptor based on blockSize
    if (blockSize > 1) {
        // Use BSR format for blocked matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateBsr(
            &m_matrixDescriptor,
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
            CUSPARSE_ORDER_ROW
        ));
    } else {
        // Use CSR format for scalar matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
            &m_matrixDescriptor,
            m_numberOfRows,
            m_numberOfRows,
            m_numberOfNonzeroBlocks,
            m_rowIndices.data(),
            m_columnIndices.data(),
            m_nonZeroElements.data(),
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO,
            getDataType()
        ));
    }

    // Preprocess SpMV operation
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
    , m_numberOfRows(detail::to_int(rowIndices.dim()-1))
    , m_blockSize(detail::to_int(blockSize))
    , m_matrixDescriptor(nullptr)
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
    , m_buffer(nullptr)
    , m_bufferSize(0)
{
    // Create matrix descriptor based on blockSize
    if (blockSize > 1) {
        // Use BSR format for blocked matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateBsr(
            &m_matrixDescriptor,
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
            CUSPARSE_ORDER_ROW
        ));
    } else {
        // Use CSR format for scalar matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
            &m_matrixDescriptor,
            m_numberOfRows,
            m_numberOfRows,
            m_numberOfNonzeroBlocks,
            m_rowIndices.data(),
            m_columnIndices.data(),
            m_nonZeroElements.data(),
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO,
            getDataType()
        ));
    }

    // Preprocess SpMV operation
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
    , m_matrixDescriptor(nullptr)
    , m_cusparseHandle(detail::CuSparseHandle::getInstance())
    , m_buffer(nullptr)
    , m_bufferSize(0)
{
    // Create a new matrix descriptor with the same properties
    if (other.blockSize() > 1) {
        // Use BSR format for blocked matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateBsr(
            &m_matrixDescriptor,
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
            CUSPARSE_ORDER_ROW
        ));
    } else {
        // Use CSR format for scalar matrices
        OPM_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
            &m_matrixDescriptor,
            m_numberOfRows,
            m_numberOfRows,
            m_numberOfNonzeroBlocks,
            m_rowIndices.data(),
            m_columnIndices.data(),
            m_nonZeroElements.data(),
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO,
            getDataType()
        ));
    }

    // Preprocess SpMV operation
    preprocessSpMV();
}

template <class T>
void GpuSparseMatrixGeneric<T>::preprocessSpMV()
{
    // Initialize buffer for SpMV operations and preprocess
    size_t vecSize = m_numberOfRows * m_blockSize;

    // Create temporary vectors and descriptors for preprocessing
    GpuVector<T> tempX(vecSize);
    GpuVector<T> tempY(vecSize);
    cusparseDnVecDescr_t tempVecX, tempVecY;

    // Create vector descriptors for preprocessing
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
        &tempVecX, vecSize, tempX.data(), getDataType()));
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
        &tempVecY, vecSize, tempY.data(), getDataType()));

    // Determine buffer size for SpMV
    T alpha = 1.0;
    T beta = 0.0;
    // TODO: Check buffersize for different alpha, beta values
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV_bufferSize(
        m_cusparseHandle.get(),
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        m_matrixDescriptor,
        tempVecX,
        &beta,
        tempVecY,
        getDataType(),
        CUSPARSE_SPMV_ALG_DEFAULT,
        &m_bufferSize
    ));

    // Allocate buffer
    OPM_GPU_SAFE_CALL(cudaMalloc(&m_buffer, m_bufferSize));

    // Preprocess SpMV operation to optimize for this sparsity pattern
    // TODO: Does the value of the alpha and beta matter in this preprocessing?
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV_preprocess(
        m_cusparseHandle.get(),
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        m_matrixDescriptor,
        tempVecX,
        &beta,
        tempVecY,
        getDataType(),
        CUSPARSE_SPMV_ALG_DEFAULT,
        m_buffer
    ));

    // Clean up temporary descriptors
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(tempVecX));
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(tempVecY));
}

template <class T>
GpuSparseMatrixGeneric<T>::~GpuSparseMatrixGeneric()
{
    if (m_matrixDescriptor != nullptr) {
        OPM_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(m_matrixDescriptor));
    }

    if (m_buffer != nullptr) {
        OPM_GPU_SAFE_CALL(cudaFree(m_buffer));
    }
}

template <typename T>
template <typename MatrixType>
GpuSparseMatrixGeneric<T>
GpuSparseMatrixGeneric<T>::fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    // TODO: Do we need this intermediate storage? Or this shuffling of data?
    std::vector<int> columnIndices;
    std::vector<int> rowIndices;

    rowIndices.push_back(0);

    // We must find the pointer to the first element in the matrix
    // Iterate until we find an element, we can get the blocksize from the element
    // TODO: Can this be done more cleanly in the DUNE api to access the raw data more directly?
    constexpr size_t blockSizeTmp = MatrixType::block_type::rows;
    T* nonZeroElementsTmp = nullptr;
    for (auto rowIt = matrix.begin(); rowIt != matrix.end(); ++rowIt){
        auto colIt = rowIt->begin();
        if (colIt != rowIt->end()){
            nonZeroElementsTmp = const_cast<T*>(&((*colIt)[0][0]));
            break;
        }
    }

    OPM_ERROR_IF(nonZeroElementsTmp == nullptr, "error converting DUNE matrix to CuSparse matrix");

    const size_t blockSize = blockSizeTmp;
    const size_t numberOfRows = matrix.N();
    const size_t numberOfNonzeroBlocks = matrix.nonzeroes();

    // Extract column indices and row pointers
    columnIndices.reserve(numberOfNonzeroBlocks);
    rowIndices.reserve(numberOfRows + 1);
    for (auto& row : matrix) {
        for (auto columnIterator = row.begin(); columnIterator != row.end(); ++columnIterator) {
            columnIndices.push_back(columnIterator.index());
        }
        rowIndices.push_back(detail::to_int(columnIndices.size()));
    }

    // Sanity check
    // h_rows and h_cols could be changed to 'unsigned int', but cusparse expects 'int'
    OPM_ERROR_IF(rowIndices[matrix.N()] != detail::to_int(matrix.nonzeroes()),
                 "Error size of rows do not sum to number of nonzeroes in GpuSparseMatrixGeneric.");
    OPM_ERROR_IF(rowIndices.size() != numberOfRows + 1, "Row indices do not match for GpuSparseMatrixGeneric.");
    OPM_ERROR_IF(columnIndices.size() != numberOfNonzeroBlocks, "Column indices do not match for GpuSparseMatrixGeneric.");

    if (copyNonZeroElementsDirectly) {
        const T* nonZeroElements = nonZeroElementsTmp;
        return GpuSparseMatrixGeneric<T>(
            nonZeroElements, rowIndices.data(), columnIndices.data(), numberOfNonzeroBlocks, blockSize, numberOfRows);
    } else {
        auto nonZeroElementData = extractNonzeroValues<T>(matrix);
        return GpuSparseMatrixGeneric<T>(nonZeroElementData.data(),
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
GpuSparseMatrixGeneric<T>::updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    OPM_ERROR_IF(nonzeroes() != matrix.nonzeroes(), "Matrix does not have the same number of non-zero elements.");
    OPM_ERROR_IF(matrix[0][0].N() != blockSize(), "Matrix does not have the same blocksize.");
    OPM_ERROR_IF(matrix.N() != N(), "Matrix does not have the same number of rows.");

    if (!copyNonZeroElementsDirectly) {
        auto nonZeroElementsData = extractNonzeroValues<T>(matrix);
        m_nonZeroElements.copyFromHost(nonZeroElementsData.data(), nonzeroes() * blockSize() * blockSize());
    } else {
        const T* newNonZeroElements = static_cast<const T*>(&((matrix[0][0][0][0])));
        m_nonZeroElements.copyFromHost(newNonZeroElements, nonzeroes() * blockSize() * blockSize());
    }

    // Update the matrix values in the descriptor
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMatSetValues(
        m_matrixDescriptor,
        m_nonZeroElements.data()
    ));
}

template <class T>
void
GpuSparseMatrixGeneric<T>::updateNonzeroValues(const GpuSparseMatrixGeneric<T>& matrix)
{
    OPM_ERROR_IF(nonzeroes() != matrix.nonzeroes(), "Matrix does not have the same number of non-zero elements.");
    OPM_ERROR_IF(matrix.N() != N(), "Matrix does not have the same number of rows.");

    m_nonZeroElements.copyFromDeviceToDevice(matrix.getNonZeroValues());

    // Update the matrix values in the descriptor
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMatSetValues(
        m_matrixDescriptor,
        m_nonZeroElements.data()
    ));
}

template <typename T>
void GpuSparseMatrixGeneric<T>::createVectorDescriptors(
    const GpuVector<T>& x, GpuVector<T>& y,
    cusparseDnVecDescr_t& vecX, cusparseDnVecDescr_t& vecY) const 
{
    // Create dense vector descriptors
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
        &vecX,
        x.dim(),
        const_cast<T*>(x.data()),
        getDataType()
    ));
    OPM_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
        &vecY,
        y.dim(),
        y.data(),
        getDataType()
    ));
}

template <typename T>
void GpuSparseMatrixGeneric<T>::spMV(
    T alpha, const GpuVector<T>& x, T beta, GpuVector<T>& y) const 
{
    assertSameSize(x);
    assertSameSize(y);

    // Create vector descriptors
    cusparseDnVecDescr_t vecX, vecY;
    createVectorDescriptors(x, y, vecX, vecY);

    // Execute SpMV operation - buffer already prepared and preprocessed in constructor
    OPM_CUSPARSE_SAFE_CALL(cusparseSpMV(
        m_cusparseHandle.get(),
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        m_matrixDescriptor,
        vecX,
        &beta,
        vecY,
        getDataType(),
        CUSPARSE_SPMV_ALG_DEFAULT,
        m_buffer
    ));

    // Destroy vector descriptors
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY));
}

template <typename T>
void GpuSparseMatrixGeneric<T>::mv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    T alpha = 1.0;
    T beta = 0.0;
    spMV(alpha, x, beta, y);
}

template <typename T>
void GpuSparseMatrixGeneric<T>::umv(const GpuVector<T>& x, GpuVector<T>& y) const
{
    T alpha = 1.0;
    T beta = 1.0;
    spMV(alpha, x, beta, y);
}

template <typename T>
void GpuSparseMatrixGeneric<T>::usmv(T alpha, const GpuVector<T>& x, GpuVector<T>& y) const
{
    T beta = 1.0;
    spMV(alpha, x, beta, y);
}

template <class T>
template <class VectorType>
void GpuSparseMatrixGeneric<T>::assertSameSize(const VectorType& vector) const
{
    if (vector.dim() != blockSize() * N()) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Size mismatch. Input vector has {} elements, while we have {} rows.",
                              vector.dim(),
                              blockSize() * N()));
    }
}

#define INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(realtype, blockdim)                                     \
    template GpuSparseMatrixGeneric<realtype> GpuSparseMatrixGeneric<realtype>::fromMatrix(                            \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template GpuSparseMatrixGeneric<realtype> GpuSparseMatrixGeneric<realtype>::fromMatrix(                            \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool);                                \
    template void GpuSparseMatrixGeneric<realtype>::updateNonzeroValues(                                               \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template void GpuSparseMatrixGeneric<realtype>::updateNonzeroValues(                                               \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool)

template class GpuSparseMatrixGeneric<float>;
template class GpuSparseMatrixGeneric<double>;

INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 1);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 2);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 3);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 4);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 5);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(double, 6);

INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 1);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 2);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 3);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 4);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 5);
INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(float, 6);

} // namespace Opm::gpuistl
