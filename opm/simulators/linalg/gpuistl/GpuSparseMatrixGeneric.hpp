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
#ifndef OPM_GPUSPARSEMATRIXGENERIC_HPP
#define OPM_GPUSPARSEMATRIXGENERIC_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_utilities.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <cstddef>
#include <cusparse.h>

namespace Opm::gpuistl
{

/**
 * @brief The GpuSparseMatrixGeneric class uses cuSPARSE Generic API for sparse matrix operations.
 *
 * @note We support raw primitives for T (double and float). Block size is handled through the
 * blockSize parameter.
 *
 * @tparam T the type to store. Can be either float or double.
 *
 * @note we only support square matrices.
 */
template <typename T>
class GpuSparseMatrixGeneric
{
public:
    using field_type = T;

    /**
    * @brief Maximum block size supported by this implementation.
    *
    * This constant defines an upper bound on the block size to ensure reasonable compilation times.
    * While this class itself could support larger values, functions that call dispatchOnBlocksize()
    * might have limitations. This value can be increased if needed, but will increase compilation time
    * due to template instantiations.
    */
    static constexpr int max_block_size = 6;

    //! Create the sparse matrix specified by the raw data.
    //!
    //! \param[in] nonZeroElements the non-zero values of the matrix
    //! \param[in] rowIndices      the row indices of the non-zero elements
    //! \param[in] columnIndices   the column indices of the non-zero elements
    //! \param[in] numberOfNonzeroBlocks number of nonzero elements
    //! \param[in] blockSize size of each block matrix (typically 3)
    //! \param[in] numberOfRows the number of rows
    //!
    GpuSparseMatrixGeneric(const T* nonZeroElements,
                           const int* rowIndices,
                           const int* columnIndices,
                           std::size_t numberOfNonzeroBlocks,
                           std::size_t blockSize,
                           std::size_t numberOfRows);

    //! Create a sparse matrix by copying the sparsity structure of another matrix, not filling in the values
    //!
    //! \param[in] rowIndices      the row indices of the non-zero elements
    //! \param[in] columnIndices   the column indices of the non-zero elements
    //! \param[in] blockSize size of each block matrix (typically 3)
    //!
    GpuSparseMatrixGeneric(const GpuVector<int>& rowIndices, const GpuVector<int>& columnIndices, std::size_t blockSize);

    GpuSparseMatrixGeneric(const GpuSparseMatrixGeneric&);

    /**
     * @brief Preprocess SpMV operation to optimize for sparsity pattern
     *
     * This function preprocesses the sparsity pattern of the matrix to optimize for the SpMV operation.
     */
    void preprocessSpMV();

    // We want to have this as non-mutable as possible, that is we do not want
    // to deal with changing matrix sizes and sparsity patterns.
    GpuSparseMatrixGeneric& operator=(const GpuSparseMatrixGeneric&) = delete;

    virtual ~GpuSparseMatrixGeneric();

    /**
     * @brief fromMatrix creates a new matrix with the same block size and values as the given matrix
     * @param matrix the matrix to copy from
     * @param copyNonZeroElementsDirectly if true will do a memcpy from matrix[0][0][0][0], otherwise will build up the
     * non-zero elements by looping over the matrix. Note that setting this to true will yield a performance increase,
     * but might not always yield correct results depending on how the matrix has been initialized. If unsure,
     * leave it as false.
     * @tparam MatrixType is assumed to be a Dune::BCRSMatrix compatible matrix.
     */
    template <class MatrixType>
    static GpuSparseMatrixGeneric<T> fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly = false);

    /**
     * @brief N returns the number of rows (which is equal to the number of columns)
     */
    std::size_t N() const
    {
        return detail::to_size_t(m_numberOfRows);
    }

    /**
     * @brief nonzeroes behaves as the Dune::BCRSMatrix::nonzeros() function and returns the number of non zero blocks
     * @return number of non zero blocks.
     */
    std::size_t nonzeroes() const
    {
        // Technically this safe conversion is not needed since we enforce these to be
        // non-negative in the constructor, but keeping them for added sanity for now.
        //
        // We don't believe this will yield any performance penality (it's used too far away from the inner loop),
        // but should that be false, they can be removed.
        return detail::to_size_t(m_numberOfNonzeroBlocks);
    }

    /**
     * @brief getNonZeroValues returns the GPU vector containing the non-zero values (ordered by block)
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<T>& getNonZeroValues()
    {
        return m_nonZeroElements;
    }

    /**
     * @brief getNonZeroValues returns the GPU vector containing the non-zero values (ordered by block)
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<T>& getNonZeroValues() const
    {
        return m_nonZeroElements;
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<int>& getRowIndices()
    {
        return m_rowIndices;
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<int>& getRowIndices() const
    {
        return m_rowIndices;
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<int>& getColumnIndices()
    {
        return m_columnIndices;
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<int>& getColumnIndices() const
    {
        return m_columnIndices;
    }

    /**
     * @brief dim returns the dimension of the vector space on which this matrix acts
     *
     * This is equivalent to matrix.N() * matrix.blockSize()
     * @return matrix.N() * matrix.blockSize()
     */
    std::size_t dim() const
    {
        // Technically this safe conversion is not needed since we enforce these to be
        // non-negative in the constructor, but keeping them for added sanity for now.
        //
        // We don't believe this will yield any performance penality (it's used too far away from the inner loop),
        // but should that be false, they can be removed.
        return detail::to_size_t(m_blockSize) * detail::to_size_t(m_numberOfRows);
    }

    /**
     * @brief blockSize size of the blocks
     */
    std::size_t blockSize() const
    {
        // Technically this safe conversion is not needed since we enforce these to be
        // non-negative in the constructor, but keeping them for added sanity for now.
        //
        // We don't believe this will yield any performance penality (it's used too far away from the inner loop),
        // but should that be false, they can be removed.
        return detail::to_size_t(m_blockSize);
    }

    /**
     * @brief mv performs matrix vector multiply y = Ax
     * @param[in] x the vector to multiply the matrix with
     * @param[out] y the output vector
     */
    virtual void mv(const GpuVector<T>& x, GpuVector<T>& y) const;

    /**
     * @brief umv computes y=Ax+y
     * @param[in] x the vector to multiply with A
     * @param[inout] y the vector to add and store the output in
     */
    virtual void umv(const GpuVector<T>& x, GpuVector<T>& y) const;


    /**
     * @brief umv computes y=alpha * Ax + y
     * @param[in] alpha The scaling factor for the matrix-vector product
     * @param[in] x the vector to multiply with A
     * @param[inout] y the vector to add and store the output in
     */
    virtual void usmv(T alpha, const GpuVector<T>& x, GpuVector<T>& y) const;

    /**
     * @brief updateNonzeroValues updates the non-zero values by using the non-zero values of the supplied matrix
     * @param matrix the matrix to extract the non-zero values from
     * @param copyNonZeroElementsDirectly if true will do a memcpy from matrix[0][0][0][0], otherwise will build up the
     * non-zero elements by looping over the matrix. Note that setting this to true will yield a performance increase,
     * but might not always yield correct results depending on how the matrix matrix has been initialized. If unsure,
     * leave it as false.
     * @note This assumes the given matrix has the same sparsity pattern.
     * @tparam MatrixType is assumed to be a Dune::BCRSMatrix compatible matrix.
     */
    template <class MatrixType>
    void updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly = false);

    /**
     * @brief updateNonzeroValues updates the non-zero values by using the non-zero values of the supplied matrix
     * @param matrix the matrix to extract the non-zero values from
     * @note This assumes the given matrix has the same sparsity pattern.
     */
    void updateNonzeroValues(const GpuSparseMatrixGeneric<T>& matrix);

    /**
     * @brief setToZero resets the matrix to zero values.
     */
    void setToZero();

    /**
     * @brief Dispatches a function based on the block size of the matrix.
     *
     * This method allows executing different code paths depending on the block size
     * of the matrix, up to the maximum block size specified by max_block_size.
     *
     * Use this function if you need the block size to be known at compile time.
     *
     * @tparam FunctionType Type of the function to be dispatched
     * @param function The function to be executed based on the block size
     * @return The result of the function execution
     *
     * You can use this function as
     *
     * \code{.cpp}
     * matrix.dispatchOnBlocksize([](auto val) {
     *    constexpr int blockSize = decltype(val)::value;
     * });
     * \endcode
     */
     template<class FunctionType>
     auto dispatchOnBlocksize(FunctionType function) const
     {
        return dispatchOnBlocksizeImpl<max_block_size>(function);
     }

private:
    GpuVector<T> m_nonZeroElements;
    GpuVector<int> m_columnIndices;
    GpuVector<int> m_rowIndices;

    // Notice that we store these three as int to make sure we are cusparse compatible.
    //
    // This gives the added benefit of checking the size constraints at construction of the matrix
    // rather than in some call to cusparse.
    const int m_numberOfNonzeroBlocks;
    const int m_numberOfRows;
    const int m_blockSize;

    // Generic API descriptors
    decltype(detail::makeSafeMatrixDescriptor()) m_matrixDescriptor;
    detail::CuSparseHandle& m_cusparseHandle;

    // Cached buffer for operations
    mutable GpuBuffer<std::byte> m_buffer;

    // Helper methods for SpMV operations
    void spMV(T alpha, const GpuVector<T>& x, T beta, GpuVector<T>& y) const;

    // Initialize matrix descriptor based on block size
    void initializeMatrixDescriptor();

    template <class VectorType>
    void assertSameSize(const VectorType& vector) const;

    // Helper to get cuSPARSE data type from C++ type
    constexpr cudaDataType getDataType() const
    {
        if constexpr (std::is_same_v<T, float>) {
            return CUDA_R_32F;
        } else if constexpr (std::is_same_v<T, double>) {
            return CUDA_R_64F;
        } else {
            static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Only float and double are supported");
            return CUDA_R_32F; // Unreachable, but needed to compile
        }
    }

    template<int blockSizeCompileTime, class FunctionType>
    auto dispatchOnBlocksizeImpl(FunctionType function) const
    {
        if (blockSizeCompileTime == m_blockSize) {
            return function(std::integral_constant<int, blockSizeCompileTime>());
        }

        if constexpr (blockSizeCompileTime > 1) {
            return dispatchOnBlocksizeImpl<blockSizeCompileTime - 1>(function);
        } else {
            OPM_THROW(std::runtime_error, fmt::format(fmt::runtime("Unsupported block size: {}"), m_blockSize));
        }
    }
};
} // namespace Opm::gpuistl
#endif
