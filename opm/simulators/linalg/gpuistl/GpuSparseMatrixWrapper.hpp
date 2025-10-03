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
#ifndef OPM_GPUSPARSEMATRIXWRAPPER_HPP
#define OPM_GPUSPARSEMATRIXWRAPPER_HPP
#include <cusparse.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixGeneric.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/gpuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#include <vector>

namespace Opm::gpuistl
{

/**
 * @brief The GpuSparseMatrixWrapper Checks CUDA/HIP version and dispatches a version either using the old or the generic CUDA API
 *
 * @note we currently only support simple raw primitives for T (double and float). Block size is handled through the
 * block size parameter
 *
 * @tparam T the type to store. Can be either float, double or int.
 *
 * @note we only support square matrices.
 *
 * @note We only support Block Compressed Sparse Row Format (BSR) for now.

 * @note This class uses the legacy cuSPARSE API, to be compatible with CuSparse's ilu0 preconditioner. However,
 *       this preconditioner is deprecated and will be removed in future versions of CuSparse. So we should migrate
 *       to the new cuSPARSE generic API in the future.
 *
 * @note To also support block size 1, we use the GpuSparseMatrixGeneric class which uses the new cuSPARSE generic API.
 *       This is a temporary solution, and we should migrate to the new API for all block sizes in the future by
 *       replacing this class with GpuSparseMatrixGeneric.
 */
template <typename T>
class GpuSparseMatrixWrapper

{
public:
    using field_type = T;

    /*
        Here is the primary function of this class.
        Since the generic API for CUDA/HIP is primaryly supported on CUDA 13 (and not yet HIP) for blocked
        matrices, places wanting to use blocked matrices can invoke this class which handles which API to use.
        Basically we just check if HIP is present, or if we are using cuda and a version prior to 13.
    */
#if USE_HIP || (!USE_HIP && CUDA_VERSION < 13000)
     using matrix_type = GpuSparseMatrix<T>;
#else
     using matrix_type = GpuSparseMatrixGeneric<T>;
#endif

    // Arrow operator overloads for direct access to the underlying matrix
    matrix_type* operator->() {
        if (!m_matrix) {
            throw std::runtime_error("GpuSparseMatrixWrapper: underlying matrix is nullptr.");
        }
        return m_matrix.get();
    }
    const matrix_type* operator->() const {
        if (!m_matrix) {
            throw std::runtime_error("GpuSparseMatrixWrapper: underlying matrix is nullptr.");
        }
        return m_matrix.get();
    }

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
    //! \note Prefer to use the constructor taking a const reference to a matrix instead.
    //!
    //! \param[in] nonZeroElements the non-zero values of the matrix
    //! \param[in] rowIndices      the row indices of the non-zero elements
    //! \param[in] columnIndices   the column indices of the non-zero elements
    //! \param[in] numberOfNonzeroBlocks number of nonzero blocks
    //! \param[in] blockSize size of each block matrix (typically 3)
    //! \param[in] numberOfRows the number of rows
    //!
    //! \note We assume numberOfNonzeroBlocks, blockSize and numberOfRows all are representable as int due to
    //!       restrictions in the current version of cusparse. This might change in future versions.
    GpuSparseMatrixWrapper(const T* nonZeroElements,
                   const int* rowIndices,
                   const int* columnIndices,
                   size_t numberOfNonzeroBlocks,
                   size_t blockSize,
                   size_t numberOfRows)
    {
        m_matrix = std::make_unique<matrix_type>(nonZeroElements,
                                                 rowIndices,
                                                 columnIndices,
                                                 numberOfNonzeroBlocks,
                                                 blockSize,
                                                 numberOfRows);
    }

    //! Create a sparse matrix by copying the sparsity structure of another matrix, not filling in the values
    //!
    //! \note Prefer to use the constructor taking a const reference to a matrix instead.
    //!
    //! \param[in] rowIndices      the row indices of the non-zero elements
    //! \param[in] columnIndices   the column indices of the non-zero elements
    //! \param[in] blockSize size of each block matrix (typically 3)
    //!
    //! \note We assume numberOfNonzeroBlocks, blockSize and numberOfRows all are representable as int due to
    //!       restrictions in the current version of cusparse. This might change in future versions.
    GpuSparseMatrixWrapper(const GpuVector<int>& rowIndices,
                   const GpuVector<int>& columnIndices,
                   size_t blockSize)
    {
        m_matrix = std::make_unique<matrix_type>(rowIndices, columnIndices, blockSize);
    }

    GpuSparseMatrixWrapper(const GpuSparseMatrixWrapper& other)
    {
        if (!other.m_matrix) {
            throw std::runtime_error("Internal error, other.m_matrix is a nullptr.");
        }
        m_matrix = std::make_unique<matrix_type>(*other.m_matrix);
    }

    // We want to have this as non-mutable as possible, that is we do not want
    // to deal with changing matrix sizes and sparsity patterns.
    GpuSparseMatrixWrapper& operator=(const GpuSparseMatrixWrapper&) = delete;

    ~GpuSparseMatrixWrapper() = default;

    GpuSparseMatrixWrapper() = default;

    const matrix_type& get() const {
        if (!m_matrix) {
            throw std::runtime_error("GpuSparseMatrixWrapper: underlying matrix is nullptr.");
        }
        return *m_matrix;
    }

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
    static GpuSparseMatrixWrapper<T> fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly = false)
    {
        GpuSparseMatrixWrapper<T> gpuSparseMatrixWrapper;
        gpuSparseMatrixWrapper.m_matrix = std::make_unique<matrix_type>(
            matrix_type::fromMatrix(matrix, copyNonZeroElementsDirectly));
        return gpuSparseMatrixWrapper;
    }

    // Only participates in overload resolution when matrix_type == GpuSparseMatrix<T>
    template <class M = matrix_type,
              typename = std::enable_if_t<std::is_same_v<M, GpuSparseMatrix<T>>>>
    void setUpperTriangular()
    {
        m_matrix->setUpperTriangular();
    }

    /**
     * @brief setLowerTriangular sets the CuSparse flag that this is an lower diagonal (with non-unit diagonal) matrix.
     */
    template <class M = matrix_type,
              typename = std::enable_if_t<std::is_same_v<M, GpuSparseMatrix<T>>>>
    void setLowerTriangular()
    {
        m_matrix->setLowerTriangular();
    }

    /**
     * @brief setUnitDiagonal sets the CuSparse flag that this has unit diagional.
     */
    template <class M = matrix_type,
              typename = std::enable_if_t<std::is_same_v<M, GpuSparseMatrix<T>>>>
    void setUnitDiagonal()
    {
        m_matrix->setUnitDiagonal();
    }

    /**
     * @brief setNonUnitDiagonal sets the CuSparse flag that this has non-unit diagional.
     */
    template <class M = matrix_type,
              typename = std::enable_if_t<std::is_same_v<M, GpuSparseMatrix<T>>>>
    void setNonUnitDiagonal()
    {
        m_matrix->setNonUnitDiagonal();
    }

    /**
     * @brief N returns the number of rows (which is equal to the number of columns)
     */
    size_t N() const
    {
        return m_matrix->N();
    }

    /**
     * @brief nonzeroes behaves as the Dune::BCRSMatrix::nonzeros() function and returns the number of non zero blocks
     * @return number of non zero blocks.
     */
    size_t nonzeroes() const
    {
        return m_matrix->nonzeroes();
    }

    /**
     * @brief getNonZeroValues returns the GPU vector containing the non-zero values (ordered by block)
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<T>& getNonZeroValues()
    {
        return m_matrix->getNonZeroValues();
    }

    /**
     * @brief getNonZeroValues returns the GPU vector containing the non-zero values (ordered by block)
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<T>& getNonZeroValues() const
    {
        return m_matrix->getNonZeroValues();
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<int>& getRowIndices()
    {
        return m_matrix->getRowIndices();
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<int>& getRowIndices() const
    {
        return m_matrix->getRowIndices();
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    GpuVector<int>& getColumnIndices()
    {
        return m_matrix->getColumnIndices();
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const GpuVector<int>& getColumnIndices() const
    {
        return m_matrix->getColumnIndices();
    }

    /**
     * @brief dim returns the dimension of the vector space on which this matrix acts
     *
     * This is equivalent to matrix.N() * matrix.blockSize()
     * @return matrix.N() * matrix.blockSize()
     */
    size_t dim() const
    {
        return m_matrix->dim();
    }

    /**
     * @brief blockSize size of the blocks
     */
    size_t blockSize() const
    {
        return m_matrix->blockSize();
    }

    /**
     * @brief getDescription the cusparse matrix description.
     *
     * This description is needed for most calls to the CuSparse library
     */
    detail::GpuSparseMatrixDescription& getDescription()
    {
        return m_matrix->getDescription();
    }

    /**
     * @brief mv performs matrix vector multiply y = Ax
     * @param[in] x the vector to multiply the matrix with
     * @param[out] y the output vector
     */
    virtual void mv(const GpuVector<T>& x, GpuVector<T>& y) const
    {
        m_matrix->mv(x, y);
    }

    /**
     * @brief umv computes y=Ax+y
     * @param[in] x the vector to multiply with A
     * @param[inout] y the vector to add and store the output in
     */
    virtual void umv(const GpuVector<T>& x, GpuVector<T>& y) const
    {
        m_matrix->umv(x, y);
    }


    /**
     * @brief umv computes y=alpha * Ax + y
     * @param[in] alpha Scaling factor
     * @param[in] x the vector to multiply with A
     * @param[in,out] y the vector to add and store the output in
     */
    virtual void usmv(T alpha, const GpuVector<T>& x, GpuVector<T>& y) const
    {
        m_matrix->usmv(alpha, x, y);
    }

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
    void updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly = false)
    {
        m_matrix->updateNonzeroValues(matrix, copyNonZeroElementsDirectly);
    }

    /**
     * @brief updateNonzeroValues updates the non-zero values by using the non-zero values of the supplied matrix
     * @param matrix the matrix to extract the non-zero values from
     * @note This assumes the given matrix has the same sparsity pattern.
     */
    void updateNonzeroValues(const GpuSparseMatrixWrapper<T>& matrix)
    {
        m_matrix->updateNonzeroValues(matrix.get());
    }


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
        return m_matrix->dispatchOnBlocksize(function);
    }

private:

    std::unique_ptr<matrix_type> m_matrix = nullptr;
};
} // namespace Opm::gpuistl
#endif
