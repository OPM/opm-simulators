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
#ifndef OPM_CUSPARSEMATRIX_HPP
#define OPM_CUSPARSEMATRIX_HPP
#include <cusparse.h>
#include <iostream>
#include <memory>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <vector>

namespace Opm::cuistl
{

/**
 * @brief The CuSparseMatrix class simple wrapper class for a CuSparse matrix.
 *
 * @note we currently only support simple raw primitives for T (double and float). Block size is handled through the
 * block size parameter
 *
 * @tparam T the type to store. Can be either float, double or int.
 *
 * @note we only support square matrices.
 *
 * @note We only support Block Compressed Sparse Row Format (BSR) for now.
 */
template <typename T>
class CuSparseMatrix
{
public:
    //! Create the sparse matrix specified by the raw data.
    //!
    //! \note Prefer to use the constructor taking a const reference to a matrix instead.
    //!
    //! \param[in] nonZeroElements the non-zero values of the matrix
    //! \param[in] rowIndices      the row indices of the non-zero elements
    //! \param[in] columnIndices   the column indices of the non-zero elements
    //! \param[in] numberOfNonzeroElements number of nonzero elements
    //! \param[in] blockSize size of each block matrix (typically 3)
    //! \param[in] numberOfRows the number of rows
    //!
    //! \note We assume numberOfNonzeroBlocks, blockSize and numberOfRows all are representable as int due to
    //!       restrictions in the current version of cusparse. This might change in future versions.
    CuSparseMatrix(const T* nonZeroElements,
                   const int* rowIndices,
                   const int* columnIndices,
                   size_t numberOfNonzeroBlocks,
                   size_t blockSize,
                   size_t numberOfRows);

    /**
     * We don't want to be able to copy this for now (too much hassle in copying the cusparse resources)
     */
    CuSparseMatrix(const CuSparseMatrix&) = delete;

    /**
     * We don't want to be able to copy this for now (too much hassle in copying the cusparse resources)
     */
    CuSparseMatrix& operator=(const CuSparseMatrix&) = delete;

    virtual ~CuSparseMatrix();

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
    static CuSparseMatrix<T> fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly = false);

    /**
     * @brief setUpperTriangular sets the CuSparse flag that this is an upper diagonal (with unit diagonal) matrix.
     */
    void setUpperTriangular();

    /**
     * @brief setLowerTriangular sets the CuSparse flag that this is an lower diagonal (with non-unit diagonal) matrix.
     */
    void setLowerTriangular();

    /**
     * @brief setUnitDiagonal sets the CuSparse flag that this has unit diagional.
     */
    void setUnitDiagonal();


    /**
     * @brief setNonUnitDiagonal sets the CuSparse flag that this has non-unit diagional.
     */
    void setNonUnitDiagonal();

    /**
     * @brief N returns the number of rows (which is equal to the number of columns)
     */
    size_t N() const
    {
        // Technically this safe conversion is not needed since we enforce these to be
        // non-negative in the constructor, but keeping them for added sanity for now.
        //
        // We don't believe this will yield any performance penality (it's used too far away from the inner loop),
        // but should that be false, they can be removed.
        return detail::to_size_t(m_numberOfRows);
    }

    /**
     * @brief nonzeroes behaves as the Dune::BCRSMatrix::nonzeros() function and returns the number of non zero blocks
     * @return number of non zero blocks.
     */
    size_t nonzeroes() const
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
    CuVector<T>& getNonZeroValues()
    {
        return m_nonZeroElements;
    }

    /**
     * @brief getNonZeroValues returns the GPU vector containing the non-zero values (ordered by block)
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const CuVector<T>& getNonZeroValues() const
    {
        return m_nonZeroElements;
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    CuVector<int>& getRowIndices()
    {
        return m_rowIndices;
    }

    /**
     * @brief getRowIndices returns the row indices used to represent the BSR structure.
     *
     * @note Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const CuVector<int>& getRowIndices() const
    {
        return m_rowIndices;
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    CuVector<int>& getColumnIndices()
    {
        return m_columnIndices;
    }

    /**
     * @brief getColumnIndices returns the column indices used to represent the BSR structure.
     *
     * @return Read the CuSPARSE documentation on Block Compressed Sparse Row Format (BSR) for the exact ordering.
     */
    const CuVector<int>& getColumnIndices() const
    {
        return m_columnIndices;
    }

    /**
     * @brief dim returns the dimension of the vector space on which this matrix acts
     *
     * This is equivalent to matrix.N() * matrix.blockSize()
     * @return matrix.N() * matrix.blockSize()
     */
    size_t dim() const
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
    size_t blockSize() const
    {
        // Technically this safe conversion is not needed since we enforce these to be
        // non-negative in the constructor, but keeping them for added sanity for now.
        //
        // We don't believe this will yield any performance penality (it's used too far away from the inner loop),
        // but should that be false, they can be removed.
        return detail::to_size_t(m_blockSize);
    }

    /**
     * @brief getDescription the cusparse matrix description.
     *
     * This description is needed for most calls to the CuSparse library
     */
    detail::CuSparseMatrixDescription& getDescription()
    {
        return *m_matrixDescription;
    }

    /**
     * @brief mv performs matrix vector multiply y = Ax
     * @param[in] x the vector to multiply the matrix with
     * @param[out] y the output vector
     *
     * @note Due to limitations of CuSparse, this is only supported for block sizes greater than 1.
     */
    virtual void mv(const CuVector<T>& x, CuVector<T>& y) const;

    /**
     * @brief umv computes y=Ax+y
     * @param[in] x the vector to multiply with A
     * @param[inout] y the vector to add and store the output in
     *
     * @note Due to limitations of CuSparse, this is only supported for block sizes greater than 1.
     */
    virtual void umv(const CuVector<T>& x, CuVector<T>& y) const;


    /**
     * @brief umv computes y=alpha * Ax + y
     * @param[in] x the vector to multiply with A
     * @param[inout] y the vector to add and store the output in
     *
     * @note Due to limitations of CuSparse, this is only supported for block sizes greater than 1.
     */
    virtual void usmv(T alpha, const CuVector<T>& x, CuVector<T>& y) const;

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

private:
    CuVector<T> m_nonZeroElements;
    CuVector<int> m_columnIndices;
    CuVector<int> m_rowIndices;

    // Notice that we store these three as int to make sure we are cusparse compatible.
    //
    // This gives the added benefit of checking the size constraints at construction of the matrix
    // rather than in some call to cusparse.
    const int m_numberOfNonzeroBlocks;
    const int m_numberOfRows;
    const int m_blockSize;

    detail::CuSparseMatrixDescriptionPtr m_matrixDescription;
    detail::CuSparseHandle& m_cusparseHandle;

    template <class VectorType>
    void assertSameSize(const VectorType& vector) const;
};
} // namespace Opm::cuistl
#endif
