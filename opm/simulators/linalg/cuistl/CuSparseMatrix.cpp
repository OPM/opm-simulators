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
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <fmt/core.h>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_wrapper.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm::cuistl
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
        // TODO: [perf] Can we avoid building nonZeroElementsData?
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
CuSparseMatrix<T>::CuSparseMatrix(const T* nonZeroElements,
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
}

template <class T>
CuSparseMatrix<T>::~CuSparseMatrix()
{
    // empty
}

template <typename T>
template <typename MatrixType>
CuSparseMatrix<T>
CuSparseMatrix<T>::fromMatrix(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
{
    // TODO: Do we need this intermediate storage? Or this shuffling of data?
    std::vector<int> columnIndices;
    std::vector<int> rowIndices;

    rowIndices.push_back(0);

    const size_t blockSize = matrix[0][0].N();
    const size_t numberOfRows = matrix.N();
    const size_t numberOfNonzeroBlocks = matrix.nonzeroes();

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
                 "Error size of rows do not sum to number of nonzeroes in CuSparseMatrix.");
    OPM_ERROR_IF(rowIndices.size() != numberOfRows + 1, "Row indices do not match for CuSparseMatrix.");
    OPM_ERROR_IF(columnIndices.size() != numberOfNonzeroBlocks, "Column indices do not match for CuSparseMatrix.");


    if (copyNonZeroElementsDirectly) {
        const T* nonZeroElements = static_cast<const T*>(&((matrix[0][0][0][0])));
        return CuSparseMatrix<T>(
            nonZeroElements, rowIndices.data(), columnIndices.data(), numberOfNonzeroBlocks, blockSize, numberOfRows);
    } else {
        auto nonZeroElementData = extractNonzeroValues<T>(matrix);
        return CuSparseMatrix<T>(nonZeroElementData.data(),
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
CuSparseMatrix<T>::updateNonzeroValues(const MatrixType& matrix, bool copyNonZeroElementsDirectly)
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
}

template <typename T>
void
CuSparseMatrix<T>::setUpperTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(m_matrixDescription->get(), CUSPARSE_FILL_MODE_UPPER));
}

template <typename T>
void
CuSparseMatrix<T>::setLowerTriangular()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(m_matrixDescription->get(), CUSPARSE_FILL_MODE_LOWER));
}

template <typename T>
void
CuSparseMatrix<T>::setUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(m_matrixDescription->get(), CUSPARSE_DIAG_TYPE_UNIT));
}

template <typename T>
void
CuSparseMatrix<T>::setNonUnitDiagonal()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(m_matrixDescription->get(), CUSPARSE_DIAG_TYPE_NON_UNIT));
}

template <typename T>
void
CuSparseMatrix<T>::mv(const CuVector<T>& x, CuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);
    if (blockSize() < 2u) {
        OPM_THROW(
            std::invalid_argument,
            "CuSparseMatrix<T>::usmv and CuSparseMatrix<T>::mv are only implemented for block sizes greater than 1.");
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
CuSparseMatrix<T>::umv(const CuVector<T>& x, CuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);
    if (blockSize() < 2u) {
        OPM_THROW(
            std::invalid_argument,
            "CuSparseMatrix<T>::usmv and CuSparseMatrix<T>::mv are only implemented for block sizes greater than 1.");
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
CuSparseMatrix<T>::usmv(T alpha, const CuVector<T>& x, CuVector<T>& y) const
{
    assertSameSize(x);
    assertSameSize(y);
    if (blockSize() < 2) {
        OPM_THROW(
            std::invalid_argument,
            "CuSparseMatrix<T>::usmv and CuSparseMatrix<T>::mv are only implemented for block sizes greater than 1.");
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
CuSparseMatrix<T>::assertSameSize(const VectorType& x) const
{
    if (x.dim() != blockSize() * N()) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Size mismatch. Input vector has {} elements, while we have {} rows.",
                              x.dim(),
                              blockSize() * N()));
    }
}



#define INSTANTIATE_CUSPARSE_DUNE_MATRIX_CONSTRUCTION_FUNTIONS(realtype, blockdim)                                     \
    template CuSparseMatrix<realtype> CuSparseMatrix<realtype>::fromMatrix(                                            \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template CuSparseMatrix<realtype> CuSparseMatrix<realtype>::fromMatrix(                                            \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool);                                \
    template void CuSparseMatrix<realtype>::updateNonzeroValues(                                                       \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template void CuSparseMatrix<realtype>::updateNonzeroValues(                                                       \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool)

template class CuSparseMatrix<float>;
template class CuSparseMatrix<double>;

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

} // namespace Opm::cuistl
