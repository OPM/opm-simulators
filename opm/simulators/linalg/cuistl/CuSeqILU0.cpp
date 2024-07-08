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
#include <cusparse.h>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/detail/fix_zero_diagonal.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

// This file is based on the guide at https://docs.nvidia.com/cuda/cusparse/index.html#csrilu02_solve ,
// it highly recommended to read that before proceeding.


namespace Opm::cuistl
{

template <class M, class X, class Y, int l>
CuSeqILU0<M, X, Y, l>::CuSeqILU0(const M& A, field_type w)
    : m_underlyingMatrix(A)
    , m_w(w)
    , m_LU(CuSparseMatrix<field_type>::fromMatrix(detail::makeMatrixWithNonzeroDiagonal(A)))
    , m_temporaryStorage(m_LU.N() * m_LU.blockSize())
    , m_descriptionL(detail::createLowerDiagonalDescription())
    , m_descriptionU(detail::createUpperDiagonalDescription())
    , m_cuSparseHandle(detail::CuSparseHandle::getInstance())
{
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_LU.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_LU.N(), A.N()));
    OPM_ERROR_IF(
        A[0][0].N() != m_LU.blockSize(),
        fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.", m_LU.blockSize(), A[0][0].N()));
    OPM_ERROR_IF(
        A.N() * A[0][0].N() != m_LU.dim(),
        fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.", m_LU.dim(), A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != m_LU.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_LU.nonzeroes(),
                             A.nonzeroes()));


    updateILUConfiguration();
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::apply(X& v, const Y& d)
{

    // We need to pass the solve routine a scalar to multiply.
    // In our case this scalar is 1.0
    const field_type one = 1.0;

    const auto numberOfRows = detail::to_int(m_LU.N());
    const auto numberOfNonzeroBlocks = detail::to_int(m_LU.nonzeroes());
    const auto blockSize = detail::to_int(m_LU.blockSize());

    auto nonZeroValues = m_LU.getNonZeroValues().data();
    auto rowIndices = m_LU.getRowIndices().data();
    auto columnIndices = m_LU.getColumnIndices().data();

    // Solve L m_temporaryStorage = d
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_solve(m_cuSparseHandle.get(),
                                                        detail::CUSPARSE_MATRIX_ORDER,
                                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                        numberOfRows,
                                                        numberOfNonzeroBlocks,
                                                        &one,
                                                        m_descriptionL->get(),
                                                        nonZeroValues,
                                                        rowIndices,
                                                        columnIndices,
                                                        blockSize,
                                                        m_infoL.get(),
                                                        d.data(),
                                                        m_temporaryStorage.data(),
                                                        CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                        m_buffer->data()));

    // Solve U v = m_temporaryStorage
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_solve(m_cuSparseHandle.get(),
                                                        detail::CUSPARSE_MATRIX_ORDER,
                                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                        numberOfRows,
                                                        numberOfNonzeroBlocks,
                                                        &one,
                                                        m_descriptionU->get(),
                                                        nonZeroValues,
                                                        rowIndices,
                                                        columnIndices,
                                                        blockSize,
                                                        m_infoU.get(),
                                                        m_temporaryStorage.data(),
                                                        v.data(),
                                                        CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                        m_buffer->data()));


    v *= m_w;
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
CuSeqILU0<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::update()
{
    m_LU.updateNonzeroValues(detail::makeMatrixWithNonzeroDiagonal(m_underlyingMatrix));
    createILU();
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::analyzeMatrix()
{

    if (!m_buffer) {
        OPM_THROW(std::runtime_error,
                  "Buffer not initialized. Call findBufferSize() then initialize with the appropiate size.");
    }
    const auto numberOfRows = detail::to_int(m_LU.N());
    const auto numberOfNonzeroBlocks = detail::to_int(m_LU.nonzeroes());
    const auto blockSize = detail::to_int(m_LU.blockSize());

    auto nonZeroValues = m_LU.getNonZeroValues().data();
    auto rowIndices = m_LU.getRowIndices().data();
    auto columnIndices = m_LU.getColumnIndices().data();
    // analysis of ilu LU decomposition
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrilu02_analysis(m_cuSparseHandle.get(),
                                                             detail::CUSPARSE_MATRIX_ORDER,
                                                             numberOfRows,
                                                             numberOfNonzeroBlocks,
                                                             m_LU.getDescription().get(),
                                                             nonZeroValues,
                                                             rowIndices,
                                                             columnIndices,
                                                             blockSize,
                                                             m_infoM.get(),
                                                             CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                             m_buffer->data()));

    // Make sure we can decompose the matrix.
    int structuralZero;
    auto statusPivot = cusparseXbsrilu02_zeroPivot(m_cuSparseHandle.get(), m_infoM.get(), &structuralZero);
    OPM_ERROR_IF(statusPivot != CUSPARSE_STATUS_SUCCESS,
                 fmt::format("Found a structural zero at A({}, {}). Could not decompose LU \\approx A.\n\n A has "
                             "dimension {}, and has {} nonzeroes.",
                             structuralZero,
                             structuralZero,
                             m_LU.N(),
                             m_LU.nonzeroes()));

    // analysis of ilu apply
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_analysis(m_cuSparseHandle.get(),
                                                           detail::CUSPARSE_MATRIX_ORDER,
                                                           CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                           numberOfRows,
                                                           numberOfNonzeroBlocks,
                                                           m_descriptionL->get(),
                                                           nonZeroValues,
                                                           rowIndices,
                                                           columnIndices,
                                                           blockSize,
                                                           m_infoL.get(),
                                                           CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                           m_buffer->data()));

    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_analysis(m_cuSparseHandle.get(),
                                                           detail::CUSPARSE_MATRIX_ORDER,
                                                           CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                           numberOfRows,
                                                           numberOfNonzeroBlocks,
                                                           m_descriptionU->get(),
                                                           nonZeroValues,
                                                           rowIndices,
                                                           columnIndices,
                                                           blockSize,
                                                           m_infoU.get(),
                                                           CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                           m_buffer->data()));
    m_analysisDone = true;
}

template <class M, class X, class Y, int l>
size_t
CuSeqILU0<M, X, Y, l>::findBufferSize()
{
    // We have three calls that need buffers:
    //   1) LU decomposition
    //   2) solve Lv = y
    //   3) solve Ux = z
    // we combine these buffers into one since it is not used across calls,
    const auto numberOfRows = detail::to_int(m_LU.N());
    const auto numberOfNonzeroBlocks = detail::to_int(m_LU.nonzeroes());
    const auto blockSize = detail::to_int(m_LU.blockSize());

    auto nonZeroValues = m_LU.getNonZeroValues().data();
    auto rowIndices = m_LU.getRowIndices().data();
    auto columnIndices = m_LU.getColumnIndices().data();

    int bufferSizeM = 0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrilu02_bufferSize(m_cuSparseHandle.get(),
                                                               detail::CUSPARSE_MATRIX_ORDER,
                                                               numberOfRows,
                                                               numberOfNonzeroBlocks,
                                                               m_LU.getDescription().get(),
                                                               nonZeroValues,
                                                               rowIndices,
                                                               columnIndices,
                                                               blockSize,
                                                               m_infoM.get(),
                                                               &bufferSizeM));
    int bufferSizeL = 0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_bufferSize(m_cuSparseHandle.get(),
                                                             detail::CUSPARSE_MATRIX_ORDER,
                                                             CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                             numberOfRows,
                                                             numberOfNonzeroBlocks,
                                                             m_descriptionL->get(),
                                                             nonZeroValues,
                                                             rowIndices,
                                                             columnIndices,
                                                             blockSize,
                                                             m_infoL.get(),
                                                             &bufferSizeL));

    int bufferSizeU = 0;
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrsv2_bufferSize(m_cuSparseHandle.get(),
                                                             detail::CUSPARSE_MATRIX_ORDER,
                                                             CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                             numberOfRows,
                                                             numberOfNonzeroBlocks,
                                                             m_descriptionL->get(),
                                                             nonZeroValues,
                                                             rowIndices,
                                                             columnIndices,
                                                             blockSize,
                                                             m_infoU.get(),
                                                             &bufferSizeU));

    OPM_ERROR_IF(bufferSizeL <= 0, fmt::format("bufferSizeL is non-positive. Given value is {}.", bufferSizeL));
    OPM_ERROR_IF(bufferSizeU <= 0, fmt::format("bufferSizeU is non-positive. Given value is {}.", bufferSizeU));
    OPM_ERROR_IF(bufferSizeM <= 0, fmt::format("bufferSizeM is non-positive. Given value is {}.", bufferSizeM));

    return size_t(std::max(bufferSizeL, std::max(bufferSizeU, bufferSizeM)));
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::createILU()
{
    OPM_ERROR_IF(!m_buffer, "Buffer not initialized. Call findBufferSize() then initialize with the appropiate size.");
    OPM_ERROR_IF(!m_analysisDone, "Analyzis of matrix not done. Call analyzeMatrix() first.");

    const auto numberOfRows = detail::to_int(m_LU.N());
    const auto numberOfNonzeroBlocks = detail::to_int(m_LU.nonzeroes());
    const auto blockSize = detail::to_int(m_LU.blockSize());

    auto nonZeroValues = m_LU.getNonZeroValues().data();
    auto rowIndices = m_LU.getRowIndices().data();
    auto columnIndices = m_LU.getColumnIndices().data();
    OPM_CUSPARSE_SAFE_CALL(detail::cusparseBsrilu02(m_cuSparseHandle.get(),
                                                    detail::CUSPARSE_MATRIX_ORDER,
                                                    numberOfRows,
                                                    numberOfNonzeroBlocks,
                                                    m_LU.getDescription().get(),
                                                    nonZeroValues,
                                                    rowIndices,
                                                    columnIndices,
                                                    blockSize,
                                                    m_infoM.get(),
                                                    CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                                    m_buffer->data()));

    // We need to do this here as well. The first call was to check that we could decompose the system A=LU
    // the second call here is to make sure we can solve LUx=y
    int structuralZero;
    // cusparseXbsrilu02_zeroPivot() calls cudaDeviceSynchronize()
    auto statusPivot = cusparseXbsrilu02_zeroPivot(m_cuSparseHandle.get(), m_infoM.get(), &structuralZero);

    OPM_ERROR_IF(
        statusPivot != CUSPARSE_STATUS_SUCCESS,
        fmt::format("Found a structucal zero at LU({}, {}). Could not solve LUx = y.", structuralZero, structuralZero));
}

template <class M, class X, class Y, int l>
void
CuSeqILU0<M, X, Y, l>::updateILUConfiguration()
{
    auto bufferSize = findBufferSize();
    if (!m_buffer || m_buffer->dim() < bufferSize) {
        m_buffer.reset(new CuVector<field_type>((bufferSize + sizeof(field_type) - 1) / sizeof(field_type)));
    }
    analyzeMatrix();
    createILU();
}
} // namespace Opm::cuistl
#define INSTANTIATE_CUSEQILU0_DUNE(realtype, blockdim)                                                                 \
    template class ::Opm::cuistl::CuSeqILU0<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,         \
                                            ::Opm::cuistl::CuVector<realtype>,                                         \
                                            ::Opm::cuistl::CuVector<realtype>>;                                        \
    template class ::Opm::cuistl::CuSeqILU0<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,          \
                                            ::Opm::cuistl::CuVector<realtype>,                                         \
                                            ::Opm::cuistl::CuVector<realtype>>


INSTANTIATE_CUSEQILU0_DUNE(double, 1);
INSTANTIATE_CUSEQILU0_DUNE(double, 2);
INSTANTIATE_CUSEQILU0_DUNE(double, 3);
INSTANTIATE_CUSEQILU0_DUNE(double, 4);
INSTANTIATE_CUSEQILU0_DUNE(double, 5);
INSTANTIATE_CUSEQILU0_DUNE(double, 6);

INSTANTIATE_CUSEQILU0_DUNE(float, 1);
INSTANTIATE_CUSEQILU0_DUNE(float, 2);
INSTANTIATE_CUSEQILU0_DUNE(float, 3);
INSTANTIATE_CUSEQILU0_DUNE(float, 4);
INSTANTIATE_CUSEQILU0_DUNE(float, 5);
INSTANTIATE_CUSEQILU0_DUNE(float, 6);
