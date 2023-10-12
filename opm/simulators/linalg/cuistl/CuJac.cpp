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
#include <opm/simulators/linalg/cuistl/CuJac.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuBlasHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_constants.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/detail/fix_zero_diagonal.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm::cuistl
{

template <class M, class X, class Y, int l>
CuJac<M, X, Y, l>::CuJac(const M& A, field_type w)
    : cpuMatrix(A)
    , relaxation_factor(w)
    , cuMatrix(CuSparseMatrix<field_type>::fromMatrix(detail::makeMatrixWithNonzeroDiagonal(A)))
    , cuVec_diagInvFlattened(cuMatrix.N() * cuMatrix.blockSize() * cuMatrix.blockSize())
{
    // Some sanity check
    OPM_ERROR_IF(A.N() != cuMatrix.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", cuMatrix.N(), A.N()));
    OPM_ERROR_IF(
        A[0][0].N() != cuMatrix.blockSize(),
        fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.", cuMatrix.blockSize(), A[0][0].N()));
    OPM_ERROR_IF(A.N() * A[0][0].N() != cuMatrix.dim(),
                 fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.",
                             cuMatrix.dim(),
                             A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != cuMatrix.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             cuMatrix.nonzeroes(),
                             A.nonzeroes()));

    // Compute the inverted diagonal of A and store it in a vector format in cuVec_diagInvFlattened
    detail::invertDiagonalAndFlatten(cuMatrix.getNonZeroValues().data(),
                                     cuMatrix.getRowIndices().data(),
                                     cuMatrix.getColumnIndices().data(),
                                     detail::to_int(cuMatrix.N()),
                                     detail::to_int(cuMatrix.blockSize()),
                                     cuVec_diagInvFlattened.data());
}

template <class M, class X, class Y, int l>
void
CuJac<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
CuJac<M, X, Y, l>::apply(X& v, const Y& d)
{   
    // Jacobi preconditioner: x_{n+1} = x_n + w * (D^-1 * (b - Ax_n) )
    // Working with defect d and update v it we only need to set v = w*(D^-1)*d

    const auto numberOfRows = detail::to_int(cuMatrix.N());
    const auto numberOfNonzeroBlocks = detail::to_int(cuMatrix.nonzeroes());
    const auto blockSize = detail::to_int(cuMatrix.blockSize());

    // Could probably be further optimized a little bit by fusing the copy with the 
    // following multiplication. To do so add a kernel that does not do the multiplication
    // in place.

    v = d; // cuda copy
    detail::blockVectorMultiplicationAtAllIndices(cuVec_diagInvFlattened.data(),
                                                  detail::to_size_t(numberOfRows),
                                                  detail::to_size_t(blockSize),
                                                  v.data());
    v *= relaxation_factor; // cuBlas axpy
}

template <class M, class X, class Y, int l>
void
CuJac<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
CuJac<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
CuJac<M, X, Y, l>::update()
{
    cuMatrix.updateNonzeroValues(detail::makeMatrixWithNonzeroDiagonal(cpuMatrix));
    detail::invertDiagonalAndFlatten(cuMatrix.getNonZeroValues().data(),
                                     cuMatrix.getRowIndices().data(),
                                     cuMatrix.getColumnIndices().data(),
                                     detail::to_int(cuMatrix.N()),
                                     detail::to_int(cuMatrix.blockSize()),
                                     cuVec_diagInvFlattened.data());
}

} // namespace Opm::cuistl
#define INSTANTIATE_CUJAC_DUNE(realtype, blockdim)                                                                     \
    template class ::Opm::cuistl::CuJac<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,             \
                                        ::Opm::cuistl::CuVector<realtype>,                                             \
                                        ::Opm::cuistl::CuVector<realtype>>;                                            \
    template class ::Opm::cuistl::CuJac<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,              \
                                        ::Opm::cuistl::CuVector<realtype>,                                             \
                                        ::Opm::cuistl::CuVector<realtype>>

INSTANTIATE_CUJAC_DUNE(double, 1);
INSTANTIATE_CUJAC_DUNE(double, 2);
INSTANTIATE_CUJAC_DUNE(double, 3);
INSTANTIATE_CUJAC_DUNE(double, 4);
INSTANTIATE_CUJAC_DUNE(double, 5);
INSTANTIATE_CUJAC_DUNE(double, 6);

INSTANTIATE_CUJAC_DUNE(float, 1);
INSTANTIATE_CUJAC_DUNE(float, 2);
INSTANTIATE_CUJAC_DUNE(float, 3);
INSTANTIATE_CUJAC_DUNE(float, 4);
INSTANTIATE_CUJAC_DUNE(float, 5);
INSTANTIATE_CUJAC_DUNE(float, 6);
