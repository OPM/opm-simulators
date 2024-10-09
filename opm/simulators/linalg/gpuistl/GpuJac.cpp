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
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/GpuJac.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/JacKernels.hpp>
#include <opm/simulators/linalg/gpuistl/detail/vector_operations.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm::gpuistl
{

template <class M, class X, class Y, int l>
GpuJac<M, X, Y, l>::GpuJac(const M& A, field_type w)
    : m_cpuMatrix(A)
    , m_relaxationFactor(w)
    , m_gpuMatrix(GpuSparseMatrix<field_type>::fromMatrix(A))
    , m_diagInvFlattened(m_gpuMatrix.N() * m_gpuMatrix.blockSize() * m_gpuMatrix.blockSize())
{
    // Some sanity check
    OPM_ERROR_IF(A.N() != m_gpuMatrix.N(),
                 fmt::format("CuSparse matrix not same size as DUNE matrix. {} vs {}.", m_gpuMatrix.N(), A.N()));
    OPM_ERROR_IF(A[0][0].N() != m_gpuMatrix.blockSize(),
                 fmt::format("CuSparse matrix not same blocksize as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.blockSize(),
                             A[0][0].N()));
    OPM_ERROR_IF(A.N() * A[0][0].N() != m_gpuMatrix.dim(),
                 fmt::format("CuSparse matrix not same dimension as DUNE matrix. {} vs {}.",
                             m_gpuMatrix.dim(),
                             A.N() * A[0][0].N()));
    OPM_ERROR_IF(A.nonzeroes() != m_gpuMatrix.nonzeroes(),
                 fmt::format("CuSparse matrix not same number of non zeroes as DUNE matrix. {} vs {}. ",
                             m_gpuMatrix.nonzeroes(),
                             A.nonzeroes()));

    // Compute the inverted diagonal of A and store it in a vector format in m_diagInvFlattened
    invertDiagonalAndFlatten();
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::pre([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
{
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::apply(X& v, const Y& d)
{
    // Jacobi preconditioner: x_{n+1} = x_n + w * (D^-1 * (b - Ax_n) )
    // Working with defect d and update v it we only need to set v = w*(D^-1)*d

    // Compute the MV product where the matrix is diagonal and therefore stored as a vector.
    // The product is thus computed as a hadamard product.
    detail::weightedDiagMV(
        m_diagInvFlattened.data(), m_gpuMatrix.N(), m_gpuMatrix.blockSize(), m_relaxationFactor, d.data(), v.data());
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::post([[maybe_unused]] X& x)
{
}

template <class M, class X, class Y, int l>
Dune::SolverCategory::Category
GpuJac<M, X, Y, l>::category() const
{
    return Dune::SolverCategory::sequential;
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::update()
{
    m_gpuMatrix.updateNonzeroValues(m_cpuMatrix);
    invertDiagonalAndFlatten();
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::invertDiagonalAndFlatten()
{
    detail::JAC::invertDiagonalAndFlatten<field_type, matrix_type::block_type::cols>(
        m_gpuMatrix.getNonZeroValues().data(),
        m_gpuMatrix.getRowIndices().data(),
        m_gpuMatrix.getColumnIndices().data(),
        m_gpuMatrix.N(),
        m_diagInvFlattened.data());
}

} // namespace Opm::gpuistl
#define INSTANTIATE_CUJAC_DUNE(realtype, blockdim)                                                                     \
    template class ::Opm::gpuistl::GpuJac<Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,             \
                                        ::Opm::gpuistl::GpuVector<realtype>,                                             \
                                        ::Opm::gpuistl::GpuVector<realtype>>;                                            \
    template class ::Opm::gpuistl::GpuJac<Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,              \
                                        ::Opm::gpuistl::GpuVector<realtype>,                                             \
                                        ::Opm::gpuistl::GpuVector<realtype>>

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
