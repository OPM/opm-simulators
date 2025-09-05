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
#include <stdexcept>

namespace Opm::gpuistl
{

template <class M, class X, class Y, int l>
GpuJac<M, X, Y, l>::GpuJac(const M& A, field_type w)
    : m_matrix(A)
    , m_relaxationFactor(w)
    , m_diagInvFlattened(m_matrix.N() * m_matrix.blockSize() * m_matrix.blockSize())
{

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
        m_diagInvFlattened.data(), m_matrix.N(), m_matrix.blockSize(), m_relaxationFactor, d.data(), v.data());
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
    invertDiagonalAndFlatten();
}

template <class M, class X, class Y, int l>
template<int blocksize>
void
GpuJac<M, X, Y, l>::dispatchInvertDiagonalAndFlatten()
{
    if (m_matrix.blockSize() != blocksize) {
        if constexpr (blocksize > 1) {
            dispatchInvertDiagonalAndFlatten<blocksize - 1>();
        } else {
            OPM_THROW(std::runtime_error, fmt::format("Block size {} is not supported.", m_matrix.blockSize()));
        }
    } else {
        detail::JAC::invertDiagonalAndFlatten<field_type, blocksize>(
            m_matrix.getNonZeroValues().data(),
            m_matrix.getRowIndices().data(),
            m_matrix.getColumnIndices().data(),
            m_matrix.N(),
            m_diagInvFlattened.data());
    }
}

template <class M, class X, class Y, int l>
void
GpuJac<M, X, Y, l>::invertDiagonalAndFlatten()
{
    dispatchInvertDiagonalAndFlatten<6>();   
}

} // namespace Opm::gpuistl
#define INSTANTIATE_CUJAC_DUNE(realtype)                                                         \
    template class ::Opm::gpuistl::GpuJac<::Opm::gpuistl::GpuSparseMatrixWrapper<realtype>,             \
                                        ::Opm::gpuistl::GpuVector<realtype>,                     \
                                        ::Opm::gpuistl::GpuVector<realtype>>
   
INSTANTIATE_CUJAC_DUNE(double);
INSTANTIATE_CUJAC_DUNE(float);
