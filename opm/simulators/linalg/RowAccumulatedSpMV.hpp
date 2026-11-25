/*
  Copyright 2026 SINTEF Digital

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

#ifndef OPM_ROW_ACCUMULATED_SPMV_HEADER_INCLUDED
#define OPM_ROW_ACCUMULATED_SPMV_HEADER_INCLUDED

#include <dune/istl/bcrsmatrix.hh>

#include <cstddef>
#include <type_traits>
#include <utility>

// Stand-in for BCRSMatrix::mv/usmv until dune-istl accumulates rows locally itself;
// then the operators can call A.mv/A.usmv again and this header can go.
namespace Opm::detail {

// BCRSMatrix or a class derived from it (e.g. OverlappingBCRSMatrix). Other matrix
// types (e.g. MixedMatrixWrapper) only provide mv/usmv, not row iterators.
template <class B, class A>
std::true_type isBCRSMatrix(const Dune::BCRSMatrix<B, A>*);
std::false_type isBCRSMatrix(const void*);

template <class M>
using IsBCRSMatrix = decltype(isBCRSMatrix(std::declval<const M*>()));

// y = A x and y += alpha A x over rows [0, rowEnd), summing each row in a local block:
// accumulating straight into y[i] makes the compiler store y[i] after every block.
// Same per-block operations in the same order as BCRSMatrix::mv/usmv.
template <class M, class X, class Y>
void rowAccumulatedMv(const M& A, const X& x, Y& y, std::size_t rowEnd)
{
    for (auto row = A.begin(); row != A.end() && row.index() < rowEnd; ++row) {
        auto acc = y[row.index()];
        acc = 0;
        for (auto col = row->begin(); col != row->end(); ++col)
            col->umv(x[col.index()], acc);
        y[row.index()] = acc;
    }
}

template <class M, class X, class Y, class F>
void rowAccumulatedUsmv(const F& alpha, const M& A, const X& x, Y& y, std::size_t rowEnd)
{
    for (auto row = A.begin(); row != A.end() && row.index() < rowEnd; ++row) {
        auto acc = y[row.index()];
        for (auto col = row->begin(); col != row->end(); ++col)
            col->usmv(alpha, x[col.index()], acc);
        y[row.index()] = acc;
    }
}

// Whole-matrix versions; fall back to the matrix's own mv/usmv for non-BCRS types.
template <class M, class X, class Y>
void mv(const M& A, const X& x, Y& y)
{
    if constexpr (IsBCRSMatrix<M>::value)
        rowAccumulatedMv(A, x, y, A.N());
    else
        A.mv(x, y);
}

template <class M, class X, class Y, class F>
void usmv(const F& alpha, const M& A, const X& x, Y& y)
{
    if constexpr (IsBCRSMatrix<M>::value)
        rowAccumulatedUsmv(alpha, A, x, y, A.N());
    else
        A.usmv(alpha, x, y);
}

} // namespace Opm::detail

#endif // OPM_ROW_ACCUMULATED_SPMV_HEADER_INCLUDED
