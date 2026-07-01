// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 */
#ifndef OPM_BLOCKSPARSEMATRIX_HEADER_INCLUDED
#define OPM_BLOCKSPARSEMATRIX_HEADER_INCLUDED

#if defined(OPM_USE_EXPERIMENTAL_IBCRSMATRIX) && OPM_USE_EXPERIMENTAL_IBCRSMATRIX
#include <dune/istl/ibcrsmatrix.hh>
#else
#include <dune/istl/bcrsmatrix.hh>
#endif

namespace Opm {

#if defined(OPM_USE_EXPERIMENTAL_IBCRSMATRIX) && OPM_USE_EXPERIMENTAL_IBCRSMATRIX
// Use the exprimental IBCRSMatrix, which has a better performance under distributed workloads.
template<class B, class A = std::allocator<B>>
using BlockSparseMatrix = Dune::IBCRSMatrix<B, uint_least32_t, A>;
#else
// Use the classic BCRSMatrix.
template<class B, class A = std::allocator<B>>
using BlockSparseMatrix = Dune::BCRSMatrix<B, A>;
#endif

} // namespace Opm

#endif
