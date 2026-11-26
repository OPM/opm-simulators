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
 * \brief Selects the BCRS matrix implementation used by OPM.
 */
#ifndef OPM_BCRSMATRIX_HEADER_INCLUDED
#define OPM_BCRSMATRIX_HEADER_INCLUDED

#if defined(OPM_USE_EXPERIMENTAL_BCRSMATRIX) && OPM_USE_EXPERIMENTAL_BCRSMATRIX
#include <dune/istl/experimental/bcrsmatrix.hh>
#else
#include <dune/istl/bcrsmatrix.hh>
#endif

namespace Opm {

#if defined(OPM_USE_EXPERIMENTAL_BCRSMATRIX) && OPM_USE_EXPERIMENTAL_BCRSMATRIX
/// \brief Block-compressed row-storage matrix used by the ISTL solvers.
///
/// The experimental implementation is selected when
/// `OPM_USE_EXPERIMENTAL_BCRSMATRIX` is enabled.
/// \tparam B Scalar or block type stored in each nonzero matrix entry.
/// \tparam A Allocator used for matrix blocks.
template<class B, class A = std::allocator<B>>
using BCRSMatrix = Dune::Experimental::BCRSMatrix<B, std::uint_least32_t, A>;
#else
/// \brief Block-compressed row-storage matrix used by the ISTL solvers.
///
/// This aliases DUNE's standard implementation when the experimental
/// implementation is not enabled.
/// \tparam B Scalar or block type stored in each nonzero matrix entry.
/// \tparam A Allocator used for matrix blocks.
template<class B, class A = std::allocator<B>>
using BCRSMatrix = Dune::BCRSMatrix<B, A>;
#endif

} // namespace Opm

#endif
