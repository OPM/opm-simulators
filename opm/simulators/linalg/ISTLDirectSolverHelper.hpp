/*
  Copyright 2026 OPM-OP

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

#ifndef OPM_ISTL_DIRECT_SOLVER_HELPER_HEADER_INCLUDED
#define OPM_ISTL_DIRECT_SOLVER_HELPER_HEADER_INCLUDED

#include <dune-istl-config.hh>
#include <dune/common/version.hh>

//! \file Utilites to handle differences in DUNE ISTL versions for direct solvers.
//! \details Newer versions of DUNE ISTL (2.12 and later) have changed the number of parameters of UMFPack and SuperLU solvers.
#if DUNE_VERSION_GT(DUNE_ISTL, 2, 11)
namespace Dune {
    template<typename M, typename X, typename Y>
    class UMFPack;

    template<typename M, typename X, typename Y>
    class SuperLU;
}

namespace Opm {
    template<typename M, typename X, typename Y>
    using UMFPack = Dune::UMFPack<M, X, Y>;

    template<typename M, typename X, typename Y>
    using SuperLU = Dune::SuperLU<M, X, Y>;
}

#else
namespace Dune {
    template<typename M>
    class UMFPack;

    template<typename M>
    class SuperLU;
}

namespace Opm {
    template<typename M, typename X, typename Y>
    using UMFPack = Dune::UMFPack<M>;

    template<typename M, typename X, typename Y>
    using SuperLU = Dune::SuperLU<M>;
}
#endif

#endif // OPM_ISTL_DIRECT_SOLVER_HELPER_HEADER_INCLUDED
