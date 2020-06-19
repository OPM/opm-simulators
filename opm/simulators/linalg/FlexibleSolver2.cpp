/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.

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

#include "config.h"

#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

// Explicit instantiations of FlexibleSolver

template <int N>
using BV = Dune::BlockVector<Dune::FieldVector<double, N>>;
template <int N>
using BM = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;
template <int N>
using OBM = Dune::BCRSMatrix<Opm::MatrixBlock<double, N, N>>;

// Variants using Dune::FieldMatrix blocks.
template class Dune::FlexibleSolver<BM<2>, BV<2>>;

// Variants using Opm::MatrixBlock blocks.
template class Dune::FlexibleSolver<OBM<2>, BV<2>>;

#if HAVE_MPI

using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;

template Dune::FlexibleSolver<OBM<2>, BV<2>>::FlexibleSolver(const MatrixType& matrix,
                                                             const Comm& comm,
                                                             const boost::property_tree::ptree& prm,
                                                             const std::function<BV<2>()>& weightsCalculator);

#endif // HAVE_MPI
