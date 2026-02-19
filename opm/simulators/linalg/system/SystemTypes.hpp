#pragma once

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm
{

inline constexpr int numResDofs = 3;
inline constexpr int numWellDofs = 4;

using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;

using ResVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;

using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                                                Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
using SystemVector = Dune::MultiTypeBlockVector<ResVector, WellVector>;

} // namespace Opm
