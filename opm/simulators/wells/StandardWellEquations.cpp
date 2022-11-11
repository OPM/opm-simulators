/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <config.h>
#include <opm/simulators/wells/StandardWellEquations.hpp>

#include <opm/common/Exceptions.hpp>

#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm
{

template<class Scalar, int numEq>
StandardWellEquations<Scalar,numEq>::
StandardWellEquations(const ParallelWellInfo& parallel_well_info)
    : parallelB_(duneB_, parallel_well_info)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::
init(const int num_cells,
     const int numWellEq,
     const int numPerfs,
     const std::vector<int>& cells)
{
    // setup sparsity pattern for the matrices
    //[A C^T    [x    =  [ res
    // B D] x_well]      res_well]
    // set the size of the matrices
    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, num_cells, numPerfs);
    duneC_.setSize(1, num_cells, numPerfs);

    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }
      // the block size is run-time determined now
    duneD_[0][0].resize(numWellEq, numWellEq);

    for (auto row = duneB_.createbegin(),
              end = duneB_.createend(); row != end; ++row) {
        for (int perf = 0 ; perf < numPerfs; ++perf) {
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0 ; perf < numPerfs; ++perf) {
        const int cell_idx = cells[perf];
        // the block size is run-time determined now
        duneB_[0][cell_idx].resize(numWellEq, numEq);
    }

         // make the C^T matrix
    for (auto row = duneC_.createbegin(),
              end = duneC_.createend(); row != end; ++row) {
        for (int perf = 0; perf < numPerfs; ++perf) {
            const int cell_idx = cells[perf];
            row.insert(cell_idx);
        }
    }

    for (int perf = 0; perf < numPerfs; ++perf) {
        const int cell_idx = cells[perf];
        duneC_[0][cell_idx].resize(numWellEq, numEq);
    }

    resWell_.resize(1);
    // the block size of resWell_ is also run-time determined now
    resWell_[0].resize(numWellEq);

    // resize temporary class variables
    Bx_.resize(duneB_.N());
    for (unsigned i = 0; i < duneB_.N(); ++i) {
        Bx_[i].resize(numWellEq);
    }

    invDrw_.resize(duneD_.N());
    for (unsigned i = 0; i < duneD_.N(); ++i) {
        invDrw_[i].resize(numWellEq);
    }
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::apply(const BVector& x, BVector& Ax) const
{
    assert(Bx_.size() == duneB_.N());
    assert(invDrw_.size() == invDuneD_.N());

    // Bx_ = duneB_ * x
    parallelB_.mv(x, Bx_);

    // invDBx = invDuneD_ * Bx_
    // TODO: with this, we modified the content of the invDrw_.
    // Is it necessary to do this to save some memory?
    auto& invDBx = invDrw_;
    invDuneD_.mv(Bx_, invDBx);

    // Ax = Ax - duneC_^T * invDBx
    duneC_.mmtv(invDBx, Ax);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::apply(BVector& r) const
{
    assert(invDrw_.size() == invDuneD_.N());

    // invDrw_ = invDuneD_ * resWell_
    invDuneD_.mv(resWell_, invDrw_);
    // r = r - duneC_^T * invDrw_
    duneC_.mmtv(invDrw_, r);
}

template<class Scalar, int numEq>
void StandardWellEquations<Scalar,numEq>::invert()
{
    try {
        invDuneD_ = duneD_; // Not strictly need if not cpr with well contributions is used
        detail::invertMatrix(invDuneD_[0][0]);
    } catch (NumericalProblem&) {
        // for singular matrices, use identity as the inverse
        invDuneD_[0][0] = 0.0;
        for (size_t i = 0; i < invDuneD_[0][0].rows(); ++i) {
            invDuneD_[0][0][i][i] = 1.0;
        }
    }
}

#define INSTANCE(N) \
template class StandardWellEquations<double,N>;

INSTANCE(1)
INSTANCE(2)
INSTANCE(3)
INSTANCE(4)
INSTANCE(5)
INSTANCE(6)

}
