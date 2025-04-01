/*
  Copyright 2024, SINTEF Digital

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

#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm {

template <typename Scalar, int numWellEq, int numEq>
CompWellEquations<Scalar, numWellEq, numEq>::
CompWellEquations()
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    duneD_.setBuildMode(DiagMatWell::row_wise);
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
init(const int num_conn,  const std::vector<std::size_t>& cells)
{
    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, num_conn, num_conn);
    duneC_.setSize(1, num_conn, num_conn);

    for (auto row = duneD_.createbegin(),
              end = duneD_.createend(); row != end; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }

    for (auto row = duneB_.createbegin(),
                 end = duneB_.createend(); row != end; ++row) {
        for (int con = 0 ; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    // make the C^T matrix
    // TODO: let us see whether we should change the naming of DuneC_
    for (auto row = duneC_.createbegin(),
                 end = duneC_.createend(); row != end; ++row) {
        for (int con = 0; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    resWell_.resize(1);

    Bx_.resize(duneB_.N());

    invDrw_.resize(duneD_.N());

    this->cells_ = cells;
    // some others in the future
}

template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
clear()
{
    duneB_ = 0.0;
    duneC_ = 0.0;
    duneD_ = 0.0;
    resWell_ = 0.0;
}

template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
solve(BVectorWell& dx_well) const
{
    invDuneD_.mv(resWell_, dx_well);
}

template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
invert()
{
    try {
        invDuneD_ = duneD_; // Not strictly need if not cpr with well contributions is used
        detail::invertMatrix(invDuneD_[0][0]);
    } catch (NumericalProblem&) {
        // for singular matrices, use identity as the inverse
        invDuneD_[0][0] = 0.0;
        for (std::size_t i = 0; i < invDuneD_[0][0].rows; ++i) {
           invDuneD_[0][0][i][i] = 1.0;
        }
    }
}


template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
apply(BVector& r) const
{
    assert(invDrw_.size() == invDuneD_.N());

    // invDrw_ = invDuneD_ * resWell_
    invDuneD_.mv(resWell_, invDrw_);
    // r = r - duneC_^T * invDrw_
    duneC_.mmtv(invDrw_, r);
}

template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    BVectorWell resWell = resWell_;
    // resWell = resWell - B * x
    duneB_.mmv(x, resWell);
    // xw = D^-1 * resWell
    invDuneD_.mv(resWell, xw);
}

} // end of namespace Opm
