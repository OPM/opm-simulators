/*
  Copyright 2024, 2026, SINTEF Digital

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
#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>

#include <cmath>

namespace Opm {

template <typename Scalar, int numWellEq, int numEq>
CompWellEquations<Scalar, numWellEq, numEq>::
CompWellEquations()
{
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template <typename Scalar, int numWellEq, int numEq>
void
CompWellEquations<Scalar, numWellEq, numEq>::
init(const int num_conn,  const std::vector<std::size_t>& cells)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise),
    duneD_.setBuildMode(DiagMatWell::row_wise);

    duneD_.setSize(1, 1, 1);
    duneB_.setSize(1, num_conn, num_conn);
    duneC_.setSize(1, num_conn, num_conn);

    auto endD = duneD_.createend();
    for (auto row = duneD_.createbegin(); row != endD; ++row) {
        // Add nonzeros for diagonal
        row.insert(row.index());
    }

    auto endB = duneB_.createend();
    for (auto row = duneB_.createbegin(); row != endB; ++row) {
        for (int con = 0 ; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    // make the C^T matrix
    // TODO: let us see whether we should change the naming of DuneC_
    auto endC = duneC_.createend();
    for (auto row = duneC_.createbegin(); row != endC; ++row) {
        for (int con = 0; con < num_conn; ++con) {
            row.insert(con);
        }
    }

    resWell_.resize(1);

    Bx_.resize(duneB_.N());

    invDrw_.resize(duneD_.N());

    this->cells_ = cells;
    this->res_scales_.assign(cells.size(), 1.0);
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
sumAndPinRows(const int last_row, const int variable_offset)
{
    auto& D = duneD_[0][0];
    auto& residual = resWell_[0];
    for (int row = 0; row < last_row; ++row) {
        D[last_row] += D[row];
        residual[last_row] += residual[row];
        D[row] = 0.0;
        D[row][row + variable_offset] = 1.0;
        residual[row] = 0.0;
    }
    for (auto& B : duneB_[0]) {
        for (int row = 0; row < last_row; ++row) {
            B[last_row] += B[row];
            B[row] = 0.0;
        }
    }
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
    // A singular well matrix has no valid Schur complement. Let the timestep
    // recovery handle it rather than treating the identity as its inverse.
    bool singular = false;
    try {
        invDuneD_ = duneD_; // Not strictly need if not cpr with well contributions is used
        detail::invertMatrix(invDuneD_[0][0]);
    } catch (const NumericalProblem&) {
        singular = true;
    } catch (const Dune::FMatrixError&) {
        singular = true;
    }

    // Catch the silent inf/NaN paths that do not throw.
    for (std::size_t i = 0; !singular && i < invDuneD_[0][0].rows; ++i) {
        for (std::size_t j = 0; !singular && j < invDuneD_[0][0].cols; ++j) {
            if (!std::isfinite(invDuneD_[0][0][i][j])) {
                singular = true;
            }
        }
    }

    if (singular) {
        throw NumericalProblem("Singular compositional well matrix");
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
    // r = r - scale * duneC_^T * invDrw_, with the per-cell scale converting
    // the well-equation units to the reservoir residual's units. r holds the
    // well's own cells, one entry per connection.
    for (auto colC = duneC_[0].begin(), endC = duneC_[0].end(); colC != endC; ++colC) {
        VectorBlockType tmp(0.0);
        (*colC).usmtv(res_scales_[colC.index()], invDrw_[0], tmp);
        r[colC.index()] -= tmp;
    }
}

template <typename Scalar, int numWellEq, int numEq>
template <class SparseMatrixAdapter>
void
CompWellEquations<Scalar, numWellEq, numEq>::
extract(SparseMatrixAdapter& jacobian) const
{
    // A -= C^T D^-1 B, following StandardWellEquations::extract().
    // B and C have one row of blocks, with a nonzero at (0, j) only if the
    // well has a connection in cell j.
    for (auto colC = duneC_[0].begin(), endC = duneC_[0].end(); colC != endC; ++colC) {
        const auto row_index = cells_[colC.index()];
        for (auto colB = duneB_[0].begin(), endB = duneB_[0].end(); colB != endB; ++colB) {
            const auto col_index = cells_[colB.index()];
            // tmp = D^-1 B
            OffDiagMatrixBlockWellType tmp;
            detail::multMatrixImpl(invDuneD_[0][0], *colB, tmp, std::true_type());
            // block = -scale * C^T tmp, with the same unit conversion as the
            // residual correction in apply()
            typename SparseMatrixAdapter::MatrixBlock tmpMat;
            detail::negativeMultMatrixTransposed(*colC, tmp, tmpMat);
            tmpMat *= res_scales_[colC.index()];
            jacobian.addToBlock(row_index, col_index, tmpMat);
        }
    }
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
