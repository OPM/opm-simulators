/*
Copyright Equinor ASA 2026

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

#ifndef OPM_SYSTEMMATRIXWELLOPERATOR_HEADER_INCLUDED
#define OPM_SYSTEMMATRIXWELLOPERATOR_HEADER_INCLUDED

#include "SystemTypes.hpp"

#include <opm/common/Exceptions.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace Opm
{

// LinearOperatorExtra bridge that exposes the B/C/D sub-blocks of a
// SystemMatrix through the interface expected by WellModelMatrixAdapter and
// PressureBhpTransferPolicy (used by the CPRW preconditioner).
//
// Follows SystemMatrix's block naming (same as SystemPreconditioner::apply()):
//   S_.B  =  A_wr  (well-row, res-col, 4x3 blocks)  == S[_1][_0]
//   S_.C  =  A_rw  (res-row,  well-col, 3x4 blocks)  == S[_0][_1]
//   S_.D  =  A_ww  (well-well, 4x4 blocks)           == S[_1][_1]
//
// The apply() method is intentionally a no-op: well DOFs are updated
// explicitly by wellSolver_ in SystemPreconditioner, so the CPRW fine-level
// smoother (operating through WellModelMatrixAdapter) only sees A_rr.
template<class Scalar>
class SystemMatrixWellOperator
    : public LinearOperatorExtra<ResVector<Scalar>, ResVector<Scalar>>
{
public:
    using X = ResVector<Scalar>;
    using Y = ResVector<Scalar>;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<Scalar, 1, 1>>;

    SystemMatrixWellOperator(const SystemMatrix<Scalar>& S, int pressureIndex)
        : S_(S)
        , pressureIndex_(pressureIndex)
    {}

    // -----------------------------------------------------------------------
    // Dune::LinearOperator interface
    // -----------------------------------------------------------------------

    void apply(const X&, Y&) const override {}

    void applyscaleadd(Scalar, const X&, Y&) const override {}

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    // -----------------------------------------------------------------------
    // LinearOperatorExtra interface
    // -----------------------------------------------------------------------

    int getNumberOfExtraEquations() const override
    {
        return static_cast<int>(S_.D->N());
    }

    // Add sparsity entries for the well rows/columns of the coarse pressure
    // matrix.  Called once during coarse-matrix structure setup.
    // Called while the coarse pressure matrix is still in implicit build mode
    // (PressureBhpTransferPolicy calls compress() after this).  Must use
    // entry() to insert new entries rather than operator[][], which only works
    // on already-compressed rows.
    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const override
    {
        const std::size_t nCells = S_.A->N();

        // A_wr (S_.B): well-row, res-col  →  jacobian[nCells+w][cell]
        for (auto rowIt = S_.B->begin(); rowIt != S_.B->end(); ++rowIt) {
            const std::size_t well_w = rowIt.index();
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                jacobian.entry(nCells + well_w, colIt.index()) = 0.0;
            }
        }

        // A_rw (S_.C): res-row, well-col  →  jacobian[cell][nCells+w]
        for (auto rowIt = S_.C->begin(); rowIt != S_.C->end(); ++rowIt) {
            const std::size_t cell = rowIt.index();
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                jacobian.entry(cell, nCells + colIt.index()) = 0.0;
            }
        }

        // A_ww (S_.D): well-row, well-col  →  jacobian[nCells+w][nCells+w']
        for (auto rowIt = S_.D->begin(); rowIt != S_.D->end(); ++rowIt) {
            const std::size_t well_w = rowIt.index();
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                jacobian.entry(nCells + well_w, nCells + colIt.index()) = 0.0;
            }
        }
    }

    // Fill the well rows/columns of the coarse pressure matrix.
    // Implements the CPRW extended pressure system from the paper:
    //   A_cprw = [C_p^T W A_rr C_p     C_p^T W A_rw C_bhp  ]
    //            [C_bhp^T Ww A_wr C_p   C_bhp^T Ww A_ww C_bhp]
    //
    // use_well_weights=false: well weights = cell IMPES weights.
    // use_well_weights=true:  well weights = inv(D)^T * e_bhp (quasi-IMPES),
    //                         matching the formula in StandardWellEquations.
    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const X& weights,
                                  bool use_well_weights) const override
    {
        const std::size_t nCells = S_.A->N();
        constexpr int bhpIdx = numWellDofs - 1;

        // Upper-right block: jacobian[cell][nCells + well_w]
        // = sum_k  A_rw[cell][well_w][k][bhpIdx] * weights[cell][k]   (same for both modes)
        for (auto rowIt = S_.C->begin(); rowIt != S_.C->end(); ++rowIt) {
            const std::size_t cell = rowIt.index();
            const auto& wt = weights[cell];
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                const auto& block = *colIt;  // 3x4 block
                Scalar val = 0;
                for (int k = 0; k < numResDofs; ++k) {
                    val += block[k][bhpIdx] * wt[k];
                }
                jacobian[cell][nCells + colIt.index()] = val;
            }
        }

        // Pre-compute per-well weights (bwt) and diagonal values.
        // use_well_weights=false: bwt from first connected cell's IMPES weights
        // use_well_weights=true:  bwt = inv(D[w][w])^T[:,bhp] normalised
        const std::size_t nWells = S_.D->N();
        using BwtArray = std::array<Scalar, numWellDofs>;
        std::vector<BwtArray> wellBwt(nWells);
        std::vector<Scalar>   wellDiag(nWells, 1.0);

        for (auto rowIt = S_.D->begin(); rowIt != S_.D->end(); ++rowIt) {
            const std::size_t w = rowIt.index();
            BwtArray& bwt = wellBwt[w];
            bwt.fill(Scalar(0));

            const auto setDefaultWellWeights = [&]() {
                // Use averaged reservoir IMPES weights as a robust fallback when
                // quasi-IMPES well weights are unavailable.
                std::vector<Scalar> totalWeights(numResDofs, 0.0);
                int nperf = 0;
                for (auto colIt = (*S_.B)[w].begin(); colIt != (*S_.B)[w].end(); ++colIt) {
                    const auto& wt = weights[colIt.index()];
                    for (int k = 0; k < numResDofs; ++k) {
                        totalWeights[k] += wt[k];
                    }
                    nperf++;
                }
                if (nperf > 0) {
                    for (int k = 0; k < numResDofs; ++k) {
                        bwt[k] = totalWeights[k] / nperf;
                    }
                } else {
                    // no perforations on this rank: use default weights for regularization
                    for (int k = 0; k < numResDofs; ++k) {
                        bwt[k] = Scalar(1);
                    }
                }
                // diagonal: sum_k D[w][w][k][bhp] * bwt[k]  (k < numWellDofs-1)
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                    if (colIt.index() == w) {
                        const auto& block = *colIt;
                        Scalar val = 0;
                        for (int k = 0; k < numWellDofs - 1; ++k) {
                            val += block[k][bhpIdx] * bwt[k];
                        }
                        wellDiag[w] = (val != 0) ? val : Scalar(1);
                        break;
                    }
                }
            };

            if (use_well_weights) {
                // bwt = inv(D[w][w])^T * e_bhp  →  row bhp of inv(D[w][w])
                try {
                    auto invD = (*S_.D)[w][w];   // 4x4 copy
                    detail::invertMatrix(invD);
                    Scalar abs_max = 0;
                    for (int i = 0; i < numWellDofs; ++i) {
                        bwt[i] = invD[bhpIdx][i];
                        abs_max = std::max(abs_max, std::abs(bwt[i]));
                    }
                    if (abs_max > 0) {
                        for (int i = 0; i < numWellDofs; ++i) {
                            bwt[i] /= abs_max;
                        }
                        wellDiag[w] = Scalar(1) / abs_max;
                    } else {
                        wellDiag[w] = Scalar(1);
                    }
                } catch (const NumericalProblem&) {
                    setDefaultWellWeights();
                }
            } else {
                setDefaultWellWeights();
            }

            // Set diagonal entry
            jacobian[nCells + w][nCells + w] = wellDiag[w];
        }

        // Lower-left block: jacobian[nCells + well_w][cell]
        // = sum_k  A_wr[well_w][cell][k][pressureIndex_] * bwt[well_w][k]
        for (auto rowIt = S_.B->begin(); rowIt != S_.B->end(); ++rowIt) {
            const std::size_t well_w = rowIt.index();
            const BwtArray& bwt = wellBwt[well_w];
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                const auto& block = *colIt;  // 4x3 block
                Scalar val = 0;
                for (int k = 0; k < numWellDofs; ++k) {
                    val += block[k][pressureIndex_] * bwt[k];
                }
                jacobian[nCells + well_w][colIt.index()] = val;
            }
        }
    }

private:
    const SystemMatrix<Scalar>& S_;
    int pressureIndex_;
};

} // namespace Opm

#endif // OPM_SYSTEMMATRIXWELLOPERATOR_HEADER_INCLUDED
