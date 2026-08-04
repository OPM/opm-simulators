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
#ifndef OPM_SYSTEMCPRWPRESSURESTAGE_HEADER_INCLUDED
#define OPM_SYSTEMCPRWPRESSURESTAGE_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <dune/istl/paamg/pinfo.hh>

#include <opm/simulators/linalg/MatrixMarketSpecializations.hpp>

#include <dune/istl/matrixmarket.hh>

#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace Opm
{

// --------------------------------------------------------------------------
// CPRW pressure stage for the coupled (reservoir, well) system.
//
// Builds and solves the scalar pressure system of dimension
//
//     Nres + nWells
//
// obtained by contracting the full system matrix
//
//     S = [ A  C ]
//         [ B  D ]
//
// with a restriction R and a prolongation P:
//
//     R = blockdiag( w0^T ;  sum_{wb in well j} w1[wb]^T )
//     P = blockdiag( e_p  ;  e_q placed on the top block row of well j )
//
// where w0/w1 are the reservoir/well weights supplied from the outer layer,
// p is the reservoir pressure variable and q is the well pressure variable
// (bhp, or top segment pressure for a multisegment well).  Every well
// contributes exactly one coarse unknown, as in the classic CPRW.
//
// The point of doing it this way is that the assembly reads nothing but the
// four sparse blocks, the weights and the WellDofLayout.  No part of the well
// model is visible here.
//
// How the well part of the fine vectors takes part in the transfer.  The
// coarse matrix is the same in every case; only the vector transfers differ.
//
// The classic PressureBhpTransferPolicy has no well unknowns on the fine level
// at all (the wells are Schur-eliminated in the operator), so it can neither
// restrict a well residual nor apply a coarse bhp correction.  Selecting
// Classic here reproduces that, which makes the system solver and the classic
// cprw differ only in numerics rather than in formulation.
enum class WellTransfer
{
    Full,           // restrict the well residual and prolong the bhp correction
    NoProlongation, // restrict the well residual, discard the bhp correction
    Classic,        // neither -- as in PressureBhpTransferPolicy
};

// How the coarse diagonal of a well equation is formed.  Classic CPRW uses
// two different conventions: StandardWellEquations contracts D, while
// MultisegmentWellEquations sets the diagonal to minus the sum of the well
// row's reservoir entries and never reads D at all.
enum class WellCoarseDiagonal
{
    Auto,      // contract D for single-block wells, row sum for multi-block: as classic
    ContractD, // always contract D
    RowSum,    // always minus the row sum
};

inline WellCoarseDiagonal wellCoarseDiagonalFromString(const std::string& name)
{
    if (name == "auto") { return WellCoarseDiagonal::Auto; }
    if (name == "contract_d") { return WellCoarseDiagonal::ContractD; }
    if (name == "row_sum") { return WellCoarseDiagonal::RowSum; }
    OPM_THROW(std::invalid_argument,
              "Unknown well_coarse_diagonal '" + name
                  + "'. Valid values are 'auto', 'contract_d' and 'row_sum'.");
}

inline WellTransfer wellTransferFromString(const std::string& name)
{
    if (name == "full") {
        return WellTransfer::Full;
    }
    if (name == "no_prolongation") {
        return WellTransfer::NoProlongation;
    }
    if (name == "classic") {
        return WellTransfer::Classic;
    }
    OPM_THROW(std::invalid_argument,
              "Unknown well_transfer '" + name
                  + "'. Valid values are 'full', 'no_prolongation' and 'classic'.");
}

// --------------------------------------------------------------------------
template <class Scalar, class Comm = Dune::Amg::SequentialInformation>
class SystemCprwPressureStage
{
public:
    static constexpr bool isParallel = !std::is_same_v<Comm, Dune::Amg::SequentialInformation>;

    using CoarseOperator = Details::CoarseOperatorType<Scalar, Comm>;
    using CoarseMatrix = typename CoarseOperator::matrix_type;
    using CoarseVector = Details::PressureVectorType<Scalar>;
    using CoarseSolver = Dune::FlexibleSolver<CoarseOperator>;

    SystemCprwPressureStage(const SystemMatrix<Scalar>& S,
                            const PropertyTree& coarseSolverPrm,
                            const int pressureIndex,
                            const WellTransfer wellTransfer = WellTransfer::Full,
                            const Comm* comm = nullptr,
                            const WellCoarseDiagonal diagonal = WellCoarseDiagonal::ContractD,
                            const int verbosity = 0)
        : S_(S)
        , prm_(coarseSolverPrm)
        , pressureIndex_(pressureIndex)
        , wellTransfer_(wellTransfer)
        , comm_(comm)
        , diagonal_(diagonal)
        , verbosity_(verbosity)
    {
    }

    // Whether a coarse correction reaches the well unknowns at all.  The
    // caller can skip the C and D defect updates when it does not.
    bool prolongatesWellPressure() const
    {
        return wellTransfer_ == WellTransfer::Full;
    }

    // (Re)create the coarse sparsity pattern, communication and entries.
    // Separate from the coarse solver so that the assembly can be exercised on
    // its own.
    void buildCoarseSystem(const SystemVector<Scalar>& weights)
    {
        OPM_TIMEBLOCK(systemCprwBuildCoarseSystem);
        const auto& layout = wellLayout();
        const std::size_t numRes = S_.A->N();
        const std::size_t numWells = layout.numWells();

        buildCoarsePattern(numRes, numWells);
        buildCoarseCommunication(numWells);
        assembleCoarseMatrix(weights);

        coarseRhs_.resize(coarseMatrix_->N());
        coarseSol_.resize(coarseMatrix_->M());
        dumpCoarseMatrix();
    }

    // (Re)create the coarse system and the solver acting on it.  Must be
    // called whenever the well structure changes.
    void buildStructure(const SystemVector<Scalar>& weights)
    {
        OPM_TIMEBLOCK(systemCprwBuildStructure);
        buildCoarseSystem(weights);

        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
        OperatorArgs oargs(coarseMatrix_, *coarseComm_);
        coarseOperator_ = Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs);

        std::function<CoarseVector()> noWeights;
        if constexpr (isParallel) {
            coarseSolver_ = std::make_unique<CoarseSolver>(*coarseOperator_, *coarseComm_,
                                                           prm_, noWeights, /*pressureIndex=*/1);
        } else {
            coarseSolver_ = std::make_unique<CoarseSolver>(*coarseOperator_, prm_,
                                                           noWeights, /*pressureIndex=*/1);
        }
    }

    // Recompute the coarse entries from the current system matrix and weights.
    void update(const SystemVector<Scalar>& weights)
    {
        OPM_TIMEBLOCK(systemCprwUpdate);
        assembleCoarseMatrix(weights);
        dumpCoarseMatrix();
        coarseSolver_->preconditioner().update();
    }

    // One pressure-stage application:  restrict, coarse solve, prolong.
    void apply(const ResVector<Scalar>& dRes,
               const WellVector<Scalar>& dWell,
               const SystemVector<Scalar>& weights,
               ResVector<Scalar>& vRes,
               WellVector<Scalar>& vWell)
    {
        OPM_TIMEBLOCK(systemCprwApply);
        moveToCoarseLevel(dRes, dWell, weights);

        dumpCoarseRhs();
        coarseSol_ = 0.0;
        Dune::InverseOperatorResult result;
        coarseSolver_->apply(coarseSol_, coarseRhs_, result);

        moveToFineLevel(vRes, vWell);
    }

    // Restriction:  coarseRhs = R * (dRes, dWell).  Unlike the classic CPRW
    // transfer policy, the well residual really is carried to the coarse
    // level rather than dropped.
    void moveToCoarseLevel(const ResVector<Scalar>& dRes,
                           const WellVector<Scalar>& dWell,
                           const SystemVector<Scalar>& weights)
    {
        restrictInto(dRes, dWell, weights, coarseRhs_);
    }

    // Prolongation:  (vRes, vWell) = P * coarseSol.  The coarse well
    // correction lands on the pressure unknown of each well's top block --
    // the classic CPRW throws it away instead.
    void moveToFineLevel(ResVector<Scalar>& vRes, WellVector<Scalar>& vWell) const
    {
        prolongFrom(coarseSol_, vRes, vWell);
    }

    const CoarseMatrix& coarseMatrix() const
    {
        return *coarseMatrix_;
    }

    // Handles needed when the coarse level is driven from outside, e.g. by
    // SystemPressureBhpTransferPolicy.
    const std::shared_ptr<CoarseMatrix>& coarseMatrixPtr() const
    {
        return coarseMatrix_;
    }

    const Comm& coarseCommunication() const
    {
        return *coarseComm_;
    }

    void assembleCoarseEntries(const SystemVector<Scalar>& weights)
    {
        assembleCoarseMatrix(weights);
    }

    // Transfer forms writing into a caller-supplied coarse vector, so that a
    // transfer policy can use its own storage without an extra copy.
    void moveToCoarseLevel(const ResVector<Scalar>& dRes,
                           const WellVector<Scalar>& dWell,
                           const SystemVector<Scalar>& weights,
                           CoarseVector& out) const
    {
        restrictInto(dRes, dWell, weights, out);
    }

    void moveToFineLevel(const CoarseVector& in,
                         ResVector<Scalar>& vRes,
                         WellVector<Scalar>& vWell) const
    {
        prolongFrom(in, vRes, vWell);
    }

    const CoarseVector& coarseRhs() const
    {
        return coarseRhs_;
    }

    CoarseVector& coarseSolution()
    {
        return coarseSol_;
    }

private:
    const WellDofLayout& wellLayout() const
    {
        if (S_.wellLayout == nullptr) {
            OPM_THROW(std::logic_error,
                      "SystemCprwPressureStage requires a WellDofLayout on the system matrix. "
                      "It is filled by ISTLSolverSystem; a null layout means the CPRW pressure "
                      "stage was constructed outside that path.");
        }
        return *S_.wellLayout;
    }

    // Pattern: the reservoir block keeps A's pattern; each well j adds one row
    // and one column, coupled to exactly the cells its B/C rows touch, plus a
    // diagonal.  The merged D is block diagonal by well, so the well-well part
    // of the coarse system is diagonal.
    void buildCoarsePattern(const std::size_t numRes, const std::size_t numWells)
    {
        const auto& A = *S_.A;
        const auto& B = *S_.B;
        const auto& C = *S_.C;
        const auto& layout = wellLayout();

        const std::size_t dim = numRes + numWells;
        const std::size_t averageEntriesPerRow
            = static_cast<std::size_t>(std::ceil(static_cast<double>(A.nonzeroes()) / A.N()));
        const double overflowFraction = 1.2;
        coarseMatrix_ = std::make_shared<CoarseMatrix>(dim, dim,
                                                       averageEntriesPerRow,
                                                       overflowFraction,
                                                       CoarseMatrix::implicit);

        // Reservoir-reservoir: A's pattern.
        for (auto row = A.begin(), rowEnd = A.end(); row != rowEnd; ++row) {
            for (auto col = row->begin(), colEnd = row->end(); col != colEnd; ++col) {
                coarseMatrix_->entry(row.index(), col.index()) = 0.0;
            }
        }

        for (std::size_t j = 0; j < numWells; ++j) {
            const std::size_t wdof = numRes + j;
            coarseMatrix_->entry(wdof, wdof) = 0.0;

            // Well row: the cells reached by any block row of this well.
            for (std::size_t wb = layout.firstBlock(j); wb < layout.endBlock(j); ++wb) {
                for (auto col = B[wb].begin(), colEnd = B[wb].end(); col != colEnd; ++col) {
                    coarseMatrix_->entry(wdof, col.index()) = 0.0;
                }
            }
        }

        // Well column: every cell whose C row references any block of the well.
        for (std::size_t c = 0; c < numRes; ++c) {
            for (auto col = C[c].begin(), colEnd = C[c].end(); col != colEnd; ++col) {
                const auto j = wellOfBlock(col.index());
                if (j.has_value()) {
                    coarseMatrix_->entry(c, numRes + *j) = 0.0;
                }
            }
        }

        coarseMatrix_->compress();
    }

    void buildCoarseCommunication([[maybe_unused]] const std::size_t numWells)
    {
        if constexpr (isParallel) {
            coarseComm_ = std::make_shared<Comm>(comm_->communicator(), comm_->category(), false);
            // Well DOFs are rank local and owned, appended after the reservoir
            // DOFs -- the same convention the classic CPRW coarse system uses.
            extendCommunicatorWithWells(*comm_, coarseComm_, static_cast<int>(numWells));
        } else {
            coarseComm_ = std::make_shared<Comm>();
        }
    }

    void assembleCoarseMatrix(const SystemVector<Scalar>& weights)
    {
        using namespace Dune::Indices;
        OPM_TIMEBLOCK(systemCprwAssemble);

        const auto& A = *S_.A;
        const auto& B = *S_.B;
        const auto& C = *S_.C;
        const auto& D = *S_.D;
        const auto& layout = wellLayout();
        const auto& w0 = weights[_0];
        const auto& w1 = weights[_1];

        const std::size_t numRes = A.N();
        const int p = pressureIndex_;
        const int q = layout.pressureDofIndex;

        *coarseMatrix_ = 0.0;

        // Reservoir rows, reservoir columns:  sum_i w0[c][i] * A[c][c'][i][p]
        for (auto row = A.begin(), rowEnd = A.end(); row != rowEnd; ++row) {
            const auto& bw = w0[row.index()];
            for (auto col = row->begin(), colEnd = row->end(); col != colEnd; ++col) {
                Scalar el = 0.0;
                for (std::size_t i = 0; i < bw.size(); ++i) {
                    el += (*col)[i][p] * bw[i];
                }
                (*coarseMatrix_)[row.index()][col.index()] = el;
            }
        }

        // Reservoir rows, well columns:  sum over every block row of well j of
        //   sum_i w0[c][i] * C[c][wb][i][q].
        //
        // Summing over all of a well's blocks is the Galerkin column for a
        // prolongation that spreads a well's coarse unknown over all of its
        // segment pressures, which is what MultisegmentWellEquations::
        // extractCPRPressureMatrix does (it accumulates over every segment
        // row).  Taking the top block alone instead loses every segment but
        // the first: on Norne with one segment per connection that is most of
        // the well, and it is what made the coarse system far weaker than the
        // classic cprw one for multisegment wells.
        for (std::size_t c = 0; c < numRes; ++c) {
            const auto& bw = w0[c];
            for (auto col = C[c].begin(), colEnd = C[c].end(); col != colEnd; ++col) {
                const auto j = wellOfBlock(col.index());
                if (!j.has_value() || layout.isPressureControlled(*j)) {
                    continue;
                }
                Scalar el = 0.0;
                for (std::size_t i = 0; i < bw.size(); ++i) {
                    el += (*col)[i][q] * bw[i];
                }
                (*coarseMatrix_)[c][numRes + *j] += el;
            }
        }

        // Well rows.
        for (std::size_t j = 0; j < layout.numWells(); ++j) {
            const std::size_t wdof = numRes + j;
            if (layout.isPressureControlled(j)) {
                // A pressure-controlled well has a trivial coarse equation:
                // its bhp is prescribed, so the coarse system carries dp = 0
                // rather than a contracted well equation.
                (*coarseMatrix_)[wdof][wdof] = 1.0;
                continue;
            }
            const bool rowSumDiag
                = (diagonal_ == WellCoarseDiagonal::RowSum)
                || (diagonal_ == WellCoarseDiagonal::Auto
                    && layout.endBlock(j) - layout.firstBlock(j) > 1);
            Scalar rowSum = 0.0;

            for (std::size_t wb = layout.firstBlock(j); wb < layout.endBlock(j); ++wb) {
                const auto& lw = w1[wb];

                // Well row, reservoir columns:
                //   sum_{wb in j} sum_i w1[wb][i] * B[wb][c][i][p]
                for (auto col = B[wb].begin(), colEnd = B[wb].end(); col != colEnd; ++col) {
                    Scalar el = 0.0;
                    for (std::size_t i = 0; i < lw.size(); ++i) {
                        el += lw[i] * (*col)[i][p];
                    }
                    (*coarseMatrix_)[wdof][col.index()] += el;
                    rowSum += el;
                }

                if (rowSumDiag) {
                    // The classic multisegment convention takes the diagonal
                    // from the row sum and never reads D.
                    continue;
                }

                // Well row, well columns:
                //   sum_{wb in j} sum_{wb' in k} sum_i w1[wb][i] * D[wb][wb'][i][q]
                // Summed over all of well k's blocks, to match the column
                // convention above.  Because the merged D is block diagonal by
                // well this only ever contributes to k == j, but it now picks
                // up the full segment-to-segment coupling rather than just the
                // top segment's column.
                for (auto col = D[wb].begin(), colEnd = D[wb].end(); col != colEnd; ++col) {
                    const auto k = wellOfBlock(col.index());
                    if (!k.has_value()) {
                        continue;
                    }
                    Scalar el = 0.0;
                    for (std::size_t i = 0; i < lw.size(); ++i) {
                        el += lw[i] * (*col)[i][q];
                    }
                    (*coarseMatrix_)[wdof][numRes + *k] += el;
                }
            }

            // A well whose contraction cancels exactly would leave a zero on
            // the coarse diagonal and make the pressure system singular. That
            // can happen for a well with no local perforations, or when the
            // weights annihilate the pressure column. Regularise to a unit row
            // rather than handing a singular system to AMG.
            if (rowSumDiag) {
                (*coarseMatrix_)[wdof][wdof] = -rowSum;
            }
            auto& diag = (*coarseMatrix_)[wdof][wdof][0][0];
            if (!(std::abs(diag) > 0.0)) {
                diag = 1.0;
            }
        }
    }

    // Restriction:  out = R * (dRes, dWell).  Unlike the classic CPRW transfer
    // policy the well residual really is carried to the coarse level, unless
    // the classic transfer was asked for.
    void restrictInto(const ResVector<Scalar>& dRes,
                      const WellVector<Scalar>& dWell,
                      const SystemVector<Scalar>& weights,
                      CoarseVector& out) const
    {
        using namespace Dune::Indices;
        const auto& layout = wellLayout();
        const auto& w0 = weights[_0];
        const auto& w1 = weights[_1];
        const std::size_t numRes = dRes.size();

        out = 0.0;

        for (std::size_t c = 0; c < numRes; ++c) {
            const auto& bw = w0[c];
            Scalar el = 0.0;
            for (std::size_t i = 0; i < bw.size(); ++i) {
                el += dRes[c][i] * bw[i];
            }
            out[c] = el;
        }

        if (wellTransfer_ == WellTransfer::Classic) {
            // The classic policy leaves the well rows of the coarse right-hand
            // side at zero; keep them zero so that the two formulations agree.
            return;
        }

        for (std::size_t j = 0; j < layout.numWells(); ++j) {
            Scalar el = 0.0;
            for (std::size_t wb = layout.firstBlock(j); wb < layout.endBlock(j); ++wb) {
                const auto& lw = w1[wb];
                for (std::size_t i = 0; i < lw.size(); ++i) {
                    el += lw[i] * dWell[wb][i];
                }
            }
            out[numRes + j] = el;
        }
    }

    // Prolongation:  (vRes, vWell) = P * in.  The coarse well correction lands
    // on the pressure unknown of each well's top block; the classic policy
    // throws it away instead.
    void prolongFrom(const CoarseVector& in,
                     ResVector<Scalar>& vRes,
                     WellVector<Scalar>& vWell) const
    {
        const auto& layout = wellLayout();
        const std::size_t numRes = vRes.size();
        const int q = layout.pressureDofIndex;

        vRes = 0.0;
        for (std::size_t c = 0; c < numRes; ++c) {
            vRes[c][pressureIndex_] = in[c][0];
        }

        vWell = 0.0;
        if (!prolongatesWellPressure()) {
            // The coarse bhp correction is computed but discarded; it acts only
            // through its influence on the reservoir pressure, as in the
            // classic policy.  A following well solve corrects the wells.
            return;
        }
        // Spread the well's coarse value over all of its segment pressures by
        // a constant.  This is the P the coarse matrix is assembled for; only
        // the segment pressures are set, the other well unknowns (rates,
        // compositions) are left to the well solve that follows.
        for (std::size_t j = 0; j < layout.numWells(); ++j) {
            for (std::size_t wb = layout.firstBlock(j); wb < layout.endBlock(j); ++wb) {
                vWell[wb][q] = in[numRes + j][0];
            }
        }
    }

    // Same convention as the classic path: verbosity above 10 writes the
    // coarse system out so the two can be compared entry by entry.
    void dumpCoarseMatrix() const
    {
        if (verbosity_ <= 10) {
            return;
        }
        static int counter = 0;
        std::ofstream out("system_cprw_coarse_" + std::to_string(counter++) + ".mm");
        if (out) {
            Dune::writeMatrixMarket(*coarseMatrix_, out);
        }
    }

    void dumpCoarseRhs() const
    {
        if (verbosity_ <= 10) {
            return;
        }
        static int counter = 0;
        std::ofstream out("system_cprw_rhs_" + std::to_string(counter++) + ".mm");
        if (out) {
            Dune::writeMatrixMarket(coarseRhs_, out);
        }
    }

    // Merged well block row -> well index.  Linear scan is fine: the offsets
    // are sorted and this is only used while walking sparse rows of the well
    // blocks, which are short.
    std::optional<std::size_t> wellOfBlock(const std::size_t blockRow) const
    {
        const auto& offsets = wellLayout().wellBlockOffsets;
        const auto it = std::upper_bound(offsets.begin(), offsets.end(), blockRow);
        if (it == offsets.begin() || it == offsets.end()) {
            return std::nullopt;
        }
        return static_cast<std::size_t>(std::distance(offsets.begin(), it) - 1);
    }

    const SystemMatrix<Scalar>& S_;
    PropertyTree prm_;
    int pressureIndex_ = 0;
    WellTransfer wellTransfer_ = WellTransfer::Full;
    const Comm* comm_ = nullptr;
    WellCoarseDiagonal diagonal_ = WellCoarseDiagonal::ContractD;
    int verbosity_ = 0;

    std::shared_ptr<Comm> coarseComm_;
    std::shared_ptr<CoarseMatrix> coarseMatrix_;
    std::shared_ptr<CoarseOperator> coarseOperator_;
    std::unique_ptr<CoarseSolver> coarseSolver_;

    CoarseVector coarseRhs_;
    CoarseVector coarseSol_;
};

} // namespace Opm

#endif // OPM_SYSTEMCPRWPRESSURESTAGE_HEADER_INCLUDED
