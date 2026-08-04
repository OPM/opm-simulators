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
#include <config.h>
#define BOOST_TEST_MODULE OPM_test_SystemCprwPressureStage
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/system/SystemCprwPressureStage.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

using Scalar = double;
using RRMatrix = Opm::RRMatrix<Scalar>;
using RWMatrix = Opm::RWMatrix<Scalar>;
using WRMatrix = Opm::WRMatrix<Scalar>;
using WWMatrix = Opm::WWMatrix<Scalar>;
using SystemVector = Opm::SystemVector<Scalar>;
using Stage = Opm::SystemCprwPressureStage<Scalar>;

constexpr int numRes = Opm::numResDofs;   // 3
constexpr int numWell = Opm::numWellDofs; // 4
constexpr int pressureIndex = 0;
constexpr int wellPressureIndex = numWell - 1;

// The test fixture below models 4 reservoir cells and 2 wells:
//   well 0 - a standard well, 1 block row, perforating cells 0 and 1
//   well 1 - a multisegment well, 3 block rows (segments), perforating
//            cells 2 and 3 from different segments
// so that the per-well aggregation over block rows is actually exercised.
constexpr std::size_t numCells = 4;
constexpr std::size_t numWellBlocks = 4; // 1 + 3

struct BlockSpec
{
    std::size_t column;
    Scalar base;
};

using MatrixPattern = std::vector<std::vector<BlockSpec>>;

// Distinct, non-symmetric values so that a transposed or mis-indexed
// contraction cannot accidentally produce the right answer.
template <class Block>
Block makeBlock(const Scalar base)
{
    Block block;
    for (int row = 0; row < Block::rows; ++row) {
        for (int col = 0; col < Block::cols; ++col) {
            block[row][col] = base + 3.0 * row + 7.0 * col + 0.25 * row * col;
        }
    }
    return block;
}

template <class Matrix>
Matrix buildMatrix(const std::size_t rows, const std::size_t cols, const MatrixPattern& pattern)
{
    std::size_t nonzeroes = 0;
    for (const auto& row : pattern) {
        nonzeroes += row.size();
    }

    Matrix matrix(rows, cols, nonzeroes, Matrix::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (const auto& e : pattern[row.index()]) {
            row.insert(e.column);
        }
    }
    for (std::size_t row = 0; row < rows; ++row) {
        for (const auto& e : pattern[row]) {
            matrix[row][e.column] = makeBlock<typename Matrix::block_type>(e.base);
        }
    }
    return matrix;
}

struct Fixture
{
    RRMatrix A;
    RWMatrix C;
    WRMatrix B;
    WWMatrix D;
    Opm::WellDofLayout layout;
    Opm::SystemMatrix<Scalar> S;
    SystemVector weights;

    Fixture()
    {
        // A: tridiagonal over the 4 cells.
        A = buildMatrix<RRMatrix>(numCells, numCells,
                                  {{{0, 1.0}, {1, 2.0}},
                                   {{0, 3.0}, {1, 4.0}, {2, 5.0}},
                                   {{1, 6.0}, {2, 7.0}, {3, 8.0}},
                                   {{2, 9.0}, {3, 10.0}}});

        // C: cell -> well block.  Cells 0,1 see well 0 (block 0); cells 2,3
        // see well 1 through segments 1 and 2 (blocks 2 and 3).  Only the top
        // block of each well carries a coarse unknown, so the entries on
        // blocks 2 and 3 must be ignored by the assembly.
        C = buildMatrix<RWMatrix>(numCells, numWellBlocks,
                                  {{{0, 11.0}},
                                   {{0, 12.0}},
                                   {{1, 13.0}, {2, 14.0}},
                                   {{1, 15.0}, {3, 16.0}}});

        // B: well block -> cell.
        B = buildMatrix<WRMatrix>(numWellBlocks, numCells,
                                  {{{0, 17.0}, {1, 18.0}},
                                   {{2, 19.0}},
                                   {{2, 20.0}},
                                   {{3, 21.0}}});

        // D: block diagonal by well.  Well 0 is a single block; well 1 is a
        // 3x3 segment coupling with blocks 1..3.
        D = buildMatrix<WWMatrix>(numWellBlocks, numWellBlocks,
                                  {{{0, 22.0}},
                                   {{1, 23.0}, {2, 24.0}},
                                   {{1, 25.0}, {2, 26.0}, {3, 27.0}},
                                   {{2, 28.0}, {3, 29.0}}});

        layout.wellBlockOffsets = {0, 1, 4};
        layout.pressureDofIndex = wellPressureIndex;

        S.A = &A;
        S.C = &C;
        S.B = &B;
        S.D = &D;
        S.wellLayout = &layout;

        weights[Dune::Indices::_0].resize(numCells);
        for (std::size_t c = 0; c < numCells; ++c) {
            for (int i = 0; i < numRes; ++i) {
                weights[Dune::Indices::_0][c][i] = 0.5 + 0.1 * c + 0.3 * i;
            }
        }
        weights[Dune::Indices::_1].resize(numWellBlocks);
        for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
            for (int i = 0; i < numWell; ++i) {
                weights[Dune::Indices::_1][wb][i] = 0.2 + 0.7 * wb - 0.15 * i;
            }
        }
    }
};

// Independent brute-force reference: densify S, build R and P explicitly and
// form R*S*P.  Deliberately written without reusing any production helper.
std::vector<std::vector<Scalar>> referenceCoarseMatrix(const Fixture& f)
{
    const std::size_t nWells = f.layout.numWells();
    const std::size_t fineDim = numCells * numRes + numWellBlocks * numWell;
    const std::size_t coarseDim = numCells + nWells;

    // Dense fine system.
    std::vector<std::vector<Scalar>> S(fineDim, std::vector<Scalar>(fineDim, 0.0));
    const auto resOff = [](const std::size_t c, const int i) { return c * numRes + i; };
    const auto wellOff = [](const std::size_t wb, const int i) {
        return numCells * numRes + wb * numWell + i;
    };

    for (auto row = f.A.begin(); row != f.A.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (int i = 0; i < numRes; ++i) {
                for (int j = 0; j < numRes; ++j) {
                    S[resOff(row.index(), i)][resOff(col.index(), j)] = (*col)[i][j];
                }
            }
        }
    }
    for (auto row = f.C.begin(); row != f.C.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (int i = 0; i < numRes; ++i) {
                for (int j = 0; j < numWell; ++j) {
                    S[resOff(row.index(), i)][wellOff(col.index(), j)] = (*col)[i][j];
                }
            }
        }
    }
    for (auto row = f.B.begin(); row != f.B.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (int i = 0; i < numWell; ++i) {
                for (int j = 0; j < numRes; ++j) {
                    S[wellOff(row.index(), i)][resOff(col.index(), j)] = (*col)[i][j];
                }
            }
        }
    }
    for (auto row = f.D.begin(); row != f.D.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (int i = 0; i < numWell; ++i) {
                for (int j = 0; j < numWell; ++j) {
                    S[wellOff(row.index(), i)][wellOff(col.index(), j)] = (*col)[i][j];
                }
            }
        }
    }

    // R (coarseDim x fineDim) and P (fineDim x coarseDim).
    std::vector<std::vector<Scalar>> R(coarseDim, std::vector<Scalar>(fineDim, 0.0));
    std::vector<std::vector<Scalar>> P(fineDim, std::vector<Scalar>(coarseDim, 0.0));

    for (std::size_t c = 0; c < numCells; ++c) {
        for (int i = 0; i < numRes; ++i) {
            R[c][resOff(c, i)] = f.weights[Dune::Indices::_0][c][i];
        }
        P[resOff(c, pressureIndex)][c] = 1.0;
    }
    for (std::size_t j = 0; j < nWells; ++j) {
        for (std::size_t wb = f.layout.firstBlock(j); wb < f.layout.endBlock(j); ++wb) {
            for (int i = 0; i < numWell; ++i) {
                R[numCells + j][wellOff(wb, i)] = f.weights[Dune::Indices::_1][wb][i];
            }
        }
        P[wellOff(f.layout.firstBlock(j), wellPressureIndex)][numCells + j] = 1.0;
    }

    std::vector<std::vector<Scalar>> coarse(coarseDim, std::vector<Scalar>(coarseDim, 0.0));
    for (std::size_t r = 0; r < coarseDim; ++r) {
        for (std::size_t c = 0; c < coarseDim; ++c) {
            Scalar sum = 0.0;
            for (std::size_t k = 0; k < fineDim; ++k) {
                if (R[r][k] == 0.0) {
                    continue;
                }
                for (std::size_t l = 0; l < fineDim; ++l) {
                    if (P[l][c] != 0.0) {
                        sum += R[r][k] * S[k][l] * P[l][c];
                    }
                }
            }
            coarse[r][c] = sum;
        }
    }
    return coarse;
}

Scalar coarseEntry(const Stage& stage, const std::size_t row, const std::size_t col)
{
    const auto& m = stage.coarseMatrix();
    if (!m.exists(row, col)) {
        return 0.0;
    }
    return m[row][col][0][0];
}

} // anonymous namespace

// The coarse system must be exactly R*S*P.  This fails if the C or B
// contraction is dropped, if the wrong block row is taken as a well's coarse
// unknown, or if the reservoir/well pressure index is wrong.
BOOST_AUTO_TEST_CASE(CoarseMatrixEqualsRestrictedSystem)
{
    const Fixture f;
    Stage stage(f.S, Opm::PropertyTree(), pressureIndex);
    stage.buildCoarseSystem(f.weights);

    const auto expected = referenceCoarseMatrix(f);
    const std::size_t coarseDim = numCells + f.layout.numWells();

    BOOST_REQUIRE_EQUAL(stage.coarseMatrix().N(), coarseDim);
    BOOST_REQUIRE_EQUAL(stage.coarseMatrix().M(), coarseDim);

    for (std::size_t r = 0; r < coarseDim; ++r) {
        for (std::size_t c = 0; c < coarseDim; ++c) {
            BOOST_CHECK_CLOSE(coarseEntry(stage, r, c), expected[r][c], 1e-10);
        }
    }
}

// The well coupling must actually be present: every well row/column pair that
// the fixture perforates has to be non-zero.  Guards against a coarse system
// that is merely the reservoir block padded with an identity.
BOOST_AUTO_TEST_CASE(WellCouplingIsPresent)
{
    const Fixture f;
    Stage stage(f.S, Opm::PropertyTree(), pressureIndex);
    stage.buildCoarseSystem(f.weights);

    // Well 0 perforates cells 0 and 1, well 1 cells 2 and 3.
    const std::vector<std::vector<std::size_t>> perfs = {{0, 1}, {2, 3}};
    for (std::size_t j = 0; j < perfs.size(); ++j) {
        const std::size_t wdof = numCells + j;
        BOOST_CHECK_NE(coarseEntry(stage, wdof, wdof), 0.0);
        for (const auto c : perfs[j]) {
            BOOST_CHECK_NE(coarseEntry(stage, wdof, c), 0.0);
            BOOST_CHECK_NE(coarseEntry(stage, c, wdof), 0.0);
        }
    }

    // Cells belonging to one well must not couple to the other well.
    BOOST_CHECK_EQUAL(coarseEntry(stage, numCells + 0, 2), 0.0);
    BOOST_CHECK_EQUAL(coarseEntry(stage, numCells + 1, 0), 0.0);
}

// The reservoir block of the coarse system must be the plain CPR contraction,
// i.e. adding the wells must not perturb the reservoir rows.
BOOST_AUTO_TEST_CASE(ReservoirBlockIsPlainCprContraction)
{
    const Fixture f;
    Stage stage(f.S, Opm::PropertyTree(), pressureIndex);
    stage.buildCoarseSystem(f.weights);

    for (auto row = f.A.begin(); row != f.A.end(); ++row) {
        const auto& bw = f.weights[Dune::Indices::_0][row.index()];
        for (auto col = row->begin(); col != row->end(); ++col) {
            Scalar expected = 0.0;
            for (int i = 0; i < numRes; ++i) {
                expected += (*col)[i][pressureIndex] * bw[i];
            }
            BOOST_CHECK_CLOSE(coarseEntry(stage, row.index(), col.index()), expected, 1e-10);
        }
    }
}

// Restriction and prolongation must be transposes of each other in the sense
// that R*P is the identity when the weights select exactly the unknowns the
// prolongation writes.
BOOST_AUTO_TEST_CASE(RestrictOfProlongIsIdentityForSelectingWeights)
{
    Fixture f;
    // Weights that pick out precisely the prolonged components.
    for (std::size_t c = 0; c < numCells; ++c) {
        f.weights[Dune::Indices::_0][c] = 0.0;
        f.weights[Dune::Indices::_0][c][pressureIndex] = 1.0;
    }
    for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
        f.weights[Dune::Indices::_1][wb] = 0.0;
    }
    for (std::size_t j = 0; j < f.layout.numWells(); ++j) {
        f.weights[Dune::Indices::_1][f.layout.firstBlock(j)][wellPressureIndex] = 1.0;
    }

    Stage stage(f.S, Opm::PropertyTree(), pressureIndex);
    stage.buildCoarseSystem(f.weights);

    const std::size_t coarseDim = numCells + f.layout.numWells();
    auto& coarseSol = stage.coarseSolution();
    coarseSol.resize(coarseDim);
    for (std::size_t i = 0; i < coarseDim; ++i) {
        coarseSol[i] = 1.0 + 2.0 * i;
    }

    Opm::ResVector<Scalar> vRes(numCells);
    Opm::WellVector<Scalar> vWell(numWellBlocks);
    stage.moveToFineLevel(vRes, vWell);
    stage.moveToCoarseLevel(vRes, vWell, f.weights);

    for (std::size_t i = 0; i < coarseDim; ++i) {
        BOOST_CHECK_CLOSE(stage.coarseRhs()[i][0], 1.0 + 2.0 * i, 1e-10);
    }
}

// well_transfer only changes the vector transfers -- the coarse matrix is the
// same in every mode.
BOOST_AUTO_TEST_CASE(WellTransferModeDoesNotChangeCoarseMatrix)
{
    const Fixture f;
    const std::size_t coarseDim = numCells + f.layout.numWells();

    Stage full(f.S, Opm::PropertyTree(), pressureIndex, Opm::WellTransfer::Full);
    Stage classic(f.S, Opm::PropertyTree(), pressureIndex, Opm::WellTransfer::Classic);
    full.buildCoarseSystem(f.weights);
    classic.buildCoarseSystem(f.weights);

    for (std::size_t r = 0; r < coarseDim; ++r) {
        for (std::size_t c = 0; c < coarseDim; ++c) {
            BOOST_CHECK_EQUAL(coarseEntry(full, r, c), coarseEntry(classic, r, c));
        }
    }
}

// Classic mode must leave the well rows of the coarse right-hand side at zero
// and must not write any well correction, matching PressureBhpTransferPolicy.
BOOST_AUTO_TEST_CASE(ClassicTransferDropsWellResidualAndCorrection)
{
    const Fixture f;
    const std::size_t coarseDim = numCells + f.layout.numWells();

    Opm::ResVector<Scalar> dRes(numCells);
    for (std::size_t c = 0; c < numCells; ++c) {
        for (int i = 0; i < numRes; ++i) {
            dRes[c][i] = 1.0 + c + i;
        }
    }
    Opm::WellVector<Scalar> dWell(numWellBlocks);
    for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
        for (int i = 0; i < numWell; ++i) {
            dWell[wb][i] = 2.0 + wb - i;
        }
    }

    for (const auto mode : {Opm::WellTransfer::Full,
                            Opm::WellTransfer::NoProlongation,
                            Opm::WellTransfer::Classic}) {
        Stage stage(f.S, Opm::PropertyTree(), pressureIndex, mode);
        stage.buildCoarseSystem(f.weights);
        stage.moveToCoarseLevel(dRes, dWell, f.weights);

        const bool restricts = (mode != Opm::WellTransfer::Classic);
        for (std::size_t j = 0; j < f.layout.numWells(); ++j) {
            if (restricts) {
                BOOST_CHECK_NE(stage.coarseRhs()[numCells + j][0], 0.0);
            } else {
                BOOST_CHECK_EQUAL(stage.coarseRhs()[numCells + j][0], 0.0);
            }
        }
        // The reservoir rows are restricted identically in every mode.
        for (std::size_t c = 0; c < numCells; ++c) {
            Scalar expected = 0.0;
            for (int i = 0; i < numRes; ++i) {
                expected += dRes[c][i] * f.weights[Dune::Indices::_0][c][i];
            }
            BOOST_CHECK_CLOSE(stage.coarseRhs()[c][0], expected, 1e-10);
        }

        auto& coarseSol = stage.coarseSolution();
        coarseSol.resize(coarseDim);
        for (std::size_t i = 0; i < coarseDim; ++i) {
            coarseSol[i] = 1.0 + i;
        }
        Opm::ResVector<Scalar> vRes(numCells);
        Opm::WellVector<Scalar> vWell(numWellBlocks);
        stage.moveToFineLevel(vRes, vWell);

        const bool prolongs = (mode == Opm::WellTransfer::Full);
        BOOST_CHECK_EQUAL(stage.prolongatesWellPressure(), prolongs);
        for (std::size_t wb = 0; wb < numWellBlocks; ++wb) {
            for (int i = 0; i < numWell; ++i) {
                if (!prolongs) {
                    BOOST_CHECK_EQUAL(vWell[wb][i], 0.0);
                }
            }
        }
        if (prolongs) {
            for (std::size_t j = 0; j < f.layout.numWells(); ++j) {
                BOOST_CHECK_NE(vWell[f.layout.firstBlock(j)][wellPressureIndex], 0.0);
            }
        }
        // The reservoir prolongation is the same in every mode.
        for (std::size_t c = 0; c < numCells; ++c) {
            BOOST_CHECK_CLOSE(vRes[c][pressureIndex], 1.0 + c, 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(WellTransferFromStringRejectsUnknownValues)
{
    BOOST_CHECK(Opm::wellTransferFromString("full") == Opm::WellTransfer::Full);
    BOOST_CHECK(Opm::wellTransferFromString("no_prolongation") == Opm::WellTransfer::NoProlongation);
    BOOST_CHECK(Opm::wellTransferFromString("classic") == Opm::WellTransfer::Classic);
    BOOST_CHECK_THROW(Opm::wellTransferFromString("nonsense"), std::invalid_argument);
}

// A layout with no wells must reproduce a reservoir-only coarse system.
BOOST_AUTO_TEST_CASE(NoWellsGivesReservoirOnlyCoarseSystem)
{
    Fixture f;
    Opm::WellDofLayout emptyLayout;
    emptyLayout.wellBlockOffsets = {0};
    f.S.wellLayout = &emptyLayout;

    Stage stage(f.S, Opm::PropertyTree(), pressureIndex);
    stage.buildCoarseSystem(f.weights);

    BOOST_CHECK_EQUAL(stage.coarseMatrix().N(), numCells);
    BOOST_CHECK_EQUAL(stage.coarseMatrix().M(), numCells);
}
