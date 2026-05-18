#include <config.h>

#define BOOST_TEST_MODULE OPM_test_WellMatrixMerger
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace {

using Scalar = double;
using WRMatrix = Opm::WRMatrixT<Scalar>;
using RWMatrix = Opm::RWMatrixT<Scalar>;
using WWMatrix = Opm::WWMatrixT<Scalar>;
using WRBlock = WRMatrix::block_type;
using RWBlock = RWMatrix::block_type;
using WWBlock = WWMatrix::block_type;

// WellMatrixMerger assembles the global coupled well part of
//
//     [ A  C ]
//     [ B  D ]
//
// from the per-well blocks B_j, C_j and D_j. It preserves each well's local
// sparsity pattern and only does two structural operations: concatenate the
// well blocks and remap perforation-related rows/columns through the list of
// perforated reservoir cells for each well.

// Give each block a distinctive value pattern so it is easy to see where it
// ended up after merging.
template<class Block>
Block makeBlock(const Scalar base)
{
    Block block;
    for (int row = 0; row < Block::rows; ++row) {
        for (int col = 0; col < Block::cols; ++col) {
            block[row][col] = base + row * 10.0 + col;
        }
    }
    return block;
}

struct BlockSpec
{
    int column;
    Scalar base;
};

BlockSpec entry(const int column, const Scalar base)
{
    return {column, base};
}

using MatrixRow = std::vector<BlockSpec>;
using MatrixPattern = std::vector<std::vector<BlockSpec>>;

// Small helper for hand-written BCRS matrices used throughout the test.
template<class Matrix>
Matrix buildMatrix(const std::size_t rows,
                   const std::size_t cols,
                   const MatrixPattern& pattern)
{
    std::size_t nonzeroes = 0;
    for (const auto& row : pattern) {
        nonzeroes += row.size();
    }

    Matrix matrix(rows, cols, nonzeroes, Matrix::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (const auto& entry : pattern[row.index()]) {
            row.insert(entry.column);
        }
    }

    for (std::size_t row = 0; row < rows; ++row) {
        for (const auto& entry : pattern[row]) {
            matrix[row][entry.column] = makeBlock<typename Matrix::block_type>(entry.base);
        }
    }

    return matrix;
}

template<class Block>
void checkBlockEqual(const Block& actual, const Block& expected)
{
    for (int row = 0; row < Block::rows; ++row) {
        for (int col = 0; col < Block::cols; ++col) {
            BOOST_CHECK_EQUAL(actual[row][col], expected[row][col]);
        }
    }
}

template<class Matrix>
void checkMatrixEqual(const Matrix& actual, const Matrix& expected)
{
    BOOST_REQUIRE_EQUAL(actual.N(), expected.N());
    BOOST_REQUIRE_EQUAL(actual.M(), expected.M());

    for (std::size_t row = 0; row < actual.N(); ++row) {
        auto actualIt = actual[row].begin();
        auto expectedIt = expected[row].begin();
        for (; actualIt != actual[row].end() && expectedIt != expected[row].end(); ++actualIt, ++expectedIt) {
            BOOST_CHECK_EQUAL(actualIt.index(), expectedIt.index());
            checkBlockEqual(*actualIt, *expectedIt);
        }
        BOOST_CHECK(actualIt == actual[row].end());
        BOOST_CHECK(expectedIt == expected[row].end());
    }
}

struct TestMatrices
{
    std::vector<WRMatrix> bMatrices;
    std::vector<RWMatrix> cMatrices;
    std::vector<WWMatrix> dMatrices;
    std::vector<std::vector<int>> wellCells;
};

struct MergedMatrices
{
    WRMatrix b;
    RWMatrix c;
    WWMatrix d;
};

// Two wells:
// - well 0 is a toy multisegment well with two segment blocks and two
//   perforated reservoir cells, 1 and 3. Each segment sees exactly one
//   perforation, and D contains the segment-to-segment coupling.
// - well 1 is a toy standard well with one well block and two perforated
//   reservoir cells, 0 and 4. Its single well block couples to both
//   perforations, so B is dense across the two perforation columns and C is
//   dense down the two perforation rows.
// The B/C/D matrices below are the local per-well matrices before merging.
TestMatrices buildTestMatrices()
{
    TestMatrices matrices;

    // wellCells[well][local perforation] = global reservoir cell index.
    // The first well perforates cells 1 and 3, the second well perforates
    // cells 0 and 4.
    matrices.wellCells = {{1, 3}, {0, 4}};

    // Toy MSW: two segment rows/columns and one local perforation attached to
    // each segment. B and C therefore have one entry per row. D has the full
    // 2x2 segment-coupling pattern, which is the part that matters for the
    // structure-cache test later.
    const MatrixPattern multisegmentB{{
        MatrixRow{entry(0, 10.0)},
        MatrixRow{entry(1, 20.0)}
    }};
    const MatrixPattern multisegmentC{{
        MatrixRow{entry(0, 30.0)},
        MatrixRow{entry(1, 40.0)}
    }};
    const MatrixPattern multisegmentD{{
        MatrixRow{entry(0, 50.0), entry(1, 60.0)},
        MatrixRow{entry(0, 70.0), entry(1, 80.0)}
    }};

    // Toy standard well: one well row/column block coupled to both local
    // perforations. B is therefore dense across its two perforation columns,
    // C has one entry in each perforated-cell row, and D is just 1x1.
    const MatrixPattern standardWellB{{
        MatrixRow{entry(0, 90.0), entry(1, 100.0)}
    }};
    const MatrixPattern standardWellC{{
        MatrixRow{entry(0, 110.0)},
        MatrixRow{entry(0, 120.0)}
    }};
    const MatrixPattern standardWellD{{
        MatrixRow{entry(0, 130.0)}
    }};

    // The vectors are ordered as [multisegment well, standard well].
    matrices.bMatrices.emplace_back(buildMatrix<WRMatrix>(
        2,
        2,
        multisegmentB));
    matrices.bMatrices.emplace_back(buildMatrix<WRMatrix>(
        1,
        2,
        standardWellB));

    matrices.cMatrices.emplace_back(buildMatrix<RWMatrix>(
        2,
        2,
        multisegmentC));
    matrices.cMatrices.emplace_back(buildMatrix<RWMatrix>(
        2,
        1,
        standardWellC));

    matrices.dMatrices.emplace_back(buildMatrix<WWMatrix>(
        2,
        2,
        multisegmentD));
    matrices.dMatrices.emplace_back(buildMatrix<WWMatrix>(
        1,
        1,
        standardWellD));

    return matrices;
}

WRMatrix expectedMergedB()
{
    // After merging, the well rows are concatenated in well order. The local
    // perforation columns are replaced by the corresponding global perforated
    // reservoir cells from wellCells: {0,1} -> {1,3} for the MSW and
    // {0,1} -> {0,4} for the standard well.
    return buildMatrix<WRMatrix>(
        3,
        5,
        MatrixPattern{{
            MatrixRow{entry(1, 10.0)},
            MatrixRow{entry(3, 20.0)},
            MatrixRow{entry(0, 90.0), entry(4, 100.0)}
        }});
}

RWMatrix expectedMergedC()
{
    // C writes to global reservoir rows. Only the perforated reservoir cells
    // get entries, and the well columns are appended in well order.
    return buildMatrix<RWMatrix>(
        5,
        3,
        MatrixPattern{{
            MatrixRow{entry(2, 110.0)},
            MatrixRow{entry(0, 30.0)},
            {},
            MatrixRow{entry(1, 40.0)},
            MatrixRow{entry(2, 120.0)}
        }});
}

WWMatrix expectedMergedD()
{
    // D is the block-diagonal concatenation of the per-well well/well blocks.
    return buildMatrix<WWMatrix>(
        3,
        3,
        MatrixPattern{{
            MatrixRow{entry(0, 50.0), entry(1, 60.0)},
            MatrixRow{entry(0, 70.0), entry(1, 80.0)},
            MatrixRow{entry(2, 130.0)}
        }});
}

// Change one existing entry in each local matrix so updateValues() can be
// checked without changing the sparsity pattern.
void updateSourceValues(TestMatrices& matrices)
{
    matrices.bMatrices[0][0][0] = makeBlock<WRBlock>(210.0);
    matrices.cMatrices[1][1][0] = makeBlock<RWBlock>(220.0);
    matrices.dMatrices[0][1][1] = makeBlock<WWBlock>(230.0);
}

MergedMatrices buildMergedMatrices(const TestMatrices& matrices,
                                   const std::size_t numResDof)
{
    MergedMatrices merged;
    const Opm::WellMatrixMerger<Scalar> merger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    merger.buildMatrices(merged.b, merged.c, merged.d);
    return merged;
}

void checkMergedMatrices(const MergedMatrices& merged)
{
    checkMatrixEqual(merged.b, expectedMergedB());
    checkMatrixEqual(merged.c, expectedMergedC());
    checkMatrixEqual(merged.d, expectedMergedD());
}

void checkUpdatedMergedMatrices(const MergedMatrices& merged)
{
    auto expectedB = expectedMergedB();
    auto expectedC = expectedMergedC();
    auto expectedD = expectedMergedD();
    expectedB[0][1] = makeBlock<WRBlock>(210.0);
    expectedC[4][2] = makeBlock<RWBlock>(220.0);
    expectedD[1][1] = makeBlock<WWBlock>(230.0);

    checkMatrixEqual(merged.b, expectedB);
    checkMatrixEqual(merged.c, expectedC);
    checkMatrixEqual(merged.d, expectedD);
}

void checkEmptyMergedMatrices(const MergedMatrices& merged,
                              const std::size_t numResDof)
{
    BOOST_CHECK_EQUAL(merged.b.N(), std::size_t{0});
    BOOST_CHECK_EQUAL(merged.b.M(), numResDof);
    BOOST_CHECK_EQUAL(merged.c.N(), numResDof);
    BOOST_CHECK_EQUAL(merged.c.M(), std::size_t{0});
    BOOST_CHECK_EQUAL(merged.d.N(), std::size_t{0});
    BOOST_CHECK_EQUAL(merged.d.M(), std::size_t{0});
}

} // namespace

BOOST_AUTO_TEST_CASE(MergeHandlesEmptyWellSet)
{
    constexpr std::size_t numResDof = 5;
    const TestMatrices matrices;
    const Opm::WellMatrixMerger<Scalar> merger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const auto structure = merger.buildStructure();

    BOOST_CHECK_EQUAL(structure.numResDofs, numResDof);
    BOOST_CHECK_EQUAL(structure.totalWellDofs, std::size_t{0});

    const auto merged = buildMergedMatrices(matrices, numResDof);
    checkEmptyMergedMatrices(merged, numResDof);
}

BOOST_AUTO_TEST_CASE(MergeBuildsExpectedMatricesAndStructure)
{
    constexpr std::size_t numResDof = 5;

    const auto matrices = buildTestMatrices();
    const Opm::WellMatrixMerger<Scalar> merger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const auto structure = merger.buildStructure();

    // The cached structure records both the perforated-cell mapping and the
    // exact sparsity of every per-well B, C and D block.
    BOOST_CHECK_EQUAL(structure.totalWellDofs, std::size_t{3});
    BOOST_CHECK(structure.wellCells == matrices.wellCells);
    BOOST_REQUIRE_EQUAL(structure.bPatterns.size(), matrices.bMatrices.size());
    BOOST_REQUIRE_EQUAL(structure.cPatterns.size(), matrices.cMatrices.size());
    BOOST_REQUIRE_EQUAL(structure.dPatterns.size(), matrices.dMatrices.size());
    BOOST_CHECK(structure.bPatterns[0] == Opm::captureMatrixSparsity(matrices.bMatrices[0]));
    BOOST_CHECK(structure.cPatterns[1] == Opm::captureMatrixSparsity(matrices.cMatrices[1]));
    BOOST_CHECK(structure.dPatterns[0] == Opm::captureMatrixSparsity(matrices.dMatrices[0]));
    BOOST_CHECK(merger.hasSameStructure(structure));

    const auto merged = buildMergedMatrices(matrices, numResDof);
    checkMergedMatrices(merged);
}

BOOST_AUTO_TEST_CASE(UpdateValuesReusesMergedStructure)
{
    constexpr std::size_t numResDof = 5;

    auto matrices = buildTestMatrices();
    const Opm::WellMatrixMerger<Scalar> merger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    auto merged = buildMergedMatrices(matrices, numResDof);

    const auto mergedBPattern = Opm::captureMatrixSparsity(merged.b);
    const auto mergedCPattern = Opm::captureMatrixSparsity(merged.c);
    const auto mergedDPattern = Opm::captureMatrixSparsity(merged.d);

    updateSourceValues(matrices);
    const Opm::WellMatrixMerger<Scalar> updatedMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    updatedMerger.updateValues(merged.b, merged.c, merged.d);

    BOOST_CHECK(mergedBPattern == Opm::captureMatrixSparsity(merged.b));
    BOOST_CHECK(mergedCPattern == Opm::captureMatrixSparsity(merged.c));
    BOOST_CHECK(mergedDPattern == Opm::captureMatrixSparsity(merged.d));
    BOOST_CHECK(merger.hasSameStructure(updatedMerger.buildStructure()));
    checkUpdatedMergedMatrices(merged);
}

BOOST_AUTO_TEST_CASE(StructureChangesWhenWellPatternChanges)
{
    constexpr std::size_t numResDof = 5;

    const auto matrices = buildTestMatrices();
    auto changedDMatrices = matrices.dMatrices;
    changedDMatrices[0] = buildMatrix<WWMatrix>(
        2,
        2,
        MatrixPattern{{
            MatrixRow{entry(0, 50.0)},
            MatrixRow{entry(1, 80.0)}
        }});

    // Mimic an MSW topology change: same number of segments, but the segment
    // connections disappear from D. The cached well structure must change
    // even though the matrix dimensions stay the same.
    const Opm::WellMatrixMerger<Scalar> referenceMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const Opm::WellMatrixMerger<Scalar> changedMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        changedDMatrices,
        matrices.wellCells);
    const auto reference = referenceMerger.buildStructure();
    const auto changed = changedMerger.buildStructure();

    BOOST_CHECK(referenceMerger.hasSameStructure(reference));
    BOOST_CHECK(!changedMerger.hasSameStructure(reference));
    BOOST_CHECK_EQUAL(reference.totalWellDofs, changed.totalWellDofs);
    BOOST_CHECK(reference != changed);
    BOOST_CHECK(reference.dPatterns[0] != changed.dPatterns[0]);
}

BOOST_AUTO_TEST_CASE(StructureChangesWhenCouplingPatternChanges)
{
    constexpr std::size_t numResDof = 5;

    const auto matrices = buildTestMatrices();
    auto changedBMatrices = matrices.bMatrices;
    // Remove one of the standard well's B couplings while keeping the block
    // dimensions unchanged. That should only change the cached B pattern for
    // well 1, leaving the corresponding C and D patterns untouched.
    changedBMatrices[1] = buildMatrix<WRMatrix>(
        1,
        2,
        MatrixPattern{{
            MatrixRow{entry(0, 90.0)}
        }});

    const Opm::WellMatrixMerger<Scalar> referenceMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const Opm::WellMatrixMerger<Scalar> changedMerger(
        numResDof,
        changedBMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const auto reference = referenceMerger.buildStructure();
    const auto changed = changedMerger.buildStructure();

    BOOST_CHECK(!changedMerger.hasSameStructure(reference));
    BOOST_CHECK(reference.bPatterns[1] != changed.bPatterns[1]);
    BOOST_CHECK(reference.cPatterns[1] == changed.cPatterns[1]);
    BOOST_CHECK(reference.dPatterns[1] == changed.dPatterns[1]);
}

BOOST_AUTO_TEST_CASE(StructureChangesWhenPerforationMappingChanges)
{
    constexpr std::size_t numResDof = 5;

    const auto matrices = buildTestMatrices();
    auto changedWellCells = matrices.wellCells;
    // Keep the local B/C/D sparsity unchanged but swap which global reservoir
    // cells the standard well perforates. The structure key must change
    // because the merged B and C entries move to different reservoir slots.
    changedWellCells[1] = {4, 0};

    const Opm::WellMatrixMerger<Scalar> referenceMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        matrices.wellCells);
    const Opm::WellMatrixMerger<Scalar> changedMerger(
        numResDof,
        matrices.bMatrices,
        matrices.cMatrices,
        matrices.dMatrices,
        changedWellCells);
    const auto reference = referenceMerger.buildStructure();
    const auto changed = changedMerger.buildStructure();

    BOOST_CHECK(!changedMerger.hasSameStructure(reference));
    BOOST_CHECK(reference != changed);
    BOOST_CHECK(reference.wellCells != changed.wellCells);
}