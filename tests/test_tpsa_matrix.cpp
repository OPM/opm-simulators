/*
  Copyright 2025 NORCE AS

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

#define BOOST_TEST_MODULE TpsaMatrixTest

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/tpsa/TpsaMatrix.hpp>
#include <opm/simulators/linalg/tpsa/TpsaVector.hpp>

#include <dune/common/indices.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

//! \brief Scalar types every test case is instantiated for.
#if FLOW_INSTANTIATE_FLOAT
using ScalarTypes = boost::mpl::list<double, float>;
#else
using ScalarTypes = boost::mpl::list<double>;
#endif // FLOW_INSTANTIATE_FLOAT

template <class Scalar>
using Matrix = Opm::Linear::TpsaMatrix<Scalar>;
template <class Scalar>
using Vector = Opm::Linear::TpsaVector<Scalar>;
template <class Scalar>
using MatrixBlock = typename Matrix<Scalar>::MatrixBlock;
template <class Scalar>
using EqVector = typename Vector<Scalar>::EqVector;

constexpr int numEq = Opm::Linear::numTpsaEq;
constexpr int numFields = Opm::Linear::numTpsaFields;
constexpr std::size_t numCells = 5;

/*!
 * \brief Relative tolerance, in percent, of the inexact comparisons.
 *
 * The assembled values are exactly representable in both precisions, so only the results of
 * accumulation and scaling need a tolerance. float carries roughly seven significant digits,
 * which the double tolerance is far below.
 *
 * \tparam Scalar Field type the comparison is made in.
 */
template <class Scalar>
struct Tolerance;

template <>
struct Tolerance<double>
{
    static constexpr double percent = 1e-10;
};

template <>
struct Tolerance<float>
{
    static constexpr float percent = 1e-3;
};

/*!
 * \brief Tri-diagonal stencil, as produced by TpsaLinearizer for a 1-D grid.
 *
 * \return One entry per cell, holding that row's column indices in ascending order.
 */
std::vector<std::set<unsigned> >
makePattern()
{
    std::vector<std::set<unsigned> > pattern(numCells);
    for (unsigned i = 0; i < numCells; ++i) {
        pattern[i].insert(i);
        if (i > 0) {
            pattern[i].insert(i - 1);
        }
        if (i + 1 < numCells) {
            pattern[i].insert(i + 1);
        }
    }

    return pattern;
}

/*!
 * \brief A distinct value for every scalar of every block.
 *
 * \tparam Scalar Field type of the value.
 *
 * \param[in] row Block row index.
 * \param[in] col Block column index.
 * \param[in] eqIdx Equation index within the block.
 * \param[in] pvIdx Primary variable index within the block.
 *
 * \return The value the given scalar is expected to hold.
 */
template <class Scalar>
Scalar
entry(std::size_t row, std::size_t col, int eqIdx, int pvIdx)
{
    return Scalar(1000) * static_cast<Scalar>(row) + Scalar(100) * static_cast<Scalar>(col)
        + Scalar(10) * static_cast<Scalar>(eqIdx) + static_cast<Scalar>(pvIdx) + Scalar(1);
}

/*!
 * \brief Build the dense 7x7 block belonging to one (row, col) position.
 *
 * \tparam Scalar Field type of the block entries.
 *
 * \param[in] row Block row index.
 * \param[in] col Block column index.
 *
 * \return Dense block whose entries are filled from entry().
 */
template <class Scalar>
MatrixBlock<Scalar>
makeBlock(std::size_t row, std::size_t col)
{
    MatrixBlock<Scalar> b(Scalar(0));
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            b[eqIdx][pvIdx] = entry<Scalar>(row, col, eqIdx, pvIdx);
        }
    }

    return b;
}

/*!
 * \brief True for the six displacement-displacement couplings that are dropped.
 *
 * \param[in] eqIdx Equation index within the block.
 * \param[in] pvIdx Primary variable index within the block.
 *
 * \return Whether that entry is a displacement-displacement off-diagonal
 */
bool
isDropped(int eqIdx, int pvIdx)
{
    return eqIdx < 3 && pvIdx < 3 && eqIdx != pvIdx;
}

/*!
 * \brief Field a given equation belongs to.
 *
 * \param[in] idx Equation or primary variable index within a block.
 *
 * \return Field index: 0/1/2 for u_x/u_y/u_z, 3 for the rotation and 4 for the solid pressure.
 */
std::size_t
fieldOf(int idx)
{
    if (idx < 3) {
        return static_cast<std::size_t>(idx);
    }

    return (idx < 6) ? 3 : 4;
}

/*!
 * \brief Assemble the test matrix the way the linearizer would.
 *
 * Adds makeBlock(row, col) to every block of the pattern, going through the same
 * BlockAddress interface the linearizer uses.
 *
 * \tparam Scalar Field type of the matrix entries.
 *
 * \param[in,out] m Matrix the blocks are added to. Must already be reserved with \p pattern.
 * \param[in] pattern Sparsity pattern the matrix was reserved with.
 */
template <class Scalar>
void
assemble(Matrix<Scalar>& m, const std::vector<std::set<unsigned> >& pattern)
{
    for (std::size_t row = 0; row < numCells; ++row) {
        for (const auto col : pattern[row]) {
            auto address = m.blockAddress(row, col);
            const auto b = makeBlock<Scalar>(row, col);
            // Same syntax the linearizer uses on a MatrixBlock*.
            *address += b;
        }
    }
}

/*!
 * \brief A distinct value for every equation of every cell.
 *
 * \tparam Scalar Field type of the value.
 *
 * \param[in] dofIdx Index of the cell.
 * \param[in] eqIdx Equation index within the cell.
 *
 * \return The value the given entry is expected to hold.
 */
template <class Scalar>
Scalar
vectorEntry(std::size_t dofIdx, int eqIdx)
{
    return Scalar(0.5) * static_cast<Scalar>(dofIdx) + Scalar(0.25) * static_cast<Scalar>(eqIdx)
        + Scalar(1);
}

/*!
 * \brief Fill a vector with vectorEntry(), optionally alternating the sign.
 *
 * \tparam Scalar Field type of the vector entries.
 *
 * \param[out] v Vector that is filled. Must already be sized for numCells.
 * \param[in] alternateSign Whether to flip the sign of every other entry, so that the norms
 *                          and the plain sums differ.
 */
template <class Scalar>
void
fillVector(Vector<Scalar>& v, bool alternateSign = false)
{
    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar sign = (alternateSign && ((i + eqIdx) % 2 == 1)) ? Scalar(-1) : Scalar(1);
            v[i][eqIdx] = sign * vectorEntry<Scalar>(i, eqIdx);
        }
    }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(ScatterMatchesTheFieldSplit, Scalar, ScalarTypes)
{
    using namespace Dune::Indices;

    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);
    const auto& view = m.istlMatrix();

    for (std::size_t row = 0; row < numCells; ++row) {
        for (const auto col : pattern[row]) {
            // Diagonal displacement blocks: 1x1, one per component.
            BOOST_CHECK_EQUAL(view[_0][_0][row][col][0][0], entry<Scalar>(row, col, 0, 0));
            BOOST_CHECK_EQUAL(view[_1][_1][row][col][0][0], entry<Scalar>(row, col, 1, 1));
            BOOST_CHECK_EQUAL(view[_2][_2][row][col][0][0], entry<Scalar>(row, col, 2, 2));

            for (int j = 0; j < 3; ++j) {
                // Displacement-rotation: 1x3
                BOOST_CHECK_EQUAL(view[_0][_3][row][col][0][j], entry<Scalar>(row, col, 0, 3 + j));
                BOOST_CHECK_EQUAL(view[_1][_3][row][col][0][j], entry<Scalar>(row, col, 1, 3 + j));
                BOOST_CHECK_EQUAL(view[_2][_3][row][col][0][j], entry<Scalar>(row, col, 2, 3 + j));

                // Rotation-displacement: 3x1
                BOOST_CHECK_EQUAL(view[_3][_0][row][col][j][0], entry<Scalar>(row, col, 3 + j, 0));
                BOOST_CHECK_EQUAL(view[_3][_1][row][col][j][0], entry<Scalar>(row, col, 3 + j, 1));
                BOOST_CHECK_EQUAL(view[_3][_2][row][col][j][0], entry<Scalar>(row, col, 3 + j, 2));

                // Rotation-solid pressure and solid pressure-rotation
                BOOST_CHECK_EQUAL(view[_3][_4][row][col][j][0], entry<Scalar>(row, col, 3 + j, 6));
                BOOST_CHECK_EQUAL(view[_4][_3][row][col][0][j], entry<Scalar>(row, col, 6, 3 + j));

                for (int i = 0; i < 3; ++i) {
                    BOOST_CHECK_EQUAL(view[_3][_3][row][col][i][j],
                                      entry<Scalar>(row, col, 3 + i, 3 + j));
                }
            }

            // Displacement-solid pressure and solid pressure-displacement
            BOOST_CHECK_EQUAL(view[_0][_4][row][col][0][0], entry<Scalar>(row, col, 0, 6));
            BOOST_CHECK_EQUAL(view[_1][_4][row][col][0][0], entry<Scalar>(row, col, 1, 6));
            BOOST_CHECK_EQUAL(view[_2][_4][row][col][0][0], entry<Scalar>(row, col, 2, 6));
            BOOST_CHECK_EQUAL(view[_4][_0][row][col][0][0], entry<Scalar>(row, col, 6, 0));
            BOOST_CHECK_EQUAL(view[_4][_1][row][col][0][0], entry<Scalar>(row, col, 6, 1));
            BOOST_CHECK_EQUAL(view[_4][_2][row][col][0][0], entry<Scalar>(row, col, 6, 2));

            BOOST_CHECK_EQUAL(view[_4][_4][row][col][0][0], entry<Scalar>(row, col, 6, 6));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(GatherDropsDisplacementCrossCouplings, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    // Reading a block back gives the assembled values, except for the six
    // displacement-displacement off-diagonals, which are not stored at all, and should return 0.
    for (std::size_t row = 0; row < numCells; ++row) {
        for (const auto col : pattern[row]) {
            MatrixBlock<Scalar> b;
            m.block(row, col, b);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                    const Scalar expected = isDropped(eqIdx, pvIdx)
                        ? Scalar(0)
                        : entry<Scalar>(row, col, eqIdx, pvIdx);
                    BOOST_CHECK_EQUAL(b[eqIdx][pvIdx], expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MatrixVectorProductMatchesDenseReference, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    Vector<Scalar> x(numCells);
    fillVector(x);

    Vector<Scalar> y(numCells);
    y = Scalar(0);
    m.istlMatrix().mv(x.istlVector(), y.istlVector());

    // Reference: dense multiply with the dropped entries set to zero.
    for (std::size_t row = 0; row < numCells; ++row) {
        EqVector<Scalar> expected(Scalar(0));
        for (const auto col : pattern[row]) {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                    if (isDropped(eqIdx, pvIdx)) {
                        continue;
                    }
                    expected[eqIdx] += entry<Scalar>(row, col, eqIdx, pvIdx)
                        * vectorEntry<Scalar>(col, pvIdx);
                }
            }
        }

        const EqVector<Scalar> actual = static_cast<const Vector<Scalar>&>(y)[row];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            BOOST_CHECK_CLOSE(actual[eqIdx], expected[eqIdx], Tolerance<Scalar>::percent);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ClearAndClearRow, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    constexpr std::size_t row = 2;
    m.clearRow(row, Scalar(1));

    MatrixBlock<Scalar> b;
    for (const auto col : pattern[row]) {
        m.block(row, col, b);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                const bool onDiagonal = col == row && eqIdx == pvIdx;
                BOOST_CHECK_EQUAL(b[eqIdx][pvIdx], onDiagonal ? Scalar(1) : Scalar(0));
            }
        }
    }

    // Other rows are untouched.
    m.block(1, 1, b);
    BOOST_CHECK_EQUAL(b[0][0], entry<Scalar>(1, 1, 0, 0));

    m.clear();
    m.block(1, 1, b);
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            BOOST_CHECK_EQUAL(b[eqIdx][pvIdx], Scalar(0));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(HypreCompatibleDiagonalBlocks, Scalar, ScalarTypes)
{
    // Check that displacement and solid-pressure diagonal blocks are 1x1
    using DD00 = std::decay_t<decltype(std::declval<Matrix<Scalar>&>().dd00())>;
    using SPSP = std::decay_t<decltype(std::declval<Matrix<Scalar>&>().spsp())>;
    static_assert(DD00::block_type::rows == 1 && DD00::block_type::cols == 1);
    static_assert(SPSP::block_type::rows == 1 && SPSP::block_type::cols == 1);

    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    const Scalar* values = &m.dd00()[0][0][0][0];
    std::size_t k = 0;
    for (std::size_t row = 0; row < numCells; ++row) {
        for (const auto col : pattern[row]) {
            BOOST_CHECK_EQUAL(values[k], entry<Scalar>(row, col, 0, 0));
            ++k;
        }
    }
    BOOST_CHECK_EQUAL(k, m.nonzeroes());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SetAndAddBlock, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);

    // (2, 3) is on the tri-diagonal, so it is part of the pattern.
    const auto b = makeBlock<Scalar>(2, 3);
    m.setBlock(2, 3, b);

    MatrixBlock<Scalar> read;
    m.block(2, 3, read);
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            const Scalar expected = isDropped(eqIdx, pvIdx)
                ? Scalar(0)
                : entry<Scalar>(2, 3, eqIdx, pvIdx);
            BOOST_CHECK_EQUAL(read[eqIdx][pvIdx], expected);
        }
    }

    // addToBlock accumulates on top of what is already stored.
    m.addToBlock(2, 3, b);
    m.block(2, 3, read);
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            const Scalar expected = isDropped(eqIdx, pvIdx)
                ? Scalar(0)
                : Scalar(2) * entry<Scalar>(2, 3, eqIdx, pvIdx);
            BOOST_CHECK_EQUAL(read[eqIdx][pvIdx], expected);
        }
    }

    // setBlock overwrites rather than accumulating.
    m.setBlock(2, 3, b);
    m.block(2, 3, read);
    BOOST_CHECK_EQUAL(read[0][0], entry<Scalar>(2, 3, 0, 0));

    // Assigning a scalar through the block handle sets every stored entry.
    auto address = m.blockAddress(2, 3);
    *address = Scalar(7);
    m.block(2, 3, read);
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            const Scalar expected = isDropped(eqIdx, pvIdx) ? Scalar(0) : Scalar(7);
            BOOST_CHECK_EQUAL(read[eqIdx][pvIdx], expected);
        }
    }

    // Neighbouring blocks are untouched.
    m.block(2, 2, read);
    BOOST_CHECK_EQUAL(read[0][0], Scalar(0));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MatrixScaleFields, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    // Distinct factors, so that a sub-matrix picking up the wrong field shows up.
    const std::array<Scalar, numFields> rowFac{
        {Scalar(2), Scalar(3), Scalar(5), Scalar(7), Scalar(11)}};
    const std::array<Scalar, numFields> colFac{
        {Scalar(13), Scalar(17), Scalar(19), Scalar(23), Scalar(29)}};
    m.scaleFields(rowFac, colFac);

    const auto checkScaled = [&]() {
        for (std::size_t row = 0; row < numCells; ++row) {
            for (const auto col : pattern[row]) {
                MatrixBlock<Scalar> b;
                m.block(row, col, b);
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                        if (isDropped(eqIdx, pvIdx)) {
                            BOOST_CHECK_EQUAL(b[eqIdx][pvIdx], Scalar(0));
                            continue;
                        }

                        const Scalar expected = entry<Scalar>(row, col, eqIdx, pvIdx)
                            * rowFac[fieldOf(eqIdx)] * colFac[fieldOf(pvIdx)];
                        BOOST_CHECK_CLOSE(b[eqIdx][pvIdx], expected, Tolerance<Scalar>::percent);
                    }
                }
            }
        }
    };
    checkScaled();

    // Unit factors take the early-out branch and must leave the matrix alone.
    const std::array<Scalar, numFields> ones{
        {Scalar(1), Scalar(1), Scalar(1), Scalar(1), Scalar(1)}};
    m.scaleFields(ones, ones);
    checkScaled();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(UmvAndUsmvAccumulate, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    Vector<Scalar> x(numCells);
    fillVector(x);

    // mv is checked against a dense reference elsewhere, so it can serve as the
    // reference for the accumulating variants.
    Vector<Scalar> product(numCells);
    product = Scalar(0);
    m.istlMatrix().mv(x.istlVector(), product.istlVector());

    Vector<Scalar> initial(numCells);
    fillVector(initial, /*alternateSign =*/ true);

    Vector<Scalar> y = initial;
    m.istlMatrix().umv(x.istlVector(), y.istlVector());
    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            BOOST_CHECK_CLOSE(y[i][eqIdx],
                              initial[i][eqIdx] + product[i][eqIdx],
                              Tolerance<Scalar>::percent);
        }
    }

    const Scalar alpha = Scalar(-2.5);
    Vector<Scalar> z = initial;
    m.istlMatrix().usmv(alpha, x.istlVector(), z.istlVector());
    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            BOOST_CHECK_CLOSE(z[i][eqIdx],
                              initial[i][eqIdx] + alpha * product[i][eqIdx],
                              Tolerance<Scalar>::percent);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MakeOverlapRowsInvalid, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);
    assemble(m, pattern);

    const std::vector<int> overlapRows{1, 3};
    m.makeOverlapRowsInvalid(overlapRows);

    for (std::size_t row = 0; row < numCells; ++row) {
        const bool isOverlap = (row == 1) || (row == 3);
        for (const auto col : pattern[row]) {
            MatrixBlock<Scalar> b;
            m.block(row, col, b);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                    Scalar expected = Scalar(0);
                    if (isDropped(eqIdx, pvIdx)) {
                        expected = Scalar(0);
                    }
                    else if (!isOverlap) {
                        expected = entry<Scalar>(row, col, eqIdx, pvIdx);
                    }
                    else if (col == row && eqIdx == pvIdx) {
                        expected = Scalar(1);
                    }

                    BOOST_CHECK_EQUAL(b[eqIdx][pvIdx], expected);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ThrowsOnInvalidUse, Scalar, ScalarTypes)
{
    const auto pattern = makePattern();
    Matrix<Scalar> m(numCells, numCells);
    m.reserve(pattern);

    // Row 0 of the tri-diagonal stencil only couples to columns 0 and 1.
    MatrixBlock<Scalar> b;
    BOOST_CHECK_THROW(m.blockAddress(0, 3), std::logic_error);
    BOOST_CHECK_THROW(m.block(0, 3, b), std::logic_error);
    BOOST_CHECK_THROW(m.setBlock(0, 3, b), std::logic_error);

    // A pattern that does not have one entry per matrix row.
    Matrix<Scalar> mismatched(numCells, numCells);
    const std::vector<std::set<unsigned> > shortPattern(numCells - 1);
    BOOST_CHECK_THROW(mismatched.reserve(shortPattern), std::logic_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(VectorRoundTrip, Scalar, ScalarTypes)
{
    using namespace Dune::Indices;

    Vector<Scalar> v(numCells);
    v = Scalar(0);

    for (std::size_t i = 0; i < numCells; ++i) {
        EqVector<Scalar> contribution;
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            contribution[eqIdx] = Scalar(10) * static_cast<Scalar>(i) + static_cast<Scalar>(eqIdx);
        }
        v[i] += contribution;
        v[i] += contribution;
    }

    const Vector<Scalar>& constView = v;
    for (std::size_t i = 0; i < numCells; ++i) {
        const EqVector<Scalar> block = constView[i];
        BOOST_CHECK_EQUAL(block.size(), static_cast<std::size_t>(numEq));
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar expected
                = Scalar(2) * (Scalar(10) * static_cast<Scalar>(i) + static_cast<Scalar>(eqIdx));
            BOOST_CHECK_CLOSE(block[eqIdx], expected, Tolerance<Scalar>::percent);
        }
    }

    // Sub-vectors carry the fields in the expected order.
    BOOST_CHECK_CLOSE(v.istlVector()[_0][3][0], Scalar(2) * Scalar(30), Tolerance<Scalar>::percent);
    BOOST_CHECK_CLOSE(v.istlVector()[_3][3][1],
                      Scalar(2) * Scalar(30 + 4),
                      Tolerance<Scalar>::percent);
    BOOST_CHECK_CLOSE(v.istlVector()[_4][3][0],
                      Scalar(2) * Scalar(30 + 6),
                      Tolerance<Scalar>::percent);

    Scalar expectedOneNorm = Scalar(0);
    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            expectedOneNorm
                += Scalar(2) * (Scalar(10) * static_cast<Scalar>(i) + static_cast<Scalar>(eqIdx));
        }
    }
    BOOST_CHECK_CLOSE(v.one_norm(), expectedOneNorm, Tolerance<Scalar>::percent);

    v[2] = Scalar(0);
    BOOST_CHECK_EQUAL(constView[2][4], Scalar(0));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(VectorArithmetic, Scalar, ScalarTypes)
{
    Vector<Scalar> a(numCells);
    fillVector(a);

    Vector<Scalar> b(numCells);
    fillVector(b, /*alternateSign =*/ true);

    Vector<Scalar> sum = a;
    sum += b;

    Vector<Scalar> difference = a;
    difference -= b;

    Vector<Scalar> scaled = a;
    scaled *= Scalar(3);

    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar av = a[i][eqIdx];
            const Scalar bv = b[i][eqIdx];
            BOOST_CHECK_CLOSE(sum[i][eqIdx], av + bv, Tolerance<Scalar>::percent);
            BOOST_CHECK_CLOSE(difference[i][eqIdx], av - bv, Tolerance<Scalar>::percent);
            BOOST_CHECK_CLOSE(scaled[i][eqIdx], Scalar(3) * av, Tolerance<Scalar>::percent);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(VectorNorms, Scalar, ScalarTypes)
{
    Vector<Scalar> v(numCells);
    fillVector(v, /*alternateSign =*/ true);

    Scalar expectedOneNorm = Scalar(0);
    Scalar expectedTwoNorm2 = Scalar(0);
    Scalar expectedInfinityNorm = Scalar(0);
    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar value = std::abs(vectorEntry<Scalar>(i, eqIdx));
            expectedOneNorm += value;
            expectedTwoNorm2 += value * value;
            expectedInfinityNorm = std::max(expectedInfinityNorm, value);
        }
    }

    BOOST_CHECK_CLOSE(v.one_norm(), expectedOneNorm, Tolerance<Scalar>::percent);
    BOOST_CHECK_CLOSE(v.two_norm2(), expectedTwoNorm2, Tolerance<Scalar>::percent);
    BOOST_CHECK_CLOSE(v.two_norm(), std::sqrt(expectedTwoNorm2), Tolerance<Scalar>::percent);
    BOOST_CHECK_CLOSE(v.infinity_norm(), expectedInfinityNorm, Tolerance<Scalar>::percent);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(VectorScaleFields, Scalar, ScalarTypes)
{
    Vector<Scalar> v(numCells);
    fillVector(v);

    const std::array<Scalar, numFields> factors{
        {Scalar(2), Scalar(3), Scalar(5), Scalar(7), Scalar(11)}};
    v.scaleFields(factors);

    for (std::size_t i = 0; i < numCells; ++i) {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar expected = vectorEntry<Scalar>(i, eqIdx) * factors[fieldOf(eqIdx)];
            BOOST_CHECK_CLOSE(v[i][eqIdx], expected, Tolerance<Scalar>::percent);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(EntryProxyOperations, Scalar, ScalarTypes)
{
    BOOST_CHECK_EQUAL(Vector<Scalar>::EntryProxy::size(), static_cast<std::size_t>(numEq));

    Vector<Scalar> v(numCells);
    v = Scalar(0);

    EqVector<Scalar> value;
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        value[eqIdx] = vectorEntry<Scalar>(1, eqIdx);
    }

    // Scatter a dense block, then gather it back through the conversion operator.
    v[1] = value;
    const EqVector<Scalar> gathered = v[1];
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_CLOSE(gathered[eqIdx], value[eqIdx], Tolerance<Scalar>::percent);
    }

    v[1] *= Scalar(4);
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_CLOSE(v[1][eqIdx], Scalar(4) * value[eqIdx], Tolerance<Scalar>::percent);
    }

    v[1] -= value;
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_CLOSE(v[1][eqIdx], Scalar(3) * value[eqIdx], Tolerance<Scalar>::percent);
    }

    // Only cell 1 was touched.
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_EQUAL(v[0][eqIdx], Scalar(0));
        BOOST_CHECK_EQUAL(v[2][eqIdx], Scalar(0));
    }

    // Proxy-to-proxy assignment
    v[3] = v[1];
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        BOOST_CHECK_CLOSE(v[3][eqIdx], v[1][eqIdx], Tolerance<Scalar>::percent);
    }
}
