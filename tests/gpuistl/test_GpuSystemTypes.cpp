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

#define BOOST_TEST_MODULE TestGpuSystemTypes

#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/gpuistl/GpuSystemTypes.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>
#include <opm/simulators/linalg/system/WellMatrixMerger.hpp>

#include <dune/common/indices.hh>

#include <cstddef>
#include <vector>

using namespace Opm;
using namespace Opm::gpuistl;
using Dune::Indices::_0;
using Dune::Indices::_1;

namespace
{

// Fill a BCRS matrix with the given block pattern and deterministic values.
template <class Matrix>
Matrix
makeMatrix(std::size_t rows,
           std::size_t cols,
           const std::vector<std::vector<std::size_t>>& pattern,
           double offset)
{
    Matrix m(rows, cols, Matrix::row_wise);
    for (auto row = m.createbegin(); row != m.createend(); ++row) {
        for (const auto c : pattern[row.index()]) {
            row.insert(c);
        }
    }
    double value = offset;
    for (auto row = m.begin(); row != m.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (std::size_t i = 0; i < col->N(); ++i) {
                for (std::size_t j = 0; j < col->M(); ++j) {
                    (*col)[i][j] = value;
                    value += 0.25;
                }
            }
        }
    }
    return m;
}

template <class BlockVector>
void
fill(BlockVector& v, double offset)
{
    double value = offset;
    for (auto& block : v) {
        for (auto& x : block) {
            x = value;
            value -= 0.5;
        }
    }
}

template <class BlockVector>
void
checkClose(const BlockVector& expected, const gpuistl::GpuVector<double>& actual)
{
    BlockVector result(expected.size());
    actual.copyToHost(result);
    for (std::size_t i = 0; i < expected.size(); ++i) {
        for (std::size_t j = 0; j < expected[i].size(); ++j) {
            BOOST_CHECK_CLOSE(result[i][j], expected[i][j], 1e-10);
        }
    }
}

// Two reservoir cells, one well with a single block, perforating both cells.
struct SmallSystem {
    RRMatrix<double> A = makeMatrix<RRMatrix<double>>(2, 2, {{0, 1}, {0, 1}}, 10.0);
    RWMatrix<double> C = makeMatrix<RWMatrix<double>>(2, 1, {{0}, {0}}, -1.0);
    WRMatrix<double> B = makeMatrix<WRMatrix<double>>(1, 2, {{0, 1}}, 2.0);
    WWMatrix<double> D = makeMatrix<WWMatrix<double>>(1, 1, {{0}}, 5.0);
};

} // namespace

BOOST_AUTO_TEST_CASE(NonSquareScalarMatrix)
{
    // Three block columns, of which the middle one is empty.
    const auto C = makeMatrix<RWMatrix<double>>(2, 3, {{0, 2}, {2}}, 1.0);
    auto gpuC = makeScalarGpuMatrix(C);
    BOOST_CHECK_EQUAL(gpuC.N(), 2 * numResDofs);
    BOOST_CHECK_EQUAL(gpuC.M(), 3 * numWellDofs);

    WellVector<double> x(3);
    fill(x, 1.0);
    ResVector<double> y(2);
    fill(y, -3.0);

    gpuistl::GpuVector<double> gpuX(x);
    gpuistl::GpuVector<double> gpuY(y);
    gpuC.usmv(-2.0, gpuX, gpuY);
    C.usmv(-2.0, x, y);
    checkClose(y, gpuY);
}

BOOST_AUTO_TEST_CASE(NonSquareScalarMatrixUpdate)
{
    auto B = makeMatrix<WRMatrix<double>>(2, 3, {{0, 1}, {1, 2}}, 0.5);
    auto gpuB = makeScalarGpuMatrix(B);

    // New values, same pattern.
    for (auto row = B.begin(); row != B.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            *col *= -3.0;
        }
    }
    std::vector<double> scratch;
    updateScalarGpuMatrix(B, gpuB, scratch);

    ResVector<double> x(3);
    fill(x, 2.0);
    WellVector<double> y(2);
    gpuistl::GpuVector<double> gpuX(x);
    gpuistl::GpuVector<double> gpuY(2 * numWellDofs);
    gpuB.mv(gpuX, gpuY);
    B.mv(x, y);
    checkClose(y, gpuY);
}

BOOST_AUTO_TEST_CASE(SystemMatrixApply)
{
    SmallSystem s;

    SystemMatrix<double> cpuS;
    cpuS.A = &s.A;
    cpuS.B = &s.B;
    cpuS.C = &s.C;
    cpuS.D = &s.D;

    const auto gpuA = GpuRRMatrix<double>::fromMatrix(s.A);
    const auto gpuB = makeScalarGpuMatrix(s.B);
    const auto gpuC = makeScalarGpuMatrix(s.C);
    const auto gpuD = makeScalarGpuMatrix(s.D);
    GpuSystemMatrix<double> gpuS;
    gpuS.A = &gpuA;
    gpuS.B = &gpuB;
    gpuS.C = &gpuC;
    gpuS.D = &gpuD;
    *gpuS.hostD = &s.D;

    SystemVector<double> x;
    x[_0].resize(2);
    x[_1].resize(1);
    fill(x[_0], 1.0);
    fill(x[_1], 4.0);
    SystemVector<double> y = x;
    fill(y[_0], 0.3);
    fill(y[_1], -0.7);

    GpuSystemVector<double> gpuX(2 * numResDofs, numWellDofs);
    gpuX.res.copyFromHost(x[_0]);
    gpuX.well.copyFromHost(x[_1]);
    GpuSystemVector<double> gpuY(2 * numResDofs, numWellDofs);
    gpuY.res.copyFromHost(y[_0]);
    gpuY.well.copyFromHost(y[_1]);

    gpuS.usmv(0.5, gpuX, gpuY);
    cpuS.usmv(0.5, x, y);
    checkClose(y[_0], gpuY.res);
    checkClose(y[_1], gpuY.well);

    gpuS.mv(gpuX, gpuY);
    cpuS.mv(x, y);
    checkClose(y[_0], gpuY.res);
    checkClose(y[_1], gpuY.well);

    BOOST_CHECK_CLOSE(gpuX.two_norm(), x.two_norm(), 1e-10);
    BOOST_CHECK_CLOSE(gpuX.dot(gpuY), x.dot(y), 1e-10);
}

BOOST_AUTO_TEST_CASE(SystemVectorWithoutWells)
{
    GpuSystemVector<double> v(2 * numResDofs, 0);
    v = 1.0;
    BOOST_CHECK_EQUAL(v.dim(), 2 * numResDofs);
    BOOST_CHECK_CLOSE(v.two_norm() * v.two_norm(), 2.0 * numResDofs, 1e-10);

    GpuSystemVector<double> w(2 * numResDofs, 0);
    w = 2.0;
    v += w;
    v -= w;
    v *= 3.0;
    v.axpy(0.5, w);
    BOOST_CHECK_CLOSE(v.dot(w), 2 * numResDofs * 8.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SystemMatrixWithoutWells)
{
    // The merged blocks of a run without open wells, as uploaded by
    // ISTLSolverGPUSystem: B is 0 x N, C is N x 0 and D is 0 x 0.
    const std::size_t numCells = 2;
    // The merger keeps references to its inputs.
    const std::vector<WRMatrix<double>> wellB;
    const std::vector<RWMatrix<double>> wellC;
    const std::vector<WWMatrix<double>> wellD;
    const SparseTable<int> wellCells;
    const WellMatrixMerger<double> merger(numCells, wellB, wellC, wellD, wellCells);
    WRMatrix<double> B;
    RWMatrix<double> C;
    WWMatrix<double> D;
    merger.buildMatrices(B, C, D);

    const auto A = makeMatrix<RRMatrix<double>>(numCells, numCells, {{0, 1}, {0, 1}}, 10.0);
    const auto gpuA = GpuRRMatrix<double>::fromMatrix(A);
    auto gpuB = makeScalarGpuMatrix(B);
    auto gpuC = makeScalarGpuMatrix(C);
    auto gpuD = makeScalarGpuMatrix(D);
    std::vector<double> scratch;
    updateScalarGpuMatrix(B, gpuB, scratch);
    updateScalarGpuMatrix(C, gpuC, scratch);
    updateScalarGpuMatrix(D, gpuD, scratch);

    GpuSystemMatrix<double> gpuS;
    gpuS.A = &gpuA;
    gpuS.B = &gpuB;
    gpuS.C = &gpuC;
    gpuS.D = &gpuD;
    *gpuS.hostD = &D;

    ResVector<double> x(numCells);
    fill(x, 1.0);
    ResVector<double> y(numCells);
    fill(y, 0.3);
    GpuSystemVector<double> gpuX(numCells * numResDofs, 0);
    gpuX.res.copyFromHost(x);
    GpuSystemVector<double> gpuY(numCells * numResDofs, 0);
    gpuY.res.copyFromHost(y);

    gpuS.usmv(0.5, gpuX, gpuY);
    A.usmv(0.5, x, y);
    checkClose(y, gpuY.res);

    gpuS.mv(gpuX, gpuY);
    A.mv(x, y);
    checkClose(y, gpuY.res);
}
