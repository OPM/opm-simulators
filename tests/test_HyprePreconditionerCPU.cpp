/*
  Copyright 2024 SINTEF AS
  Copyright 2024 Equinor ASA

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

#define BOOST_TEST_MODULE TestHyprePreconditionerCPU
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "HypreTestHelper.hpp"
#include "MpiFixture.hpp"


#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>


using namespace HypreTestHelpers;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_FIXTURE_TEST_CASE(TestHyprePreconditioner_CpuInputCpuBackend, HypreTestFixture)
{
    constexpr int N = 100; // 100x100 grid
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);

    // Test CPU input with CPU backend
    testHyprePreconditioner<Matrix, Vector>(matrix, false);
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    int result = boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);

    return result;
}
