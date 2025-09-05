/*
  Copyright 2025 Equinor ASA

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

#define BOOST_TEST_MODULE TestHypreInterfaceGPU
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "../MpiFixture.hpp"
#include "../HypreTestHelper.hpp"

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/HypreInterface.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>



using namespace Opm::gpuistl;
using namespace HypreTestHelpers;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_FIXTURE_TEST_CASE(TestResourceManagement, HypreTestFixture)
{
    testResourceManagement(true);
}

// CPU Input + GPU Backend tests
BOOST_FIXTURE_TEST_CASE(TestVectorTransfer_CpuInputGpuBackend, HypreTestFixture)
{
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    testVectorTransfer<Vector>(true);
}

BOOST_FIXTURE_TEST_CASE(TestMatrixTransfer_CpuInputGpuBackend, HypreTestFixture)
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    testMatrixTransfer<Matrix>(true);
}

// GPU Input + GPU Backend tests
BOOST_FIXTURE_TEST_CASE(TestVectorTransfer_GpuInputGpuBackend, HypreTestFixture)
{
    testVectorTransfer<GpuVector<double>>(true);
}

BOOST_FIXTURE_TEST_CASE(TestMatrixTransfer_GpuInputGpuBackend, HypreTestFixture)
{
    testMatrixTransfer<GpuSparseMatrixWrapper<double>>(true);
}

// GPU Input + CPU Backend tests
BOOST_FIXTURE_TEST_CASE(TestVectorTransfer_GpuInputCpuBackend, HypreTestFixture)
{
    testVectorTransfer<GpuVector<double>>(false);
}

BOOST_FIXTURE_TEST_CASE(TestMatrixTransfer_GpuInputCpuBackend, HypreTestFixture)
{
    testMatrixTransfer<GpuSparseMatrixWrapper<double>>(false);
}

// Error handling test run last to avoid interfering with other tests
BOOST_FIXTURE_TEST_CASE(TestErrorHandling, HypreTestFixture)
{
    testErrorHandling(true);
}

bool
init_unit_test_func()
{
    return true;
}

int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    int result = boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);

    return result;
}
