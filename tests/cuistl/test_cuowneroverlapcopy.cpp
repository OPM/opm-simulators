/*
  Copyright 2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestCuOwnerOverlapCopy
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <memory>
#include <opm/simulators/linalg/cuistl/CuOwnerOverlapCopy.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <random>

bool
init_unit_test_func()
{
    return true;
}

int
main(int argc, char** argv)
{
    [[maybe_unused]] const auto& helper = Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}

BOOST_AUTO_TEST_CASE(TestProject)
{

    // We're going to have three points: Centered is owned, left and right is copied. We assume a periodic domain
    // ([0,1]/0~1)
    auto indexInfo = Dune::IndexInfoFromGrid<int, int>();
    indexInfo.addLocalIndex(std::make_tuple(0, 0, Dune::OwnerOverlapCopyAttributeSet::copy));
    indexInfo.addLocalIndex(std::make_tuple(1, 1, Dune::OwnerOverlapCopyAttributeSet::owner));
    indexInfo.addLocalIndex(std::make_tuple(2, 2, Dune::OwnerOverlapCopyAttributeSet::copy));

    auto ownerOverlapCopy = Dune::OwnerOverlapCopyCommunication<int>(indexInfo, MPI_COMM_WORLD);
    auto xCPU = std::vector<double> {{1.0, 2.0, 3.0}};
    auto xGPU = Opm::cuistl::CuVector<double>(xCPU);

    auto cuOwnerOverlapCopy
        = Opm::cuistl::CuOwnerOverlapCopy<double, 1, Dune::OwnerOverlapCopyCommunication<int>>(ownerOverlapCopy);

    cuOwnerOverlapCopy.project(xGPU);

    auto resultOfProject = xGPU.asStdVector();

    BOOST_CHECK_EQUAL(0.0, resultOfProject[0]);
    BOOST_CHECK_EQUAL(2.0, resultOfProject[1]);
    BOOST_CHECK_EQUAL(0.0, resultOfProject[2]);
}


BOOST_AUTO_TEST_CASE(TestDot)
{

    // We're going to have three points: Centered is owned, left and right is copied. We assume a periodic domain
    // ([0,1]/0~1)
    auto indexInfo = Dune::IndexInfoFromGrid<int, int>();
    indexInfo.addLocalIndex(std::make_tuple(0, 0, Dune::OwnerOverlapCopyAttributeSet::copy));
    indexInfo.addLocalIndex(std::make_tuple(1, 1, Dune::OwnerOverlapCopyAttributeSet::owner));
    indexInfo.addLocalIndex(std::make_tuple(2, 2, Dune::OwnerOverlapCopyAttributeSet::copy));

    indexInfo.addRemoteIndex(std::make_tuple(0, 0, Dune::OwnerOverlapCopyAttributeSet::copy));
    indexInfo.addRemoteIndex(std::make_tuple(0, 1, Dune::OwnerOverlapCopyAttributeSet::owner));
    indexInfo.addRemoteIndex(std::make_tuple(0, 2, Dune::OwnerOverlapCopyAttributeSet::copy));
    auto ownerOverlapCopy = Dune::OwnerOverlapCopyCommunication<int>(indexInfo, MPI_COMM_WORLD);
    auto xCPU = std::vector<double> {{1.0, 2.0, 3.0}};
    auto xGPU = Opm::cuistl::CuVector<double>(xCPU);

    auto cuOwnerOverlapCopy
        = Opm::cuistl::CuOwnerOverlapCopy<double, 1, Dune::OwnerOverlapCopyCommunication<int>>(ownerOverlapCopy);

    double outputDune = -1.0;
    auto xDune = xGPU.asDuneBlockVector<1>();
    ownerOverlapCopy.dot(xDune, xDune, outputDune);

    double output = -1.0;
    cuOwnerOverlapCopy.dot(xGPU, xGPU, output);


    BOOST_CHECK_EQUAL(outputDune, output);
    BOOST_CHECK_EQUAL(4.0, output);
}
