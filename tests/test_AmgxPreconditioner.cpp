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

#define BOOST_TEST_MODULE TestAmgxPreconditioner
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/AmgxPreconditioner.hpp>

#include "AmgxPreconditionerTestHelper.hpp"


//==============================================================================
// Test preconditioner solve CPU input
//==============================================================================
BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerCpuInputMatFloatVecFloat)
{
    testAmgxPreconditionerCpuInput<float, float>();
}

BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerCpuInputMatDoubleVecDouble)
{
    testAmgxPreconditionerCpuInput<double, double>();
}

// Note: Mixed precision tests disabled due to CUDA limitation
// BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerMatFloatVecDouble)
// {
//     testAmgxPreconditionerCpuInput<float, double>();
// }


//==============================================================================
// Test preconditioner update CPU input
//==============================================================================
BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerUpdateCpuInputFloat)
{
    testAmgxPreconditionerUpdateCpuInput<float, float>();
}

BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerUpdateCpuInputDouble)
{
    testAmgxPreconditionerUpdateCpuInput<double, double>();
}


//==============================================================================
// Test preconditioner solve GPU input
//==============================================================================

BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerGpuInputMatFloatVecFloat)
{
    testAmgxPreconditionerGpuInput<float, float>();
}

BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerGpuInputMatDoubleVecDouble)
{
    testAmgxPreconditionerGpuInput<double, double>();
}


//==============================================================================
// Test preconditioner update GPU input
//==============================================================================
BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerUpdateGpuInputMatFloatVecFloat)
{
    testAmgxPreconditionerUpdateGpuInput<float, float>();
}

BOOST_AUTO_TEST_CASE(TestAmgxPreconditionerUpdateGpuInputMatDoubleVecDouble)
{
    testAmgxPreconditionerUpdateGpuInput<double, double>();
}


bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    AMGX_SAFE_CALL(AMGX_initialize());

    int result = boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);

    AMGX_SAFE_CALL(AMGX_finalize());

    return result;
}
