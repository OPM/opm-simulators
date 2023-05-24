/*
  Copyright 2022-2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestCudaCheckLastError

#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_check_last_error.hpp>


BOOST_AUTO_TEST_CASE(TestNoThrowLastError)
{
    BOOST_CHECK_NO_THROW(OPM_CUDA_CHECK_LAST_ERROR;);
    BOOST_CHECK_NO_THROW(OPM_CUDA_CHECK_LAST_ERROR_IF_DEBUG;);
}


BOOST_AUTO_TEST_CASE(TestNoThrowDeviceSynchronize)
{
    BOOST_CHECK_NO_THROW(OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE;);
    BOOST_CHECK_NO_THROW(OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG;);
}
