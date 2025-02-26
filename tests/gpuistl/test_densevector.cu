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
#include <boost/test/tools/old/interface.hpp>
#include <config.h>
#include <stdexcept>

#define BOOST_TEST_MODULE TestDenseVector

#include <cuda.h>
#include <cuda_runtime.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/dense/DenseVector.hpp>
#include <opm/simulators/linalg/gpuistl/dense/FieldVector.hpp>

BOOST_AUTO_TEST_CASE(TestDenseCreation) 
{
    //auto someVec = Opm::gpuistl::dense::DenseVector<float>(10);
}
