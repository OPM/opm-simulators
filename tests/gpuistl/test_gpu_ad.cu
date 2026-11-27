/*
  Copyright 2024 SINTEF AS
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

#define BOOST_TEST_MODULE TestGpuAD

#include <boost/test/unit_test.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <cuda_runtime.h>

#include "GpuTestHelpers.hpp"

namespace{
__global__ void instantiate_ad_object(double value, Opm::DenseAd::Evaluation<float, 3>* adObj){
    *adObj = Opm::DenseAd::Evaluation<float, 3>(value, 0);
}

//TODO add more comprenehsive AD tests

} // END EMPTY NAMESPACE


BOOST_AUTO_TEST_CASE(TestInstantiateADObject)
{
    using Evaluation = Opm::DenseAd::Evaluation<float, 3>;
    double testValue = 123.456;
    Evaluation cpuMadeAd = Evaluation(testValue, 0);

    // allocate space on GPU, run kernel, and move results back to the CPU
    const auto gpuMadeAd = Opm::gpuistl::test::runSingleThreadKernel<Evaluation>(instantiate_ad_object, testValue);

    // Check that the object made in a GPU kernel is equivalent to that made on the CPU
    BOOST_CHECK(cpuMadeAd == gpuMadeAd);
}
