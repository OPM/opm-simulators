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

#define BOOST_TEST_MODULE TestHyprePreconditionerGPU
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/HyprePreconditioner.hpp>
#include <opm/simulators/linalg/gpuistl/set_device.hpp>
#include <dune/common/parallel/mpihelper.hh>

#include "MpiFixture.hpp"
#include "HyprePreconditionerTestHelper.hpp"

#if HAVE_CUDA || HAVE_HIP
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#endif

BOOST_GLOBAL_FIXTURE(MPIFixture);

// CPU input test function (specific to this test file)
inline void testHyprePreconditionerCpuInput(bool use_gpu)
{
    constexpr int N = 100; // 100x100 grid
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

    // Create matrix
    Matrix matrix;
    setupLaplace2d(N, matrix);

    testHyprePreconditionerImpl<Matrix, Vector>(matrix, use_gpu);
}

// Per-test fixture for HYPRE state isolation
struct HypreTestFixture {
    HypreTestFixture() {
        // Reset HYPRE state for each test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }

        // Re-initialize HYPRE for this test
#if HYPRE_RELEASE_NUMBER >= 22900
        HYPRE_Initialize();
#else
        HYPRE_Init();
#endif
    }

    ~HypreTestFixture() {
        // Clean state after test
        if (HYPRE_Initialized()) {
            HYPRE_Finalize();
        }
    }
};

BOOST_FIXTURE_TEST_CASE(TestHyprePreconditionerGPU_CPUInput, HypreTestFixture)
{
    // Test GPU backend with CPU input (CPU data transferred to GPU)
    testHyprePreconditionerCpuInput(true);
}

#if HAVE_CUDA || HAVE_HIP
inline void testHyprePreconditionerGpuInput(bool use_gpu)
{
    using namespace Opm::gpuistl;

    constexpr int N = 100; // 100x100 grid
    using CpuMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using GpuMatrixType = GpuSparseMatrix<double>;
    using GpuVectorType = GpuVector<double>;

    // Create matrix on CPU first
    CpuMatrix cpu_matrix;
    setupLaplace2d(N, cpu_matrix);

    // Convert to GPU matrix
    GpuMatrixType gpu_matrix = GpuMatrixType::fromMatrix(cpu_matrix);

    testHyprePreconditionerImpl<GpuMatrixType, GpuVectorType>(gpu_matrix, use_gpu);
}

BOOST_FIXTURE_TEST_CASE(TestHyprePreconditionerGPU_GPUInput, HypreTestFixture)
{
    // Test GPU backend with GPU input (GPU data used directly)
    testHyprePreconditionerGpuInput(true);
}

BOOST_FIXTURE_TEST_CASE(TestHyprePreconditionerCPU_GPUInput, HypreTestFixture)
{
    // Test CPU backend with GPU input (data transfered to CPU for use with Hypre)
    testHyprePreconditionerGpuInput(false);
}
#endif

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
