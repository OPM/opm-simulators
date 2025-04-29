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
#include <cstddef>
#include <stdexcept>

#define BOOST_TEST_MODULE TestPreconditionerFactoryGPU

#include <boost/test/unit_test.hpp>
#include <cuda.h>
#include <cuda_runtime.h>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

using Scalar = double;
constexpr static int Dim = 3;
using MatrixTypeCPU = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, Dim, Dim>>;
// using VectorType = Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>;
using CommSeq = Dune::Amg::SequentialInformation;
// using GpuOperator = Dune::MatrixAdapter<MatrixType,
//                                    VectorType,
//                                    VectorType>;

using GpuMatrixType = Opm::gpuistl::GpuSparseMatrix<Scalar>;
using GpuVectorType = Opm::gpuistl::GpuVector<Scalar>;
using GpuOperatorType = Dune::MatrixAdapter<GpuMatrixType, GpuVectorType, GpuVectorType>;
using FactoryTypeGpu = Opm::PreconditionerFactory<GpuOperatorType, CommSeq>;

BOOST_AUTO_TEST_CASE(TestMatrixAdapter)
{
    const std::size_t N = 10;
    const std::size_t nonZeroes = N * 3 - 2;
    MatrixTypeCPU matrix(N, N, nonZeroes, MatrixTypeCPU::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        // Add nonzeros for left neighbour, diagonal and right neighbour
        if (row.index() > 0) {
            row.insert(row.index() - 1);
        }
        row.insert(row.index());
        if (row.index() < matrix.N() - 1) {
            row.insert(row.index() + 1);
        }
    }

    for (std::size_t i = 0; i < N; ++i) {
        matrix[i][i] = -2;
        if (i < int(N) - 1) {
            matrix[i][i + 1] = 1;
        }

        if (i > 0) {
            matrix[i][i - 1] = 1;
        }
    }
    GpuMatrixType gpuMatrix = GpuMatrixType::fromMatrix(matrix);
    Dune::MatrixAdapter<GpuMatrixType, GpuVectorType, GpuVectorType> gpuOp(gpuMatrix);

    Opm::PropertyTree prm;
    prm.put(std::string("preconditioner"), std::string("DILU"));
    auto weightsCalculator = []() { return GpuVectorType(10); };

    // TODO: Once the GPU preconditioners are implemented, we can remove the check throw
    // and actually test the preconditioner creation.
    BOOST_CHECK_THROW(FactoryTypeGpu::create(gpuOp, prm, weightsCalculator), std::invalid_argument);
}
