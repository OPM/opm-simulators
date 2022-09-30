/*
  Copyright SINTEF AS 2022

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

#define BOOST_TEST_MODULE TestCuSeqILU0
// TODO: Check which ones are needed:


#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/amgcpr.hh>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/matrixblock.hh>

#include <boost/test/unit_test.hpp>
#include <chrono>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <memory>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/impl/cuda_safe_call.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(TestFiniteDifference1D)
{
    // Here we will test a simple 1D finite difference scheme for
    // the Laplace equation:
    //
    //    -\Delta u = f on [0,1]
    //
    // Using a central difference approximation of \Delta u, this can
    // be approximated by
    //
    //    -(u_{i+1}-2u_i+u_{i-1})/Dx^2 = f(x_i)
    //
    // giving rise to the matrix
    //
    //     -2  1  0  0 ... 0  0
    //      1 -2  1  0  0 ... 0
    //      ....
    //      0  0  0  ...1 -2  1
    //      0  0  0  ...   1 -2
    for (int k = 5; k < 25; ++k) {
        const int N = 2 << k;
        const int nonZeroes = N * 3 - 2;
        using M = Dune::FieldMatrix<double, 1, 1>;
        using SpMatrix = Dune::BCRSMatrix<M>;
        using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
        using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<double>, Opm::cuistl::CuVector<double>>;
        using C = Dune::Amg::SequentialInformation;

        SpMatrix B(N, N, nonZeroes, SpMatrix::row_wise);
        for (auto row = B.createbegin(); row != B.createend(); ++row) {
            // Add nonzeros for left neighbour, diagonal and right neighbour
            if (row.index() > 0) {
                row.insert(row.index() - 1);
            }
            row.insert(row.index());
            if (row.index() < B.N() - 1) {
                row.insert(row.index() + 1);
            }
        }
        // This might not be the most elegant way of filling in a Dune sparse matrix, but it works.
        for (int i = 0; i < N; ++i) {
            B[i][i] = -2;
            if (i < N - 1) {
                B[i][i + 1] = 1;
            }

            if (i > 0) {
                B[i][i - 1] = 1;
            }
        }


        Vector inputVector(N);
        Vector outputVector(N);
        Opm::cuistl::CuVector<double> cuInputVector(N);
        Opm::cuistl::CuVector<double> cuOutputVector(N);
        double totalDune = 0.0;
        double totalcuILUAdapt = 0.0;
        double totalcuILU = 0.0;
        const int repeats = 10;
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(1.0, 10.0);
        for (int repeat = 0; repeat < repeats; ++repeat) {
            for (int i = 0; i < N; ++i) {
                inputVector[i][0] = dist(mt);
            }
            cuInputVector.copyFromHost(inputVector);

            auto duneILU = Opm::ParallelOverlappingILU0<SpMatrix, Vector, Vector, C>(B, 0, 1.0, Opm::MILU_VARIANT::ILU);
            // Dune::SeqILU<SpMatrix, Vector, Vector>(B, 1.0);
            auto startDune = std::chrono::high_resolution_clock::now();
            for (int q = 0; q < repeats; ++q) {
                duneILU.apply(outputVector, inputVector);
            }

            auto stopDune = std::chrono::high_resolution_clock::now();
            auto durationDune = std::chrono::duration_cast<std::chrono::microseconds>(stopDune - startDune);


            auto cuILUAdapt
                = Opm::cuistl::PreconditionerAdapter<Vector, Vector, CuILU0>(std::make_shared<CuILU0>(B, 1.0));
            auto startcuILUAdapt = std::chrono::high_resolution_clock::now();
            for (int q = 0; q < repeats; ++q) {
                cuILUAdapt.apply(outputVector, inputVector);
            }
            auto stopcuILUAdapt = std::chrono::high_resolution_clock::now();
            auto durationcuILUAdapt
                = std::chrono::duration_cast<std::chrono::microseconds>(stopcuILUAdapt - startcuILUAdapt);

            auto cuILU = CuILU0(B, 1.0);
            auto startcuILU = std::chrono::high_resolution_clock::now();
            for (int q = 0; q < repeats; ++q) {
                cuILU.apply(cuOutputVector, cuInputVector);
            }
            auto stopcuILU = std::chrono::high_resolution_clock::now();
            auto durationcuILU = std::chrono::duration_cast<std::chrono::microseconds>(stopcuILU - startcuILU);
            totalDune += durationDune.count();
            totalcuILU += durationcuILU.count();
            totalcuILUAdapt += durationcuILUAdapt.count();
        }
        std::cout << N << " " << totalDune << " " << totalcuILUAdapt << " " << totalcuILU << std::endl;
    }
}
