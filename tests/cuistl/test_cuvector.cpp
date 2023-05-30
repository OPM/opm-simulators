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

#define BOOST_TEST_MODULE TestCuVector

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(TestDocumentedUsage)
{
    auto someDataOnCPU = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});

    auto dataOnGPU = ::Opm::cuistl::CuVector<double>(someDataOnCPU);

    // Multiply by 4.0:
    dataOnGPU *= 4.0;

    // Get data back on CPU in another vector:
    auto stdVectorOnCPU = dataOnGPU.asStdVector();

    std::vector<double> correctVector(someDataOnCPU.size());

    std::transform(someDataOnCPU.begin(), someDataOnCPU.end(), correctVector.begin(), [](double x) { return 4 * x; });
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVectorOnCPU.begin(), stdVectorOnCPU.end(), correctVector.begin(), correctVector.end());
}

BOOST_AUTO_TEST_CASE(TestConstructionSize)
{
    const int numberOfElements = 1234;
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(numberOfElements);
    BOOST_CHECK_EQUAL(numberOfElements, vectorOnGPU.dim());
}

BOOST_AUTO_TEST_CASE(TestCopyFromHostConstructor)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    BOOST_CHECK_EQUAL(data.size(), vectorOnGPU.dim());
    std::vector<double> buffer(data.size(), 0.0);
    vectorOnGPU.copyToHost(buffer.data(), buffer.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(buffer.begin(), buffer.end(), data.begin(), data.end());
}


BOOST_AUTO_TEST_CASE(TestCopyFromHostFunction)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.size());
    BOOST_CHECK_EQUAL(data.size(), vectorOnGPU.dim());
    vectorOnGPU.copyFromHost(data.data(), data.size());
    std::vector<double> buffer(data.size(), 0.0);
    vectorOnGPU.copyToHost(buffer.data(), buffer.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(buffer.begin(), buffer.end(), data.begin(), data.end());
}


BOOST_AUTO_TEST_CASE(TestCopyFromBvector)
{
    auto blockVector = Dune::BlockVector<Dune::FieldVector<double, 2>> {{{42, 43}, {44, 45}, {46, 47}}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(blockVector.dim());
    vectorOnGPU.copyFromHost(blockVector);
    std::vector<double> buffer(vectorOnGPU.dim());
    vectorOnGPU.copyToHost(buffer.data(), buffer.size());

    BOOST_CHECK_EQUAL_COLLECTIONS(
        buffer.begin(), buffer.end(), &blockVector[0][0], &blockVector[0][0] + blockVector.dim());
}

BOOST_AUTO_TEST_CASE(TestCopyToBvector)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7, 8, 9}};
    auto blockVector = Dune::BlockVector<Dune::FieldVector<double, 3>>(3);
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    vectorOnGPU.copyToHost(blockVector);


    BOOST_CHECK_EQUAL_COLLECTIONS(data.begin(), data.end(), &blockVector[0][0], &blockVector[0][0] + blockVector.dim());
}

BOOST_AUTO_TEST_CASE(TestDataPointer)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7, 8, 9}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());

    std::vector<double> buffer(data.size(), 0.0);
    cudaMemcpy(buffer.data(), vectorOnGPU.data(), sizeof(double) * data.size(), cudaMemcpyDeviceToHost);
    BOOST_CHECK_EQUAL_COLLECTIONS(data.begin(), data.end(), buffer.begin(), buffer.end());
}

BOOST_AUTO_TEST_CASE(TestCopyScalarMultiply)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    BOOST_CHECK_EQUAL(data.size(), vectorOnGPU.dim());
    const double scalar = 42.25;
    vectorOnGPU *= scalar;
    std::vector<double> buffer(data.size(), 0.0);
    vectorOnGPU.copyToHost(buffer.data(), buffer.size());

    for (size_t i = 0; i < buffer.size(); ++i) {
        BOOST_CHECK_EQUAL(buffer[i], scalar * data[i]);
    }
}

BOOST_AUTO_TEST_CASE(TestTwoNorm)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    auto twoNorm = vectorOnGPU.two_norm();

    double correctAnswer = 0.0;
    for (double d : data) {
        correctAnswer += d * d;
    }
    correctAnswer = std::sqrt(correctAnswer);
    BOOST_CHECK_EQUAL(correctAnswer, twoNorm);
}

BOOST_AUTO_TEST_CASE(TestDot)
{
    std::vector<double> dataA {{1, 2, 3, 4, 5, 6, 7}};
    std::vector<double> dataB {{8, 9, 10, 11, 12, 13, 14}};
    auto vectorOnGPUA = Opm::cuistl::CuVector<double>(dataA.data(), dataA.size());
    auto vectorOnGPUB = Opm::cuistl::CuVector<double>(dataB.data(), dataB.size());
    auto dot = vectorOnGPUA.dot(vectorOnGPUB);

    double correctAnswer = 0.0;
    for (size_t i = 0; i < dataA.size(); ++i) {
        correctAnswer += dataA[i] * dataB[i];
    }
    correctAnswer = correctAnswer;
    BOOST_CHECK_EQUAL(correctAnswer, dot);
}

BOOST_AUTO_TEST_CASE(Assigment)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    vectorOnGPU = 10.0;
    vectorOnGPU.copyToHost(data.data(), data.size());

    for (double x : data) {
        BOOST_CHECK_EQUAL(10.0, x);
    }
}


BOOST_AUTO_TEST_CASE(CopyAssignment)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto vectorOnGPU = Opm::cuistl::CuVector<double>(data.data(), data.size());
    vectorOnGPU.copyToHost(data.data(), data.size());
    auto vectorOnGPUB = Opm::cuistl::CuVector<double>(data.size());
    vectorOnGPUB = 4.0;
    vectorOnGPUB = vectorOnGPU;

    std::vector<double> output(data.size());
    vectorOnGPUB.copyToHost(output.data(), output.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(output.begin(), output.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(RandomVectors)
{

    using GVector = Opm::cuistl::CuVector<double>;
    std::srand(0);
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution(-100.0, 100.0);
    std::uniform_real_distribution<double> distribution01(.0, 1.0);

    const size_t N = 1000;

    const size_t retries = 100;

    for (size_t retry = 0; retry < retries; ++retry) {
        std::vector<double> a(N);
        std::vector<double> b(N);

        for (size_t i = 0; i < N; ++i) {
            a[i] = distribution(generator);
            b[i] = distribution(generator);
        }

        auto aGPU = GVector(a);
        auto bGPU = GVector(b);

        aGPU += bGPU;

        auto aOutputPlus = aGPU.asStdVector();

        for (size_t i = 0; i < N; ++i) {
            BOOST_CHECK_EQUAL(aOutputPlus[i], a[i] + b[i]);
        }

        aGPU = GVector(a);
        aGPU -= bGPU;

        auto aOutputMinus = aGPU.asStdVector();

        for (size_t i = 0; i < N; ++i) {
            BOOST_CHECK_EQUAL(aOutputMinus[i], a[i] - b[i]);
        }


        aGPU = GVector(a);
        auto scalar = distribution(generator);
        aGPU *= scalar;
        auto aOutputScalar = aGPU.asStdVector();
        for (size_t i = 0; i < N; ++i) {
            BOOST_CHECK_EQUAL(aOutputScalar[i], scalar * a[i]);
        }

        aGPU = GVector(a);
        aGPU.axpy(scalar, bGPU);
        auto aOutputSaxypy = aGPU.asStdVector();
        for (size_t i = 0; i < N; ++i) {
            BOOST_CHECK_CLOSE(aOutputSaxypy[i], a[i] + scalar * b[i], 1e-10);
        }

        aGPU = GVector(a);
        auto dotted = aGPU.dot(bGPU);
        double correct = 0.0;
        for (size_t i = 0; i < N; ++i) {
            correct += a[i] * b[i];
        }

        BOOST_CHECK_CLOSE(dotted, correct, 1e-10);

        aGPU = GVector(a);
        auto twoNorm = aGPU.two_norm();
        double correctTwoNorm = 0.0;
        for (size_t i = 0; i < N; ++i) {
            correctTwoNorm += a[i] * a[i];
        }
        correctTwoNorm = std::sqrt(correctTwoNorm);

        BOOST_CHECK_CLOSE(twoNorm, correctTwoNorm, 1e-12);

        aGPU = GVector(a);
        std::vector<int> indexSet;
        const double rejectCriteria = 0.2;
        for (size_t i = 0; i < N; ++i) {
            const auto reject = distribution01(generator);
            if (reject < rejectCriteria) {
                indexSet.push_back(i);
            }
        }
        auto indexSetGPU = Opm::cuistl::CuVector<int>(indexSet);

        aGPU.setZeroAtIndexSet(indexSetGPU);
        auto projectedA = aGPU.asStdVector();
        for (size_t i = 0; i < N; ++i) {
            // Yeah, O(N^2) so sue me
            bool found = std::find(indexSet.begin(), indexSet.end(), i) != indexSet.end();
            if (found) {
                BOOST_CHECK_EQUAL(projectedA[i], 0);
            } else {
                BOOST_CHECK_EQUAL(projectedA[i], a[i]);
            }
        }

        aGPU = GVector(a);
        auto twoNormAtIndices = aGPU.two_norm(indexSetGPU);

        double correctTwoNormAtIndices = 0.0;
        for (size_t i = 0; i < indexSet.size(); ++i) {
            correctTwoNormAtIndices += a[indexSet[i]] * a[indexSet[i]];
        }
        correctTwoNormAtIndices = std::sqrt(correctTwoNormAtIndices);

        BOOST_CHECK_CLOSE(correctTwoNormAtIndices, twoNormAtIndices, 1e-13);
    }
}
