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

#define BOOST_TEST_MODULE TestConvertToFloatAdapter
#define BOOST_TEST_NO_MAIN


#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <limits>
#include <memory>
#include <opm/simulators/linalg/cuistl/PreconditionerConvertFieldTypeAdapter.hpp>


using XDouble = Dune::BlockVector<Dune::FieldVector<double, 2>>;
using MDouble = Dune::FieldMatrix<double, 2, 2>;
using SpMatrixDouble = Dune::BCRSMatrix<MDouble>;
using XFloat = Dune::BlockVector<Dune::FieldVector<float, 2>>;
using MFloat = Dune::FieldMatrix<float, 2, 2>;
using SpMatrixFloat = Dune::BCRSMatrix<MFloat>;
namespace
{
class TestPreconditioner : Dune::PreconditionerWithUpdate<XFloat, XFloat>
{
public:
    using range_type = XFloat;
    using domain_type = XFloat;
    TestPreconditioner(const SpMatrixFloat& matrix,
                       const XDouble& expectedInput,
                       const SpMatrixDouble& expectedMatrix,
                       const XDouble& expectedOutputVector)
        : m_matrix(matrix)
        , m_expectedInput(expectedInput)
        , m_expectedMatrix(expectedMatrix)
        , m_expectedOutputVector(expectedOutputVector)
    {
    }

    virtual void pre([[maybe_unused]] XFloat& x, [[maybe_unused]] XFloat& b) override
    {
    }


    virtual void apply([[maybe_unused]] XFloat& v, const XFloat& d) override
    {
        // Make sure the correct input is copied
        for (size_t i = 0; i < d.N(); ++i) {
            for (size_t j = 0; j < d[i].N(); ++j) {
                BOOST_CHECK_EQUAL(d[i][j], float(m_expectedInput[i][j]));
                v[i][j] = float(m_expectedOutputVector[i][j]);
            }
        }

        // make sure we get the correct matrix
        BOOST_CHECK_EQUAL(m_expectedMatrix.N(), m_matrix.N());
        BOOST_CHECK_EQUAL(m_expectedMatrix.nonzeroes(), m_matrix.nonzeroes());

        for (auto row = m_matrix.begin(); row != m_matrix.end(); ++row) {
            for (auto column = row->begin(); column != row->end(); ++column) {
                for (int i = 0; i < MFloat::rows; ++i) {
                    for (int j = 0; j < MFloat::cols; ++j) {
                        BOOST_CHECK_EQUAL(float(m_expectedMatrix[i][j][i][j]), m_matrix[i][j][i][j]);
                    }
                }
            }
        }
    }

    virtual void post([[maybe_unused]] XFloat& x) override
    {
    }

    virtual void update() override
    {
    }


    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    static constexpr bool shouldCallPre()
    {
        return false;
    }
    static constexpr bool shouldCallPost()
    {
        return false;
    }


private:
    const SpMatrixFloat& m_matrix;
    const XDouble& m_expectedInput;
    const SpMatrixDouble& m_expectedMatrix;
    const XDouble& m_expectedOutputVector;
};
} // namespace

using NumericTypes = boost::mpl::list<double, float>;

BOOST_AUTO_TEST_CASE(TestFiniteDifference1D)
{


    const int N = 5;
    const int nonZeroes = N * 3 - 2;

    SpMatrixDouble B(N, N, nonZeroes, SpMatrixDouble::row_wise);
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



    // check for the standard basis {e_i}
    // (e_i=(0,...,0, 1 (i-th place), 0, ..., 0))
    for (int i = 0; i < N; ++i) {
        XDouble inputVector(N);
        XDouble outputVector(N);
        XDouble expectedOutputVector(N);
        expectedOutputVector[i][0] = 42.0;
        expectedOutputVector[i][1] = 43.0;
        inputVector[i][0] = 1.0;
        auto converter
            = Opm::cuistl::PreconditionerConvertFieldTypeAdapter<TestPreconditioner, SpMatrixDouble, XDouble, XDouble>(
                B);
        auto underlyingPreconditioner = std::make_shared<TestPreconditioner>(
            converter.getConvertedMatrix(), inputVector, B, expectedOutputVector);
        converter.setUnderlyingPreconditioner(underlyingPreconditioner);
        converter.apply(outputVector, inputVector);

        for (size_t j = 0; j < expectedOutputVector.N(); ++j) {
            for (size_t k = 0; k < expectedOutputVector[i].N(); ++k) {
                BOOST_CHECK_EQUAL(outputVector[j][k], float(expectedOutputVector[j][k]));
            }
        }
    }
}

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
