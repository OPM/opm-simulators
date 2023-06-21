/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE OPM_test_PreconditionerFactory
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/ilufirstelement.hh>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/solvers.hh>

#include <fstream>
#include <iostream>


template <class X>
class NothingPreconditioner : public Dune::Preconditioner<X, X>
{
public:
    virtual void pre(X&, X&) override
    {
    }

    virtual void apply(X& v, const X& d) override
    {
        v = d;
    }

    virtual void post(X&) override
    {
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }
};


template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
testPrec(const Opm::PropertyTree& prm, const std::string& matrix_filename, const std::string& rhs_filename)
{
    using Matrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;
    Matrix matrix;
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        using M = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
        readMatrixMarket(reinterpret_cast<M&>(matrix), mfile); // Hack to avoid hassle
    }
    Vector rhs;
    {
        std::ifstream rhsfile(rhs_filename);
        if (!rhsfile) {
            throw std::runtime_error("Could not read rhs file");
        }
        readMatrixMarket(rhs, rhsfile);
    }
    using Operator = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    Operator op(matrix);
    using PrecFactory = Opm::PreconditionerFactory<Operator,Dune::Amg::SequentialInformation>;
    bool transpose = false;

    if(prm.get<std::string>("preconditioner.type") == "cprt"){
        transpose = true;
    }
    auto wc = [&matrix, transpose]()
    {
        return Opm::Amg::getQuasiImpesWeights<Matrix, Vector>(matrix, 1, transpose);
    };

    auto prec = PrecFactory::create(op, prm.get_child("preconditioner"), wc, 1);
    Dune::BiCGSTABSolver<Vector> solver(op, *prec, prm.get<double>("tol"), prm.get<int>("maxiter"), prm.get<int>("verbosity"));
    Vector x(rhs.size());
    Dune::InverseOperatorResult res;
    solver.apply(x, rhs, res);
    return x;
}

void test1(const Opm::PropertyTree& prm)
{
    constexpr int bz = 1;
    auto sol = testPrec<bz>(prm, "matr33.txt", "rhs3.txt");
    Dune::BlockVector<Dune::FieldVector<double, bz>> expected {-1.62493,
            -1.76435e-06,
            1.86991e-10,
            -458.542,
            2.28308e-06,
            -2.45341e-07,
            -1.48005,
            -5.02264e-07,
            -1.049e-05};
    BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
        }
    }
}

void test3(const Opm::PropertyTree& prm)
{
    constexpr int bz = 3;
    auto sol = testPrec<bz>(prm, "matr33.txt", "rhs3.txt");
    Dune::BlockVector<Dune::FieldVector<double, bz>> expected {{-1.62493, -1.76435e-06, 1.86991e-10},
            {-458.542, 2.28308e-06, -2.45341e-07},
                {-1.48005, -5.02264e-07, -1.049e-05}};
    BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
        }
    }
}



BOOST_AUTO_TEST_CASE(TestDefaultPreconditionerFactory)
{
    // Read parameters.
    Opm::PropertyTree prm("options_flexiblesolver.json");

    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    test3(prm);
}


template <int bz>
using M = Dune::BCRSMatrix<Opm::MatrixBlock<double, bz, bz>>;
template <int bz>
using V = Dune::BlockVector<Dune::FieldVector<double, bz>>;
template <int bz>
using O = Dune::MatrixAdapter<M<bz>, V<bz>, V<bz>>;
template <int bz>
using PF = Opm::PreconditionerFactory<O<bz>,Dune::Amg::SequentialInformation>;


BOOST_AUTO_TEST_CASE(TestAddingPreconditioner)
{
    // Read parameters.
    Opm::PropertyTree prm("options_flexiblesolver_simple.json");

    // Test with 1x1 block solvers.
    {
        constexpr int bz = 1;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::invalid_argument);
    }

    // Test with 3x3 block solvers.
    {
        constexpr int bz = 3;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::invalid_argument);
    }


    // Add preconditioner to factory for block size 1.
    PF<1>::addCreator("nothing", [](const O<1>&, const Opm::PropertyTree&, const std::function<V<1>()>&,
                                    std::size_t) {
            return Dune::wrapPreconditioner<NothingPreconditioner<V<1>>>();
        });


    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    {
        constexpr int bz = 3;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::invalid_argument);
    }

    // Add preconditioner to factory for block size 3.
    PF<3>::addCreator("nothing", [](const O<3>&, const Opm::PropertyTree&, const std::function<V<3>()>&,
                                    std::size_t) {
            return Dune::wrapPreconditioner<NothingPreconditioner<V<3>>>();
        });

    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    test3(prm);
}


template<class Mat, class Vec>
class RepeatingOperator : public Dune::AssembledLinearOperator<Mat, Vec, Vec>
{
public:
    using matrix_type = Mat;
    using domain_type = Vec;
    using range_type = Vec;
    using field_type = typename Vec::field_type;

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    RepeatingOperator(const Mat& matrix, const int repeats)
        : matrix_(matrix)
        , repeats_(repeats)
    {
    }

    // y = A*x;
    virtual void apply(const Vec& x, Vec& y) const override
    {
        y = 0;
        applyscaleadd(1.0, x, y);
    }

    // y += \alpha * A * x
    virtual void applyscaleadd(field_type alpha, const Vec& x, Vec& y) const override
    {
        Vec temp1 = x;
        Vec temp2 = x; // For size.
        temp2 = 0.0;
        for (int rr = 0; rr < repeats_; ++rr) {
            // mv below means: temp2 = matrix_ * temp1;
            matrix_.mv(temp1, temp2);
            temp1 = temp2;
        }
        temp2 *= alpha;
        y += temp2;
    }

    virtual const matrix_type& getmat() const override
    {
        return matrix_;
    }

protected:
    const Mat& matrix_;
    const int repeats_;
};


template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
testPrecRepeating(const Opm::PropertyTree& prm, const std::string& matrix_filename, const std::string& rhs_filename)
{
    using Matrix = M<bz>;
    using Vector = V<bz>;
    Matrix matrix;
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        using M = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        readMatrixMarket(reinterpret_cast<M&>(matrix), mfile); // Hack to avoid hassle
#pragma GCC diagnostic pop
    }
    Vector rhs;
    {
        std::ifstream rhsfile(rhs_filename);
        if (!rhsfile) {
            throw std::runtime_error("Could not read rhs file");
        }
        readMatrixMarket(rhs, rhsfile);
    }
    using Operator = RepeatingOperator<Matrix, Vector>;
    Operator op(matrix, 2);
    using PrecFactory = Opm::PreconditionerFactory<Operator,Dune::Amg::SequentialInformation>;

    // Add no-oppreconditioner to factory for block size 1.
    PrecFactory::addCreator("nothing", [](const Operator&, const Opm::PropertyTree&, const std::function<Vector()>&,
                                          std::size_t) {
        return Dune::wrapPreconditioner<NothingPreconditioner<Vector>>();
    });

    auto prec = PrecFactory::create(op, prm.get_child("preconditioner"));
    Dune::BiCGSTABSolver<Vector> solver(op, *prec, prm.get<double>("tol"), prm.get<int>("maxiter"), prm.get<int>("verbosity"));
    Vector x(rhs.size());
    Dune::InverseOperatorResult res;
    solver.apply(x, rhs, res);
    return x;
}

void test1rep(const Opm::PropertyTree& prm)
{
    constexpr int bz = 1;
    auto sol = testPrecRepeating<bz>(prm, "matr33rep.txt", "rhs3rep.txt");
    Dune::BlockVector<Dune::FieldVector<double, bz>> expected {0.285714285714286,
                                                               0.285714285714286,
                                                               0.285714285714286,
                                                               -0.214285714285714,
                                                               -0.214285714285714,
                                                               -0.214285714285714,
                                                               -0.214285714285714,
                                                               -0.214285714285714,
                                                               -0.214285714285714};
    BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
        }
    }
}

void test3rep(const Opm::PropertyTree& prm)
{
    constexpr int bz = 3;
    auto sol = testPrecRepeating<bz>(prm, "matr33rep.txt", "rhs3rep.txt");
    Dune::BlockVector<Dune::FieldVector<double, bz>> expected {
        {0.285714285714286, 0.285714285714286, 0.285714285714286},
        {-0.214285714285714, -0.214285714285714, -0.214285714285714},
        {-0.214285714285714, -0.214285714285714, -0.214285714285714}
    };
    BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
        }
    }
}


BOOST_AUTO_TEST_CASE(TestWithRepeatingOperator)
{
    // Read parameters.
    Opm::PropertyTree prm("options_flexiblesolver_simple.json");

    // Test with 1x1 block solvers.
    test1rep(prm);

    // Test with 3x3 block solvers.
    test3rep(prm);
}
