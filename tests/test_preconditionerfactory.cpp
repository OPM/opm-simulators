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

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/solvers.hh>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

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
testPrec(const boost::property_tree::ptree& prm, const std::string& matrix_filename, const std::string& rhs_filename)
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;
    Matrix matrix;
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        readMatrixMarket(matrix, mfile);
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
    using PrecFactory = Dune::PreconditionerFactory<Operator>;
    auto prec = PrecFactory::create(op, prm.get_child("preconditioner"));
    Dune::BiCGSTABSolver<Vector> solver(op, *prec, prm.get<double>("tol"), prm.get<int>("maxiter"), prm.get<int>("verbosity"));
    Vector x(rhs.size());
    Dune::InverseOperatorResult res;
    solver.apply(x, rhs, res);
    return x;
}

namespace pt = boost::property_tree;

void test1(const pt::ptree& prm)
{
    const int bz = 1;
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

void test3(const pt::ptree& prm)
{
    const int bz = 3;
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
    pt::ptree prm;

    // Read parameters.
    {
        std::ifstream file("options_flexiblesolver.json");
        pt::read_json(file, prm);
    }

    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    test3(prm);
}


template <int bz>
using M = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
template <int bz>
using V = Dune::BlockVector<Dune::FieldVector<double, bz>>;
template <int bz>
using O = Dune::MatrixAdapter<M<bz>, V<bz>, V<bz>>;
template <int bz>
using PF = Dune::PreconditionerFactory<O<bz>>;


BOOST_AUTO_TEST_CASE(TestAddingPreconditioner)
{
    namespace pt = boost::property_tree;
    pt::ptree prm;

    // Read parameters.
    {
        std::ifstream file("options_flexiblesolver_simple.json"); // Requests "nothing" for preconditioner type.
        pt::read_json(file, prm);
    }

    // Test with 1x1 block solvers.
    {
        const int bz = 1;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::runtime_error);
    }

    // Test with 3x3 block solvers.
    {
        const int bz = 3;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::runtime_error);
    }


    // Add preconditioner to factory for block size 1.
    PF<1>::addCreator("nothing", [](const O<1>&, const pt::ptree&) {
            return Dune::wrapPreconditioner<NothingPreconditioner<V<1>>>();
        });


    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    {
        const int bz = 3;
        BOOST_CHECK_THROW(testPrec<bz>(prm, "matr33.txt", "rhs3.txt"), std::runtime_error);
    }

    // Add preconditioner to factory for block size 3.
    PF<3>::addCreator("nothing", [](const O<3>&, const pt::ptree&) {
            return Dune::wrapPreconditioner<NothingPreconditioner<V<3>>>();
        });

    // Test with 1x1 block solvers.
    test1(prm);

    // Test with 3x3 block solvers.
    test3(prm);
}

#else

// Do nothing if we do not have at least Dune 2.6.
BOOST_AUTO_TEST_CASE(DummyTest)
{
    BOOST_REQUIRE(true);
}

#endif
