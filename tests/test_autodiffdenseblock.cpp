/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil AS.

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE AutoDiffDenseBlockTest

#include <opm/autodiff/AutoDiffDenseBlock.hpp>

#include <boost/test/unit_test.hpp>

using namespace Opm;


typedef AutoDiffDenseBlock<double, 3> ADD;
typedef ADD::Value V;
typedef ADD::Derivative D;



BOOST_AUTO_TEST_CASE(ConstantInitialisation)
{
    V v(3);
    v << 0.2, 1.2, 13.4;

    ADD a = ADD::constant(v);
    BOOST_CHECK(a.value().matrix() == v.matrix());

    const D& da = a.derivative();
    BOOST_CHECK((da == 0.0).all());
}



BOOST_AUTO_TEST_CASE(VariableInitialisation)
{
    V v(3);
    v << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADD x = ADD::variable(FirstVar, v);

    BOOST_CHECK(x.value().matrix() == v.matrix());

    const D& dx = x.derivative();
    BOOST_CHECK((dx.col(FirstVar) == 1.0).all());
    BOOST_CHECK((dx.col(SecondVar) == 0.0).all());
    BOOST_CHECK((dx.col(ThirdVar) == 0.0).all());
}



BOOST_AUTO_TEST_CASE(FunctionInitialisation)
{
    V v(3);
    v << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    D jac = D::Zero(3, 3);
    jac(0, 0) = -1.0;
    jac(1, 2) = -1.0;
    jac(0, 2) = -1.0;

    V v_copy(v);
    D jac_copy(jac);
    ADD f = ADD::function(std::move(v_copy), std::move(jac_copy));

    BOOST_CHECK(f.value().matrix() == v.matrix());
    BOOST_CHECK(f.derivative().matrix() == jac.matrix());

    jac(0, 2) = 23.0;

    BOOST_CHECK(f.derivative().matrix() != jac.matrix());
}



BOOST_AUTO_TEST_CASE(Addition)
{
    V va(3);
    va << 0.2, 1.2, 13.4;

    V vx(3);
    vx << 1.0, 2.2, 3.4;

    enum { FirstVar = 0, SecondVar = 1, ThirdVar = 2 };

    ADD a = ADD::constant(va);
    ADD x = ADD::variable(FirstVar, vx);

    ADD xpx = x + x;

    BOOST_CHECK((xpx.value() == 2*x.value()).all());
    BOOST_CHECK((xpx.derivative() == 2*x.derivative()).all());

    V  r = 2*x.value() + a.value();
    ADD xpxpa = x + x + a;
    BOOST_CHECK(xpxpa.value().matrix() == r.matrix());
    BOOST_CHECK((xpxpa.derivative() == 2*x.derivative()).all());
}



BOOST_AUTO_TEST_CASE(AssignAddSubtractOperators)
{
    // Basic testing of += and -=.
    V vx(3);
    vx << 0.2, 1.2, 13.4;

    V vy(3);
    vy << 1.0, 2.2, 3.4;

    std::vector<V> vals{ vx, vy };
    std::vector<ADD> vars = ADD::variables(vals);

    const ADD x = vars[0];
    const ADD y = vars[1];

    ADD z = x;
    z += y;
    ADD sum = x + y;
    const double tolerance = 1e-14;
    BOOST_CHECK(z.value().isApprox(sum.value(), tolerance));
    BOOST_CHECK(z.derivative().isApprox(sum.derivative(), tolerance));
    z -= y;
    BOOST_CHECK(z.value().isApprox(x.value(), tolerance));
    BOOST_CHECK(z.derivative().isApprox(x.derivative(), tolerance));

    // Testing the case when the left hand side is constant.
    ADD yconst = ADD::constant(vy);
    z = yconst;
    z -= x;
    ADD diff = yconst - x;
    BOOST_CHECK(z.value().isApprox(diff.value(), tolerance));
    BOOST_CHECK(z.derivative().isApprox(diff.derivative(), tolerance));
    z += x;
    BOOST_CHECK(z.value().isApprox(yconst.value(), tolerance));
    BOOST_CHECK(z.derivative().isApprox(D::Zero(3, 3)));
}



BOOST_AUTO_TEST_CASE(Multiplication)
{
    V vx(3);
    vx << 0.2, 1.2, 13.4;

    V vy(3);
    vy << 1.0, 2.2, 3.4;

    std::vector<V> vals{ vx, vy };
    std::vector<ADD> vars = ADD::variables(vals);

    const ADD x = vars[0];
    const ADD y = vars[1];

    const ADD xxy = x * x * y;

    const double tolerance = 1e-14;
    BOOST_CHECK(xxy.value().isApprox(vx * vx * vy, tolerance));
    BOOST_CHECK(xxy.derivative().col(0).isApprox(2.0 * vx * vy));
    BOOST_CHECK(xxy.derivative().col(1).isApprox(vx * vx));
}



BOOST_AUTO_TEST_CASE(Division)
{
    V vx(3);
    vx << 0.2, 1.2, 13.4;

    V vy(3);
    vy << 1.0, 2.2, 3.4;

    std::vector<V> vals{ vx, vy };
    std::vector<ADD> vars = ADD::variables(vals);

    const ADD x = vars[0];
    const ADD y = vars[1];

    const ADD xxBy = x * x / y;

    const double tolerance = 1e-14;
    BOOST_CHECK(xxBy.value().isApprox(vx * vx / vy, tolerance));
    BOOST_CHECK(xxBy.derivative().col(0).isApprox(2.0 * vx / vy));
    BOOST_CHECK(xxBy.derivative().col(1).isApprox(-vx * vx / (vy * vy)));
}
