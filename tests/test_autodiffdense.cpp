/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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

#define BOOST_TEST_MODULE AutoDiffDenseTest

#include <opm/autodiff/AutoDiffDense.hpp>

#include <boost/test/unit_test.hpp>


const double tolerance = 1.0e-14;
const int NumDerivs = 2;
typedef Opm::AutoDiffDense<double, NumDerivs> AdFW;
typedef AdFW::Derivative Derivative;


void checkCloseDerivs(const Derivative& d1, const Derivative& d2, const double tol)
{
    for (int dd : { 0, 1 }) {
        BOOST_CHECK_CLOSE(d1[dd], d2[dd], tol);
    }
}


BOOST_AUTO_TEST_CASE(Initialisation)
{
    AdFW a = AdFW::variable(0, 0.0);
    AdFW b = AdFW::variable(0, 1.0);

    BOOST_CHECK_CLOSE(a.val(), 0.0, tolerance);
    BOOST_CHECK_CLOSE(b.val(), 1.0, tolerance);
    checkCloseDerivs(a.der(), b.der(), tolerance);

    a = b;
    BOOST_CHECK_EQUAL(a.val(), b.val());
    checkCloseDerivs(a.der(), b.der(), tolerance);

    AdFW c = AdFW::variable(1, 1.0);
    BOOST_CHECK_EQUAL(c.val(), 1.0);
    checkCloseDerivs(c.der(), { 0.0, 1.0 }, tolerance);
}



BOOST_AUTO_TEST_CASE(Addition)
{
    AdFW a = AdFW::variable(0, 0.0);
    AdFW b = AdFW::variable(0, 1.0);

    AdFW two_a = a + a;
    BOOST_CHECK_CLOSE(two_a.val(), 2*a.val(), tolerance);
    checkCloseDerivs(two_a.der(), 2*a.der(), tolerance);

    double av = a.val();
    Derivative ad = a.der();
    a += b;
    BOOST_CHECK_CLOSE(a.val(), av + b.val(), tolerance);
    checkCloseDerivs(a.der(), ad + b.der(), tolerance);

    av = a.val();
    ad = a.der();
    a += 1;
    BOOST_CHECK_CLOSE(a.val(), av + 1, tolerance);
    checkCloseDerivs(a.der(), ad    , tolerance);

    AdFW bpo = b + 1;           // b plus one
    BOOST_CHECK_CLOSE(bpo.val(), b.val() + 1, tolerance);
    checkCloseDerivs(bpo.der(), b.der()    , tolerance);

    AdFW opb = 1 + b;           // one plus b
    BOOST_CHECK_CLOSE(opb.val(), b.val() + 1, tolerance);
    checkCloseDerivs(opb.der(), b.der()    , tolerance);
}



BOOST_AUTO_TEST_CASE(Subtraction)
{
    AdFW a = AdFW::variable(0, 0.0);
    AdFW b = AdFW::variable(0, 1.0);

    AdFW no_a = a - a;
    BOOST_CHECK_CLOSE(no_a.val(), 0.0, tolerance);
    checkCloseDerivs(no_a.der(), { 0.0, 0.0 }, tolerance);

    AdFW amb = a - b;
    BOOST_CHECK_CLOSE(amb.val(), a.val() - b.val(), tolerance);
    checkCloseDerivs(amb.der(), a.der() - b.der(), tolerance);

    double av = a.val();
    Derivative ad = a.der();
    a -= b;
    BOOST_CHECK_CLOSE(a.val(), av - b.val(), tolerance);
    checkCloseDerivs(a.der(), ad - b.der(), tolerance);

    av = a.val();
    ad = a.der();
    a -= 1;
    BOOST_CHECK_CLOSE(a.val(), av - 1, tolerance);
    checkCloseDerivs(a.der(), ad    , tolerance);

    AdFW bmo = b - 1;           // b minus one
    BOOST_CHECK_CLOSE(bmo.val(), b.val() - 1, tolerance);
    checkCloseDerivs(bmo.der(), b.der()    , tolerance);

    AdFW omb = 1 - b;           // one minus b
    BOOST_CHECK_CLOSE(omb.val(), 1 - b.val(), tolerance);
    checkCloseDerivs(omb.der(),   - b.der(), tolerance);
}


BOOST_AUTO_TEST_CASE(Multiplication)
{
    AdFW a = AdFW::variable(0, 0.0);
    AdFW b = AdFW::variable(0, 1.0);

    AdFW no_a = a * 0;
    BOOST_CHECK_CLOSE(no_a.val(), 0.0, tolerance);
    checkCloseDerivs(no_a.der(), { 0.0, 0.0 }, tolerance);

    AdFW atb = a * b;
    BOOST_CHECK_CLOSE(atb.val(), a.val() * b.val(), tolerance);
    checkCloseDerivs(atb.der(), a.der()*b.val() + a.val()*b.der(), tolerance);

    double av = a.val();
    Derivative ad = a.der();
    a *= b;
    BOOST_CHECK_CLOSE(a.val(), av * b.val(), tolerance);
    checkCloseDerivs(a.der(), ad*b.val() + av*b.der(), tolerance);

    av = a.val();
    ad = a.der();
    a *= 1;
    BOOST_CHECK_CLOSE(a.val(), av, tolerance);
    checkCloseDerivs(a.der(), ad, tolerance);

    AdFW bto = b * 1;           // b times one
    BOOST_CHECK_CLOSE(bto.val(), b.val(), tolerance);
    checkCloseDerivs(bto.der(), b.der(), tolerance);

    AdFW otb = 1 * b;           // one times b
    BOOST_CHECK_CLOSE(otb.val(), b.val(), tolerance);
    checkCloseDerivs(otb.der(), b.der(), tolerance);
}


BOOST_AUTO_TEST_CASE(Division)
{
    AdFW a = AdFW::variable(0, 10.0);
    AdFW b = AdFW::variable(0, 1.0);

    AdFW aob = a / b;
    BOOST_CHECK_CLOSE(aob.val(), a.val() * b.val(), tolerance);
    const Derivative res = ((a.der()*b.val() - a.val()*b.der()) /
                            (b.val() * b.val()));
    checkCloseDerivs(aob.der(), res, tolerance);

    double av = a.val();
    Derivative ad = a.der();
    a /= b;
    BOOST_CHECK_CLOSE(a.val(), av * b.val(), tolerance);
    checkCloseDerivs(a.der(), res, tolerance);

    av = a.val();
    ad = a.der();
    a /= 2;
    BOOST_CHECK_CLOSE(a.val(), av / 2, tolerance);
    checkCloseDerivs(a.der(), ad / 2, tolerance);

    AdFW bot = b / 2;           // b over two
    BOOST_CHECK_CLOSE(bot.val(), b.val() / 2, tolerance);
    checkCloseDerivs(bot.der(), b.der() / 2, tolerance);

    AdFW otb = 2 / b;           // two over b
    BOOST_CHECK_CLOSE(otb.val(), 2 / b.val(), tolerance);
    checkCloseDerivs(otb.der(), -2*b.der() / (b.val() * b.val()), tolerance);
}


BOOST_AUTO_TEST_CASE(Polynomial)
{
    const AdFW x = AdFW::variable(0, 1.234e-1);

    const AdFW p0 = x * x;
    BOOST_CHECK_CLOSE(p0.val(), x.val() * x.val(), tolerance);
    checkCloseDerivs(p0.der(), 2*x.val()*x.der(), tolerance);

    const AdFW p = 10*x*x - x/2.0 + 3.0;
    BOOST_CHECK_CLOSE(p.val(), 10  *x.val()*x.val() - x.val()/2.0 + 3.0, tolerance);
    checkCloseDerivs(p.der(), 10*2*x.val()*x.der() - x.der()/2.0      , tolerance);
}
