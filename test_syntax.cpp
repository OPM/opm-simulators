#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE SyntaxTest

#include "AutoDiff.hpp"

#include <cmath>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Initialisation)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a = AdFW::variable(0.0);
    AdFW b = AdFW::variable(1.0);

    BOOST_CHECK_CLOSE(a.val(), 0.0, atol);
    BOOST_CHECK_CLOSE(b.val(), 1.0, atol);

    BOOST_CHECK_CLOSE(a.der(), b.der(), atol);

    a = b;
    BOOST_CHECK_EQUAL(a.val(), b.val());
    BOOST_CHECK_EQUAL(a.der(), b.der());
}


BOOST_AUTO_TEST_CASE(Addition)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a = AdFW::variable(0.0);
    AdFW b = AdFW::variable(1.0);

    AdFW two_a = a + a;
    BOOST_CHECK_CLOSE(two_a.val(), 2*a.val(), atol);
    BOOST_CHECK_CLOSE(two_a.der(), 2*a.der(), atol);

    double av = a.val();
    double ad = a.der();
    a += b;
    BOOST_CHECK_CLOSE(a.val(), av + b.val(), atol);
    BOOST_CHECK_CLOSE(a.der(), ad + b.der(), atol);

    av = a.val();
    ad = a.der();
    a += 1;
    BOOST_CHECK_CLOSE(a.val(), av + 1, atol);
    BOOST_CHECK_CLOSE(a.der(), ad    , atol);

    AdFW bpo = b + 1;           // b plus one
    BOOST_CHECK_CLOSE(bpo.val(), b.val() + 1, atol);
    BOOST_CHECK_CLOSE(bpo.der(), b.der()    , atol);

    AdFW opb = 1 + b;           // one plus b
    BOOST_CHECK_CLOSE(opb.val(), b.val() + 1, atol);
    BOOST_CHECK_CLOSE(opb.der(), b.der()    , atol);
}


BOOST_AUTO_TEST_CASE(Subtraction)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a = AdFW::variable(0.0);
    AdFW b = AdFW::variable(1.0);

    AdFW no_a = a - a;
    BOOST_CHECK_CLOSE(no_a.val(), 0.0, atol);
    BOOST_CHECK_CLOSE(no_a.der(), 0.0, atol);

    AdFW amb = a - b;
    BOOST_CHECK_CLOSE(amb.val(), a.val() - b.val(), atol);
    BOOST_CHECK_CLOSE(amb.der(), a.der() - b.der(), atol);

    double av = a.val();
    double ad = a.der();
    a -= b;
    BOOST_CHECK_CLOSE(a.val(), av - b.val(), atol);
    BOOST_CHECK_CLOSE(a.der(), ad - b.der(), atol);

    av = a.val();
    ad = a.der();
    a -= 1;
    BOOST_CHECK_CLOSE(a.val(), av - 1, atol);
    BOOST_CHECK_CLOSE(a.der(), ad    , atol);

    AdFW bmo = b - 1;           // b minus one
    BOOST_CHECK_CLOSE(bmo.val(), b.val() - 1, atol);
    BOOST_CHECK_CLOSE(bmo.der(), b.der()    , atol);

    AdFW omb = 1 - b;           // one minus b
    BOOST_CHECK_CLOSE(omb.val(), 1 - b.val(), atol);
    BOOST_CHECK_CLOSE(omb.der(),   - b.der(), atol);
}


BOOST_AUTO_TEST_CASE(Multiplication)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a = AdFW::variable(0.0);
    AdFW b = AdFW::variable(1.0);

    AdFW no_a = a * 0;
    BOOST_CHECK_CLOSE(no_a.val(), 0.0, atol);
    BOOST_CHECK_CLOSE(no_a.der(), 0.0, atol);

    AdFW atb = a * b;
    BOOST_CHECK_CLOSE(atb.val(), a.val() * b.val(), atol);
    BOOST_CHECK_CLOSE(atb.der(), a.der()*b.val() + a.val()*b.der(), atol);

    double av = a.val();
    double ad = a.der();
    a *= b;
    BOOST_CHECK_CLOSE(a.val(), av * b.val(), atol);
    BOOST_CHECK_CLOSE(a.der(), ad*b.val() + av*b.der(), atol);

    av = a.val();
    ad = a.der();
    a *= 1;
    BOOST_CHECK_CLOSE(a.val(), av, atol);
    BOOST_CHECK_CLOSE(a.der(), ad, atol);

    AdFW bto = b * 1;           // b times one
    BOOST_CHECK_CLOSE(bto.val(), b.val(), atol);
    BOOST_CHECK_CLOSE(bto.der(), b.der(), atol);

    AdFW otb = 1 * b;           // one times b
    BOOST_CHECK_CLOSE(otb.val(), b.val(), atol);
    BOOST_CHECK_CLOSE(otb.der(), b.der(), atol);
}


BOOST_AUTO_TEST_CASE(Division)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a = AdFW::variable(10.0);
    AdFW b = AdFW::variable(1.0);

    AdFW aob = a / b;
    BOOST_CHECK_CLOSE(aob.val(), a.val() * b.val(), atol);
    const double res = ((a.der()*b.val() - a.val()*b.der()) /
                        (b.val() * b.val()));
    BOOST_CHECK_CLOSE(aob.der(), res, atol);

    double av = a.val();
    double ad = a.der();
    a /= b;
    BOOST_CHECK_CLOSE(a.val(), av * b.val(), atol);
    BOOST_CHECK_CLOSE(a.der(), res, atol);

    av = a.val();
    ad = a.der();
    a /= 2;
    BOOST_CHECK_CLOSE(a.val(), av / 2, atol);
    BOOST_CHECK_CLOSE(a.der(), ad / 2, atol);

    AdFW bot = b / 2;           // b over two
    BOOST_CHECK_CLOSE(bot.val(), b.val() / 2, atol);
    BOOST_CHECK_CLOSE(bot.der(), b.der() / 2, atol);

    AdFW otb = 2 / b;           // two over b
    BOOST_CHECK_CLOSE(otb.val(), 2 / b.val(), atol);
    BOOST_CHECK_CLOSE(otb.der(), -2*b.der() / (b.val() * b.val()), atol);
}


BOOST_AUTO_TEST_CASE(Polynomial)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    const AdFW x = AdFW::variable(1.234e-1);

    const AdFW p0 = x * x;
    BOOST_CHECK_CLOSE(p0.val(), x.val() * x.val(), atol);
    BOOST_CHECK_CLOSE(p0.der(), 2*x.val()*x.der(), atol);

    const AdFW p = 10*x*x - x/2.0 + 3.0;
    BOOST_CHECK_CLOSE(p.val(), 10  *x.val()*x.val() - x.val()/2.0 + 3.0, atol);
    BOOST_CHECK_CLOSE(p.der(), 10*2*x.val()*x.der() - x.der()/2.0      , atol);
}


BOOST_AUTO_TEST_CASE(Cosine)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    const AdFW x = AdFW::variable(3.14159265358979323846264338327950288);

    const AdFW f = std::cos(x);
    BOOST_CHECK_CLOSE(f.val(),   std::cos(x.val()), atol);
    BOOST_CHECK_CLOSE(f.der(), - std::sin(x.val()), atol);

    const AdFW p = 10*x*x - x/2.0 + 3.0;
    const AdFW g = std::cos(p);
    BOOST_CHECK_CLOSE(g.val(),   std::cos(p.val())        , atol);
    BOOST_CHECK_CLOSE(g.der(), - std::sin(p.val())*p.der(), atol);
}


BOOST_AUTO_TEST_CASE(SquareRoot)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    const AdFW x = AdFW::variable(1.234e-5);

    const AdFW x2 = x * x;
    const AdFW g  = std::cos(x2) + x;
    const AdFW f  = std::sqrt(g) - 1.2;

    BOOST_CHECK_CLOSE(g.val(),  std::cos(x2.val())          + x.val(), atol);
    BOOST_CHECK_CLOSE(g.der(), -std::sin(x2.val())*x2.der() + x.der(), atol);

    BOOST_CHECK_CLOSE(f.val(),  std::sqrt(g.val()) - 1.2, atol);
    BOOST_CHECK_CLOSE(f.der(), 1.0/(2.0 * std::sqrt(g.val())) * g.der(), atol);
}
