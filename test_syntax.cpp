#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE SyntaxTest

#include "AutoDiff.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Initialisation)
{
    typedef AutoDiff::Forward<double> AdFW;

    const double atol = 1.0e-14;

    AdFW a(0.0), b(1.0);

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

    AdFW a(0.0), b(1.0);

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

    AdFW a(0.0), b(1.0);

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
