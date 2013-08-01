#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE ScalarMultTest

#include <opm/autodiff/AutoDiff.hpp>

#include <cmath>
#include <boost/test/unit_test.hpp>

#include <iostream>

BOOST_AUTO_TEST_CASE(ScalarMultiplication)
{


    std::cout << "er her" << std::endl;
    

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
