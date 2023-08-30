#include <opm/io/eclipse/ESmry.hpp>

#include <string>
#include <cmath>
#include <iostream>

#define BOOST_TEST_MODULE TestTuningTSINIT
#include <boost/test/unit_test.hpp>

constexpr double TOLERANCE = 1.0e-10;

inline bool is_close(double a, double b, double tol = TOLERANCE) {
    return std::abs(a-b) <= tol;
}

BOOST_AUTO_TEST_CASE(CheckTSINITAndNEXTSTEP)
{
    std::string case_name("02_TUNING_TSINIT_NEXTSTEP");

    BOOST_TEST_MESSAGE("---------------------------------------------------------------------------");
    BOOST_TEST_MESSAGE("Checking TSINIT and NEXTSTEP, see file " + case_name + ".DATA");
    BOOST_TEST_MESSAGE("---------------------------------------------------------------------------");
    Opm::EclIO::ESmry smry(case_name, false);
    smry.loadData({"TIME"});
    const auto& time = smry.get("TIME");
        
    // First time step 1 day
    BOOST_CHECK_CLOSE(time[0], 1.0, TOLERANCE);

    for (size_t i=0; i<time.size(); ++i) {
        std::cout << "######################################" << std::endl;
        std::cout << time[i] << std::endl;
                
        // Max time step 3 days initially
        if (time[i] < 14.0) BOOST_CHECK(time[i+1] <= (time[i]+3.0));
        // No short next step 
        if (is_close(time[i], 14.0)) BOOST_CHECK(time[i+1] > 15.1);
        // Persistent NEXTSTEP=0.5
        if (is_close(time[i], 31.0)) BOOST_CHECK_CLOSE(time[i+1], 31.5, TOLERANCE);
        if (is_close(time[i], 45.0)) BOOST_CHECK_CLOSE(time[i+1], 45.5, TOLERANCE);
        // Non-persistent NEXTSTEP=1.0
        if (is_close(time[i], 60.0)) BOOST_CHECK_CLOSE(time[i+1], 61.0, TOLERANCE);
        if (is_close(time[i], 74.0)) BOOST_CHECK(time[i+1] > 75.1);  
        // TSINIT=0.5
        if (is_close(time[i], 91.0)) BOOST_CHECK_CLOSE(time[i+1], 91.5, TOLERANCE);
        if (is_close(time[i], 105.0)) BOOST_CHECK(time[i+1] > 105.6);
    }  
}

