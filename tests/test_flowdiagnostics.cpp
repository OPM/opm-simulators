/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#if defined(HAVE_DYNAMIC_BOOST_TEST)
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing


#define BOOST_TEST_MODULE FlowDiagnosticsTests
#include <boost/test/unit_test.hpp>
#include <opm/core/flowdiagnostics/FlowDiagnostics.hpp>

const std::vector<double> pv(16, 18750.0);

const std::vector<double> ftof = {
    5.399999999999992e+04,
    1.139999999999997e+05,
    2.819999999999993e+05,
    8.220000000000012e+05,
    1.139999999999998e+05,
    1.774285714285711e+05,
    3.160150375939849e+05,
    8.156820645778908e+05,
    2.819999999999994e+05,
    3.160150375939841e+05,
    3.935500938781204e+05,
    7.765612369042073e+05,
    8.220000000000000e+05,
    8.156820645778894e+05,
    7.765612369042063e+05,
    8.218220225906991e+05
};

const std::vector<double> rtof = {
    8.218220225906990e+05,
    7.765612369042046e+05,
    8.156820645778881e+05,
    8.219999999999976e+05,
    7.765612369042051e+05,
    3.935500938781204e+05,
    3.160150375939846e+05,
    2.820000000000001e+05,
    8.156820645778885e+05,
    3.160150375939850e+05,
    1.774285714285714e+05,
    1.140000000000000e+05,
    8.219999999999980e+05,
    2.819999999999998e+05,
    1.140000000000000e+05,
    5.400000000000003e+04
};

const std::vector<double> F = {
    0,
    9.568875799840706e-02,
    1.913775159968141e-01,
    2.778231480508526e-01,
    3.642687801048911e-01,
    4.266515906731506e-01,
    4.890344012414101e-01,
    5.503847464649610e-01,
    6.117350916885119e-01,
    6.730854369120627e-01,
    7.344357821356134e-01,
    7.842099754925904e-01,
    8.339841688495674e-01,
    8.837583622065442e-01,
    9.335325555635212e-01,
    9.667662777817606e-01,
    1.000000000000000e+00
};

const std::vector<double> Phi = {
    0,
    6.250000000000000e-02,
    1.250000000000000e-01,
    1.875000000000000e-01,
    2.500000000000000e-01,
    3.125000000000000e-01,
    3.750000000000000e-01,
    4.375000000000000e-01,
    5.000000000000000e-01,
    5.625000000000000e-01,
    6.250000000000000e-01,
    6.875000000000000e-01,
    7.500000000000000e-01,
    8.125000000000000e-01,
    8.750000000000000e-01,
    9.375000000000000e-01,
    1.000000000000000e+00
};

const std::vector<double> Ev = {
    0,
    6.531592770912591e-01,
    6.531592770912593e-01,
    7.096322601771997e-01,
    7.096322601772002e-01,
    8.869254748464411e-01,
    8.869254748464422e-01,
    8.955406718746977e-01,
    8.955406718746983e-01,
    8.955406718746991e-01,
    8.955406718746991e-01,
    9.584612275378565e-01,
    9.584612275378565e-01,
    9.584612275378569e-01,
    9.584612275378566e-01,
    1.000000000000000e+00,
    1.000000000000000e+00
};

const std::vector<double> tD = {
    0,
    6.531592770912591e-01,
    6.531592770912593e-01,
    7.229977792392133e-01,
    7.229977792392139e-01,
    1.001878553253259e+00,
    1.001878553253261e+00,
    1.018739173712224e+00,
    1.018739173712226e+00,
    1.018739173712227e+00,
    1.018739173712227e+00,
    1.255670776053656e+00,
    1.255670776053656e+00,
    1.255670776053659e+00,
    1.255670776053656e+00,
    1.880619919417231e+00,
    1.880619919417231e+00
};

std::vector<double> wrong_length(7, 0.0);



using namespace Opm;


template <class C>
void compareCollections(const C& c1, const C& c2, const double tolerance = 1e-11)
{
    BOOST_REQUIRE(c1.size() == c2.size());
    auto c1it = c1.begin();
    auto c2it = c2.begin();
    for (; c1it != c1.end(); ++c1it, ++c2it) {
        BOOST_CHECK_CLOSE(*c1it, *c2it, tolerance);
    }
}


BOOST_AUTO_TEST_CASE(FandPhi)
{
    BOOST_CHECK_THROW(computeFandPhi(pv, ftof, wrong_length), std::runtime_error);
    auto FPhi = computeFandPhi(pv, ftof, rtof);
    compareCollections(FPhi.first, F);
    compareCollections(FPhi.second, Phi);
}




BOOST_AUTO_TEST_CASE(Lorenz)
{
    BOOST_CHECK_THROW(computeLorenz(F, wrong_length), std::runtime_error);
    const double Lc = computeLorenz(F, Phi);
    BOOST_CHECK_CLOSE(Lc, 1.645920738950826e-01, 1e-11);
}




BOOST_AUTO_TEST_CASE(Sweep)
{
    BOOST_CHECK_THROW(computeSweep(F, wrong_length), std::runtime_error);
    auto et = computeSweep(F, Phi);
    compareCollections(et.first, Ev);
    compareCollections(et.second, tD);
}
