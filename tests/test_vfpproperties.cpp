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

#define BOOST_TEST_MODULE VFPTest

#include <algorithm>
#include <filesystem>
#include <memory>
#include <map>
#include <sstream>
#include <limits>
#include <vector>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/FileSystem.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>



const double max_d_tol = 1.0e-10;
const double sad_tol = 1.0e-8;










BOOST_AUTO_TEST_SUITE( HelperTests )

BOOST_AUTO_TEST_CASE(findInterpData)
{
    std::vector<double> values = {1, 5, 7, 9, 11, 15};
    double exact = 9.0;
    double interpolate = 6.0;
    double extrapolate_left = -1.0;
    double extrapolate_right = 19;
    double first = 1;
    double last = 15;

    Opm::detail::InterpData eval0 = Opm::detail::findInterpData(exact, values);
    Opm::detail::InterpData eval1 = Opm::detail::findInterpData(interpolate, values);
    Opm::detail::InterpData eval2 = Opm::detail::findInterpData(extrapolate_left, values);
    Opm::detail::InterpData eval3 = Opm::detail::findInterpData(extrapolate_right, values);
    Opm::detail::InterpData eval4 = Opm::detail::findInterpData(first, values);
    Opm::detail::InterpData eval5 = Opm::detail::findInterpData(last, values);

    BOOST_CHECK_EQUAL(eval0.ind_[0], 2);
    BOOST_CHECK_EQUAL(eval0.ind_[1], 3);
    BOOST_CHECK_EQUAL(eval0.factor_, 1.0);

    BOOST_CHECK_EQUAL(eval1.ind_[0], 1);
    BOOST_CHECK_EQUAL(eval1.ind_[1], 2);
    BOOST_CHECK_EQUAL(eval1.factor_, 0.5);

    BOOST_CHECK_EQUAL(eval2.ind_[0], 0);
    BOOST_CHECK_EQUAL(eval2.ind_[1], 1);
    BOOST_CHECK_EQUAL(eval2.factor_, -0.25);

    BOOST_CHECK_EQUAL(eval3.ind_[0], 4);
    BOOST_CHECK_EQUAL(eval3.ind_[1], 5);
    BOOST_CHECK_EQUAL(eval3.factor_, 2.0);

    BOOST_CHECK_EQUAL(eval4.ind_[0], 0);
    BOOST_CHECK_EQUAL(eval4.ind_[1], 1);
    BOOST_CHECK_EQUAL(eval4.factor_, 0.0);

    BOOST_CHECK_EQUAL(eval5.ind_[0], 4);
    BOOST_CHECK_EQUAL(eval5.ind_[1], 5);
    BOOST_CHECK_EQUAL(eval5.factor_, 1.0);
}

BOOST_AUTO_TEST_SUITE_END() // HelperTests


/**
 * Test fixture to set up axis etc.
 * All of our axes go from 0 to 1, but with a varying number of
 * values data is given at
 */
struct TrivialFixture {
    typedef Opm::detail::VFPEvaluation VFPEvaluation;

    TrivialFixture() : table_ids(1, 1),
            thp_axis{0.0, 1.0},
            wfr_axis{0.0, 0.5, 1.0},
            gfr_axis{0.0, 0.25, 0.5, 0.75, 1},
            alq_axis{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1},
            flo_axis{0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,
                0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1},
            nx(thp_axis.size()),
            ny(wfr_axis.size()),
            nz(gfr_axis.size()),
            nu(alq_axis.size()),
            nv(flo_axis.size()),
            data(nx*ny*nz*nu*nv)
    {
    }

    ~TrivialFixture() {

    }

    /**
     * Fills our interpolation data with zeros
     */
    inline void fillData(double value) {
        for (int i=0; i<nx; ++i) {
            for (int j=0; j<ny; ++j) {
                for (int k=0; k<nz; ++k) {
                    for (int l=0; l<nu; ++l) {
                        for (int m=0; m<nv; ++m) {
                            (*this)(i,j,k,l,m) = value;
                        }
                    }
                }
            }
        }
    }

    /**
     * Fills our interpolation data with an ND plane
     */
    inline void fillDataPlane() {
        for (int i=0; i<nx; ++i) {
            double x = i / static_cast<double>(nx-1);
            for (int j=0; j<ny; ++j) {
                double y = j / static_cast<double>(ny-1);
                for (int k=0; k<nz; ++k) {
                    double z = k / static_cast<double>(nz-1);
                    for (int l=0; l<nu; ++l) {
                        double u = l / static_cast<double>(nu-1);
                        for (int m=0; m<nv; ++m) {
                            double v = m / static_cast<double>(nv-1);
                            // table[thp_idx][wfr_idx][gfr_idx][alq_idx][flo_idx];
                            (*this)(i,j,k,l,m) = x + 2*y + 3*z + 4*u + 5*v;
                        }
                    }
                }
            }
        }
    }



    /**
     * Fills our interpolation data with "random" values
     */
    inline void fillDataRandom() {
        unsigned long randx = 42;
        static double max_val = static_cast<double>(std::numeric_limits<unsigned long>::max());
        for (int i=0; i<nx; ++i) {
            for (int j=0; j<ny; ++j) {
                for (int k=0; k<nz; ++k) {
                    for (int l=0; l<nu; ++l) {
                        for (int m=0; m<nv; ++m) {
                            (*this)(i,j,k,l,m) = randx / max_val;
                            randx = (randx*1103515245 + 12345);
                        }
                    }
                }
            }
        }
    }


    inline void initProperties() {
        table.reset(new Opm::VFPProdTable(1,
                                          1000.0,
                                          Opm::VFPProdTable::FLO_TYPE::FLO_OIL,
                                          Opm::VFPProdTable::WFR_TYPE::WFR_WOR,
                                          Opm::VFPProdTable::GFR_TYPE::GFR_GOR,
                                          Opm::VFPProdTable::ALQ_TYPE::ALQ_UNDEF,
                                          flo_axis,
                                          thp_axis,
                                          wfr_axis,
                                          gfr_axis,
                                          alq_axis,
                                          data));

        properties = std::make_shared<Opm::VFPProdProperties>();
        properties->addTable( *table );
    }

    double& operator()(size_t thp_idx, size_t wfr_idx, size_t gfr_idx, size_t alq_idx, size_t flo_idx) {
    return data[thp_idx*ny*nz*nu*nv + wfr_idx*nz*nu*nv + gfr_idx*nu*nv + alq_idx*nv + flo_idx];
}




    std::shared_ptr<Opm::VFPProdProperties> properties;
    std::shared_ptr<Opm::VFPProdTable> table;
    std::vector<int> table_ids;

private:
    const std::vector<double> thp_axis;
    const std::vector<double> wfr_axis;
    const std::vector<double> gfr_axis;
    const std::vector<double> alq_axis;
    const std::vector<double> flo_axis;
    int nx;
    int ny;
    int nz;
    int nu;
    int nv;
    std::vector<double> data;
};





//Set F to be our test suite fixture for our "trivial" tests
BOOST_FIXTURE_TEST_SUITE( TrivialTests, TrivialFixture )



/**
 * Test that we can generate some dummy zero-data,
 * interpolate using doubles as input, and compare against the analytic solution
 */
BOOST_AUTO_TEST_CASE(InterpolateZero)
{
    fillData(0.0);
    initProperties();

    //Check interpolation
    double sum = 0.0;
    int n=5;
    for (int i=1; i<n; ++i) {
        const double x = i / static_cast<double>(n-1);
        for (int j=1; j<n; ++j) {
            const double y = j / static_cast<double>(n-1);
            for (int k=1; k<n; ++k) {
                const double z = k / static_cast<double>(n-1);
                for (int l=0; l<n; ++l) {
                    const double u = l / static_cast<double>(n-1);
                    for (int m=0; m<n; ++m) {
                        const double v = m / static_cast<double>(n-1);

                        //Note order of arguments!
                        sum += properties->bhp(1, v, x, y, z, u, 0.0, 0.0, false);
                    }
                }
            }
        }
    }

    BOOST_CHECK_EQUAL(sum, 0.0);
}


/**
 * Test that we can generate some dummy one-data,
 * interpolate using doubles as input, and compare against the analytic solution
 */
BOOST_AUTO_TEST_CASE(InterpolateOne)
{
    fillData(1.0);
    initProperties();

    //Check interpolation
    double sum = 0.0;
    int n=5;
    for (int i=1; i<n; ++i) {
        const double x = i / static_cast<double>(n-1);
        for (int j=1; j<n; ++j) {
            const double y = j / static_cast<double>(n-1);
            for (int k=1; k<n; ++k) {
                const double z = k / static_cast<double>(n-1);
                for (int l=0; l<n; ++l) {
                    const double u = l / static_cast<double>(n-1);
                    for (int m=0; m<n; ++m) {
                        const double v = m / static_cast<double>(n-1);

                        //Note order of arguments!
                        const double value = properties->bhp(1, v, x, y, z, u, 0, 0, false);

                        sum += value;
                    }
                }
            }
        }
    }

    double reference = (n-1)*(n-1)*(n-1)*n*n;
    BOOST_CHECK_EQUAL(sum, reference);
}


///**
// * Test that we can generate some dummy data representing an ND plane,
// * interpolate using doubles as input, and compare against the analytic solution
// */
//BOOST_AUTO_TEST_CASE(InterpolatePlane)
//{
//    const int n=5;

//    fillDataPlane();
//    initProperties();

//    //Temps used to store reference and actual variables
//    double sad = 0.0;
//    double max_d = 0.0;

//    //Check interpolation
//    for (int i=0; i<=n; ++i) {
//        const double thp = i / static_cast<double>(n);
//        for (int j=1; j<=n; ++j) {
//            const double aqua = -j / static_cast<double>(n);
//            for (int k=1; k<=n; ++k) {
//                const double vapour = -k / static_cast<double>(n);
//                for (int l=0; l<=n; ++l) {
//                    const double alq = l / static_cast<double>(n);
//                    for (int m=1; m<=n; ++m) {
//                        const double liquid = -m / static_cast<double>(n);

//                        //Find values that should be in table
//                        double flo = Opm::detail::getFlo(aqua, liquid, vapour, table->getFloType());
//                        double wfr = Opm::detail::getWFR(aqua, liquid, vapour, table->getWFRType());
//                        double gfr = Opm::detail::getGFR(aqua, liquid, vapour, table->getGFRType());

//                        //Calculate reference
//                        double reference = thp + 2*wfr + 3*gfr+ 4*alq - 5*flo;

//                        //Calculate actual
//                        //Note order of arguments: id, aqua, liquid, vapour, thp, alq
//                        double actual = properties->bhp(1, aqua, liquid, vapour, thp, alq);


//                        double abs_diff = std::abs(actual - reference);
//                        max_d = std::max(max_d, abs_diff);
//                        sad = sad + abs_diff;
//                    }
//                }
//            }
//        }
//    }

//    BOOST_CHECK_SMALL(max_d, max_d_tol);
//    BOOST_CHECK_SMALL(sad, sad_tol);
//}



///**
// * Test that we can generate some dummy data representing an ND plane,
// * interpolate using doubles as input, and compare against the analytic solution
// */
//BOOST_AUTO_TEST_CASE(ExtrapolatePlane)
//{
//    fillDataPlane();
//    initProperties();

//    //Check linear extrapolation (i.e., using values of x, y, etc. outside our interpolant domain)
//    double sum = 0.0;
//    double reference_sum = 0.0;
//    double sad = 0.0; // Sum absolute difference
//    double max_d = 0.0; // Maximum difference
//    int n=1;
//    int o=5;
//    for (int i=0; i<=n+o; ++i) {
//        const double x = i / static_cast<double>(n);
//        for (int j=1; j<=n+o; ++j) {
//            const double aqua = -j / static_cast<double>(n);
//            for (int k=1; k<=n+o; ++k) {
//                const double vapour = -k / static_cast<double>(n);
//                for (int l=0; l<=n+o; ++l) {
//                    const double u = l / static_cast<double>(n);
//                    for (int m=1; m<=n+o; ++m) {
//                        const double liquid = -m / static_cast<double>(n);

//                        //Find values that should be in table
//                        double v = Opm::detail::getFlo(aqua, liquid, vapour, table->getFloType());
//                        double y = Opm::detail::getWFR(aqua, liquid, vapour, table->getWFRType());
//                        double z = Opm::detail::getGFR(aqua, liquid, vapour, table->getGFRType());

//                        double reference = x + 2*y + 3*z+ 4*u - 5*v;
//                        reference_sum += reference;

//                        //Note order of arguments! id, aqua, liquid, vapour, thp , alq
//                        double value = properties->bhp(1, aqua, liquid, vapour, x, u);
//                        sum += value;

//                        double abs_diff = std::abs(value - reference);

//                        sad += std::abs(abs_diff);
//                        max_d = std::max(max_d, abs_diff);
//                    }
//                }
//            }
//        }
//    }

//    BOOST_CHECK_CLOSE(sum, reference_sum, 0.0001);
//    BOOST_CHECK_SMALL(max_d, max_d_tol);
//    BOOST_CHECK_SMALL(sad, sad_tol);
//}





///**
// * Test that the partial derivatives are reasonable
// */
//BOOST_AUTO_TEST_CASE(PartialDerivatives)
//{
//    const int n=5;

//    fillDataPlane();
//    initProperties();

//    //Temps used to store reference and actual variables
//    VFPEvaluation sad;
//    VFPEvaluation max_d;

//    //Check interpolation
//    for (int i=0; i<=n; ++i) {
//        const double thp = i / static_cast<double>(n);
//        for (int j=1; j<=n; ++j) {
//            const double aqua = -j / static_cast<double>(n);
//            for (int k=1; k<=n; ++k) {
//                const double vapour = -k / static_cast<double>(n);
//                for (int l=0; l<=n; ++l) {
//                    const double alq = l / static_cast<double>(n);
//                    for (int m=1; m<=n; ++m) {
//                        const double liquid = -m / static_cast<double>(n);

//                        //Find values that should be in table
//                        double flo = Opm::detail::getFlo(aqua, liquid, vapour, table->getFloType());
//                        double wfr = Opm::detail::getWFR(aqua, liquid, vapour, table->getWFRType());
//                        double gfr = Opm::detail::getGFR(aqua, liquid, vapour, table->getGFRType());

//                        //Calculate reference
//                        VFPEvaluation reference;
//                        reference.value = thp + 2*wfr + 3*gfr+ 4*alq - 5*flo;
//                        reference.dthp = 1;
//                        reference.dwfr = 2;
//                        reference.dgfr = 3;
//                        reference.dalq = 4;
//                        reference.dflo = 5;

//                        //Calculate actual
//                        //Note order of arguments: id, aqua, liquid, vapour, thp, alq
//                        VFPEvaluation actual = Opm::detail::bhp(table.get(), aqua, liquid, vapour, thp, alq);

//                        VFPEvaluation abs_diff = actual - reference;
//                        abs_diff.value = std::abs(abs_diff.value);
//                        abs_diff.dthp = std::abs(abs_diff.dthp);
//                        abs_diff.dwfr = std::abs(abs_diff.dwfr);
//                        abs_diff.dgfr = std::abs(abs_diff.dgfr);
//                        abs_diff.dalq = std::abs(abs_diff.dalq);
//                        abs_diff.dflo = std::abs(abs_diff.dflo);

//                        max_d.value = std::max(max_d.value, abs_diff.value);
//                        max_d.dthp = std::max(max_d.dthp, abs_diff.dthp);
//                        max_d.dwfr = std::max(max_d.dwfr, abs_diff.dwfr);
//                        max_d.dgfr = std::max(max_d.dgfr, abs_diff.dgfr);
//                        max_d.dalq = std::max(max_d.dalq, abs_diff.dalq);
//                        max_d.dflo = std::max(max_d.dflo, abs_diff.dflo);

//                        sad = sad + abs_diff;
//                    }
//                }
//            }
//        }
//    }

//    BOOST_CHECK_SMALL(max_d.value, max_d_tol);
//    BOOST_CHECK_SMALL(max_d.dthp, max_d_tol);
//    BOOST_CHECK_SMALL(max_d.dwfr, max_d_tol);
//    BOOST_CHECK_SMALL(max_d.dgfr, max_d_tol);
//    BOOST_CHECK_SMALL(max_d.dalq, max_d_tol);
//    BOOST_CHECK_SMALL(max_d.dflo, max_d_tol);

//    BOOST_CHECK_SMALL(sad.value, sad_tol);
//    BOOST_CHECK_SMALL(sad.dthp, sad_tol);
//    BOOST_CHECK_SMALL(sad.dwfr, sad_tol);
//    BOOST_CHECK_SMALL(sad.dgfr, sad_tol);
//    BOOST_CHECK_SMALL(sad.dalq, sad_tol);
//    BOOST_CHECK_SMALL(sad.dflo, sad_tol);
//}




BOOST_AUTO_TEST_CASE(THPToBHPAndBackPlane)
{
    fillDataPlane();
    initProperties();

    double aqua = -0.5;
    double liquid = -0.9;
    double vapour = -0.1;
    double thp = 0.5;
    double alq = 32.9;

    double bhp_val = properties->bhp(1, aqua, liquid, vapour, thp, alq, 0, 0 , false);
    double thp_val = properties->thp(1, aqua, liquid, vapour, bhp_val, alq);

    BOOST_CHECK_CLOSE(thp_val, thp, max_d_tol);
}



BOOST_AUTO_TEST_CASE(THPToBHPAndBackNonTrivial)
{
    fillDataRandom();
    initProperties();

    double aqua = -0.5;
    double liquid = -0.9;
    double vapour = -0.1;
    double thp = 0.5;
    double alq = 32.9;

    double bhp_val = properties->bhp(1, aqua, liquid, vapour, thp, alq, 0, 0, false);
    double thp_val = properties->thp(1, aqua, liquid, vapour, bhp_val, alq);

    BOOST_CHECK_CLOSE(thp_val, thp, max_d_tol);
}

BOOST_AUTO_TEST_CASE(ExplicitBhpLookup)
{
    fillDataRandom();
    initProperties();

    double aqua = -0.5;
    double liquid = -0.9;
    double vapour = -0.1;
    double thp = 0.5;
    double alq = 32.9;

    double wfr = Opm::detail::getWFR(*table, aqua, liquid, vapour);
    double gfr = Opm::detail::getGFR(*table, aqua, liquid, vapour);

    double bhp_val = properties->bhp(1, aqua, liquid, vapour, thp, alq, wfr, gfr, false);
    double bhp_val_explicit = properties->bhp(1, 2*aqua, liquid, 3*vapour, thp, alq, wfr, gfr, true);
    BOOST_CHECK_CLOSE(bhp_val, bhp_val_explicit, max_d_tol);
}


BOOST_AUTO_TEST_SUITE_END() // Trivial tests


BOOST_AUTO_TEST_SUITE(IntegrationTests)

extern const double reference[];


/**
 * Uses a VFP table that should be linear in THP vs BHP,
 * so that
 *   bhp(aqua, liquid, vapour, thp, alq) == 0.0 for thp == 0.0,
 * and
 *   bhp(aqua, liquid, vapour, thp, alq) == 1.0 for thp == 1.0,
 */
BOOST_AUTO_TEST_CASE(ParseInterpolateLine)
{
    std::string table_str = "\
-- VFP table that is basically the identity: BHP == THP   \n\
-- In our case, we simply degenerate all axes except the THP axis \n\
-- and set that bhp(flo, thp, wfr, gfr, alq) = 0 for thp = 0, and 1 for thp = 1. \n\
-- The value of flo, wfr, gfr, alq, should all be irrelevant. \n\
VFPPROD \n\
-- table_num, datum_depth, flo, wfr, gfr, pressure, alq, unit, table_vals \n\
42 7.0E+03 LIQ WCT GOR THP ' ' FIELD BHP / \n\
1.0 / flo axis \n\
0.0 1.0 / THP axis \n\
0.0 / WFR axis \n\
0.0 / GFR axis \n\
0.0 / ALQ axis \n\
-- Table itself: thp_idx wfr_idx gfr_idx alq_idx <vals> \n\
1 1 1 1 0.0 / \n\
2 1 1 1 1.0 / \n\
";

    auto units = Opm::UnitSystem::newFIELD();
    Opm::Parser parser;
    auto deck = parser.parseString(table_str);
    bool gaslift_active = false;
    Opm::VFPProdTable table(deck["VFPPROD"].front(), gaslift_active, units);
    Opm::VFPProdProperties properties;
    properties.addTable( table );

    const int n = 5; //Number of points to check per axis
    double bhp_sad = 0.0; //Sum of absolute difference
    double bhp_max_d = 0.0; //Maximum difference
    double thp_sad = 0.0;
    double thp_max_d = 0.0;
    for (int w=0; w<n; ++w) { //water
        for (int o=0; o<n; ++o) { //oil
            for (int g=0; g<n; ++g) { //gas
                for (int t=0; t<n; ++t) { //thp
                    for (int a=0; a<n; ++a) { //alq
                        double aqua = w * 52.3;
                        double liquid = o * 9.9;
                        double vapour = g * 0.1;
                        double thp = t * 456.78;
                        double alq = a * 42.24;

                        double bhp_interp = properties.bhp(42, aqua, liquid, vapour, thp, alq, 0, 0, false);
                        double bhp_ref = thp;
                        double thp_interp = properties.thp(42, aqua, liquid, vapour, bhp_ref, alq);
                        double thp_ref = thp;

                        double bhp_diff = std::abs(bhp_interp - bhp_ref);
                        bhp_sad += bhp_diff;
                        bhp_max_d = std::max(bhp_diff, bhp_max_d);

                        double thp_diff = std::abs(thp_interp - thp_ref);
                        thp_sad += thp_diff;
                        thp_max_d = std::max(thp_diff, thp_max_d);
                    }
                }
            }
        }
    }

    BOOST_CHECK_SMALL(bhp_max_d, max_d_tol);
    BOOST_CHECK_SMALL(bhp_sad, sad_tol);

    BOOST_CHECK_SMALL(thp_max_d, max_d_tol);
    BOOST_CHECK_SMALL(thp_sad, sad_tol);
}

/**
 * Tests that we can actually parse some input data, and interpolate within that space
 */
BOOST_AUTO_TEST_CASE(ParseInterpolateRealisticVFPPROD)
{
    auto units = Opm::UnitSystem::newMETRIC();

    Opm::Parser parser;
    std::filesystem::path file("VFPPROD2");

    auto deck = parser.parseFile(file.string());

    BOOST_REQUIRE(deck.hasKeyword("VFPPROD"));
    BOOST_CHECK_EQUAL(deck.count("VFPPROD"), 1);

    bool gaslift_active = false;
    Opm::VFPProdTable table(deck["VFPPROD"].front(), gaslift_active, units);
    Opm::VFPProdProperties properties;
    properties.addTable(table);

    //Do some rudimentary testing
    //Get the BHP as a function of rate, thp, wfr, gfr, alq
    double liq[] = {100, 2942.8571428571427, 5785.7142857142853, 8628.5714285714294, 11471.428571428571, 14314.285714285714, 17157.142857142859, 20000};
    double gor[] = {90, 1505.7142857142858, 2921.4285714285716, 4337.1428571428569, 5752.8571428571431, 7168.5714285714284, 8584.2857142857138, 10000};
    double wct[] = {0, 0.14285714285714285, 0.2857142857142857, 0.42857142857142855, 0.5714285714285714, 0.7142857142857143, 0.8571428571428571, 1};
    double thp[] = {16.010000000000002, 22.438571428571429, 28.867142857142859, 35.295714285714283, 41.724285714285713, 48.152857142857144, 54.581428571428575, 61.009999999999998};
    int n = sizeof(liq) / sizeof(liq[0]);

    int i = 0;
    double sad = 0.0; //Sum of absolute difference
    double max_d = 0.0; //Maximum difference
    for (int t=0; t<n; ++t) {
        for (int w=0; w<n; ++w) {
            for (int g=0; g<n; ++g) {
                //for (unsigned int a=0; a<n; ++a) { //n==1, skip this loop
                    for (int f=0; f<n; ++f) {
                        //Liq given as SM3/day => convert to SM3/second
                        double f_i = -liq[f]*1.1574074074074073e-05;

                        //THP given as BARSA => convert to Pascal
                        double t_i = thp[t]*100000.0;

                        //WCT given as fraction, SM3/SM3
                        double w_i = wct[w];

                        //GOR given as SM3 / SM3
                        double g_i = gor[g];

                        //ALQ unit not relevant in this case
                        double a_i = 0.0;

                        //Now reconstruct possible aqua, liquid,
                        //and vapour phase compositions that satisfy the above

                        //liq[f] = aqua + liquid
                        //wct[w] = aqua / (aqua + liquid)
                        //gor[g] = vapour / liquid
                        //
                        // aqua = wct[w] * liq[f]
                        // liquid = liq[f] - aqua
                        // vapour = gor[g] * liquid

                        double aqua = w_i * f_i;
                        double liquid = f_i - aqua;
                        double vapour = g_i * liquid;

                        if ((aqua + liquid) == 0.0 || liquid == 0.0) {
                            //FIXME: This skips some corner cases, when
                            //getWFR(...) and getGFR(...) are infinite
                        }
                        else {
                            //Value given as pascal, convert to barsa for comparison with reference
                            double value_i = properties.bhp(32, aqua, liquid, vapour, t_i, a_i, 0, 0, false) * 10.0e-6;

                            double abs_diff = std::abs(value_i - reference[i]);
                            sad += abs_diff;
                            max_d = std::max(max_d, abs_diff);
                        }
                        ++i;
                    }
                //}
            }
        }
    }

    BOOST_CHECK_SMALL(max_d, max_d_tol);
    BOOST_CHECK_SMALL(sad, sad_tol);
}

/**
 * Reference computed using MATLAB with the input above.
 */
const double reference[] = {
        44.850000000000001,  30.36904761904762, 40.238928571428566, 53.221714285714292, 67.885571428571438, 83.100857142857151, 99.021000000000015,             115.14, 27.986149999999999, 76.394680612244898, 146.58456352040815, 227.73038755102044, 323.29702921768705, 433.02474183673473, 519.09348404081641, 596.70497857142868, 25.483874999999998, 135.96136479591837, 276.35715114795914, 452.68607295918366, 676.76830807823114, 945.81213010204078,  1108.470853367347, 1234.3462321428572,            22.9816, 195.52804897959183, 406.12973877551008, 677.64175836734682, 1030.2395869387753, 1458.5995183673467, 1697.8482226938775, 1871.9874857142856,           21.93648, 260.13428359183678, 566.16254540816328, 972.26681151020421, 1344.8295676190476, 1774.9372779591836, 2013.0210626204084, 2185.0791771428571, 22.174319999999998,  329.1776176326531, 752.83812551020401,  1328.232636244898, 1625.1862990476188, 1918.3098995918363, 2086.7690417469385, 2212.4194514285709, 22.412159999999997, 398.22095167346936, 939.51370561224473,  1684.198460979592, 1905.5430304761903, 2061.6825212244894, 2160.5170208734694, 2239.7597257142857, 22.649999999999999, 467.26428571428573, 1126.1892857142857, 2040.1642857142861, 2185.8997619047618, 2205.0551428571425, 2234.2649999999999, 2267.0999999999999, 45.935714285714283, 31.165795918367351, 39.937755102040811,  52.72722448979593, 66.603455782312921, 81.156040816326524, 96.557926530612264, 112.18714285714286, 29.541501020408166,  68.10885813411079, 127.61792838921284, 195.14191513119539, 273.51727814382895, 362.26484685131197, 464.80925879300298, 571.05809489795934, 26.418334183673469, 118.58690218658893, 235.65489568148689, 377.46362485422748, 554.83594890670554, 764.83203053935858, 1015.1273887900876, 1276.0984719387757, 23.295167346938779, 169.06494623906704, 343.69186297376086,  559.7853345772595, 836.15461966958196, 1167.3992142274051, 1565.4455187871722, 1981.1388489795918, 21.823392653061227, 222.83265130417885, 472.45379858600586, 795.44837871720142, 1174.8842335704569, 1543.4152790437317, 1945.2266863183679, 2358.5925306122449, 21.805595102040819, 279.49675619825075, 619.46314463556848,  1078.076089212828, 1564.1616228338191, 1896.0542676618074, 2174.8576629877552, 2447.6207346938772, 21.787797551020411, 336.16086109232265, 766.47249068513111, 1360.7037997084549, 1953.4390120971816, 2248.6932562798829, 2404.4886396571433,   2536.64893877551, 21.770000000000003,  392.8249659863946, 913.48183673469384,  1643.331510204082, 2342.7164013605443, 2601.3322448979588, 2634.1196163265308, 2625.6771428571428,  47.03857142857143, 32.294666666666664, 39.655204081632654, 52.333755102040818, 66.056870748299303, 79.655795918367346, 94.904048979591849, 110.69571428571429,   31.4618112244898, 60.027269698736646, 109.31101308309036, 164.18323034985428, 226.53348555879492, 296.57316871720116, 377.61696073469392,  461.6498928571429, 27.707640306122446, 101.52058415937805, 197.02186251822158, 308.25257120991256, 443.23278785228376, 602.26389059766768,  794.2130148979594, 995.03794642857144, 23.953469387755103,  143.0138986200194, 284.73271195335275, 452.32191206997084, 659.93209014577235, 907.95461247813398, 1210.8090690612244, 1528.4259999999999, 22.096767346938773, 186.56175919922259,  385.1060212390671, 631.63031250145787,  943.5575367035957, 1305.3664132711372, 1664.3889660991258, 2019.3338204081635, 21.910702040816325, 231.91855601943632,  496.6280617784256, 841.96513350437328, 1286.1084689135077, 1783.5345340174924, 2150.5314876034986, 2472.8396897959183, 21.724636734693878, 277.27535283965011, 608.15010231778422, 1052.2999545072887, 1628.6594011234206, 2261.7026547638479, 2636.6740091078714, 2926.3455591836728,  21.53857142857143,  322.6321496598639, 719.67214285714283, 1262.6347755102042,  1971.210333333333, 2739.8707755102037, 3122.8165306122451, 3379.8514285714286, 48.172142857142859,   34.1152925170068, 39.812755102040818,   51.5600612244898, 65.615465986394554, 79.399040816326519, 94.464034693877551, 109.96857142857144, 33.996096428571427,  52.22108403790088, 92.038355357142862,  135.3683000145773, 183.62870803692908, 237.10243908163267, 298.29584687172019, 361.59696479591844, 29.556062499999996, 84.769836188046639, 160.74966900510205, 246.17866242711372, 346.31590522351797, 462.23331913265304, 600.79137528425667, 745.49898852040815, 25.116028571428572,  117.3185883381924, 229.46098265306119, 356.98902483965014, 509.00310241010675, 687.36419918367335, 903.28690369679305, 1129.4010122448979, 22.869793469387755, 151.34435734499516, 305.45451635568514, 488.96920695043741, 713.93443163848406, 980.47258414577254, 1273.7336899620996,  1572.743070204082, 22.555100408163263, 186.67057383090378, 387.85972008746353, 639.58847130029164, 956.05983537803672, 1333.4321309271136, 1704.0085286822157, 2068.4194277551019, 22.240407346938778, 221.99679031681239, 470.26492381924197, 790.20773565014576, 1198.1852391175898, 1686.3916777084546, 2134.2833674023323, 2564.0957853061223, 21.925714285714285, 257.32300680272107, 552.67012755102041, 940.82700000000011, 1440.3106428571427, 2039.3512244897956, 2564.5582061224495, 3059.7721428571431, 49.354285714285716, 36.632401360544222, 40.749489795918372, 51.042102040816339, 64.626265306122448, 78.587408163265295, 93.677473469387763, 109.11500000000001, 37.089295408163267, 45.091101822157441, 75.658611279154542, 108.27292311953356,  143.6558939091351, 182.14881000000003, 225.54025218513127, 270.28158724489799, 32.070335459183667, 68.480038083090392, 125.77860782616619, 188.00535940233243, 258.06514529276001, 337.00742500000001, 429.65019986516046, 526.05641198979606, 27.051375510204085, 91.868974344023314, 175.89860437317788, 267.73779568513123, 372.47439667638486, 491.86603999999994, 633.76014754518962,  781.8312367346939, 24.361997551020409, 116.59068105733722, 229.69027080174931,  358.0967598367348, 509.27644487949482, 684.21993831486884, 893.97807862099171,  1113.227684489796, 23.723712653061227, 142.48583272303208, 286.71467883381922, 457.81190791836747, 665.79435214188538, 909.58676159766753,  1203.596587107872, 1511.2055991836733, 23.085427755102046, 168.38098438872692, 343.73908686588919, 557.52705600000002, 822.31225940427601, 1134.9535848804665, 1513.2150955947525, 1909.1835138775511, 22.447142857142858, 194.27613605442178, 400.76349489795916, 657.24220408163274, 978.83016666666674, 1360.3204081632653, 1822.8336040816332, 2307.1614285714286, 50.617142857142852, 39.906816326530617, 42.617091836734701, 51.005673469387752, 63.234959183673467, 77.057918367346929, 92.260159183673466, 107.80571428571429, 40.776480612244896, 38.911281467444127, 60.277890889212827, 82.987050962099119, 106.68821241010689, 131.81003072886296, 159.38756669970846, 187.64996326530607, 35.347721938775507,  52.83015575801749, 92.067141217201154, 133.48505947521866, 177.94073199708453, 226.06629431486877, 280.09702088921279, 335.75996938775506, 29.918963265306125, 66.749030048590868,  123.8563915451895, 183.98306798833818, 249.19325158406212, 320.32255790087459,  400.8064750787172, 483.86997551020397, 26.750765714285716, 82.291052559766769, 157.43245227405248, 238.17180960932944, 328.07747020019428,  428.1440934110787,  543.3114419311953, 662.60384408163259, 25.572891428571427, 99.262184699708428, 192.58171988338188, 295.61007715451888,  413.6810594985422, 547.90924594752164, 705.00638577725942, 868.30065795918335, 24.395017142857142, 116.23331683965013, 227.73098749271128, 353.04834469970837, 499.28464879689005, 667.67439848396486, 866.70132962332332, 1073.9974718367343, 23.217142857142857,  133.2044489795918, 262.88025510204073, 410.48661224489791, 584.88823809523797, 787.43955102040798, 1028.3962734693878, 1279.6942857142853, 52.738571428571426, 45.757020408163264, 48.553826530612241, 55.025673469387762, 64.719925170068024, 76.397102040816335, 89.947404081632669, 103.97000000000001, 46.110853061224489, 34.359456734693879, 50.171643075801747, 68.670724868804683, 88.639938770651099, 108.26501833819242, 130.09957911953356, 152.66118265306122, 41.898091836734693, 39.401489455782318, 62.716569059766769, 88.206411880466504, 115.47130470116619, 143.54486392128285, 174.55314348396507, 206.44164030612248, 37.685330612244897, 44.443522176870758, 75.261495043731799, 107.74209889212833, 142.30267063168128, 178.82470950437326, 219.00670784839659, 260.22209795918377, 34.710051428571418, 50.933008283770675, 90.002424650145798, 131.07079234985429, 175.23425387366382, 222.83772099125378, 275.84280280000019, 330.33638285714301,           32.82432, 58.696912551992227, 106.67683752186591, 157.73905884548111, 213.53680643731784, 274.53989562682227, 343.58116377142875, 414.83187428571443, 30.938588571428571, 66.460816820213807, 123.35125039358604, 184.40732534110799, 251.83935900097188, 326.24207026239083, 411.31952474285731, 499.32736571428597, 29.052857142857142, 74.224721088435388, 140.02566326530615, 211.07559183673482, 290.14191156462596, 377.94424489795932, 479.05788571428604,  583.8228571428574, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 55.719999999999999, 56.041142857142859, 58.831428571428567,             63.686, 70.667999999999992, 79.513714285714286, 90.825885714285718,              102.8, 51.667142857142856, 38.150993197278915, 46.123418367346943, 58.862530612244903, 72.392680272108834, 86.932122448979598, 102.44057142857145, 118.19000000000001, 36.737085714285712, 78.681847376093302, 147.68909836005832, 228.42912478134116, 323.77457259475216, 433.46749090379006, 519.45052589795932, 596.94770102040832, 33.537571428571425, 137.47412621477164, 277.03675501093295, 453.19456997084558,  677.1337764820214, 946.44425262390666, 1109.0459455612245, 1234.7325484693879, 30.338057142857142, 196.26640505344994, 406.38441166180746, 677.96001516034994, 1030.4929803692903, 1459.4210143440232, 1698.6413652244896, 1872.5173959183671, 28.767510612244898, 260.47091397084557, 566.23613620991262, 972.57374535276995, 1345.0269787988339, 1775.7268427288632, 2006.8385977830908, 2169.4720379591836, 28.631197551020406,  329.4406501302235, 752.94533230320701, 1328.6857758134111, 1625.3771151175897, 1918.8726842682213, 2067.2772747125364, 2166.3408824489798, 28.494884489795918, 398.41038628960149, 939.65452839650141, 1684.7978062740526, 1905.7272514363458, 2062.0185258075799, 2127.7159516419829, 2163.2097269387759, 28.358571428571427, 467.38012244897953, 1126.3637244897959, 2040.9098367346942, 2186.0773877551019, 2205.1643673469389, 2188.1546285714289, 2160.0785714285716, 52.691836734693887, 39.207446064139951, 46.021209912536449, 58.635965014577266,  72.00586200194364, 85.730583090379014, 100.84179825072889, 116.37122448979594, 38.363955830903791, 70.753435464389852, 128.93646681591005, 195.98940140774684, 274.10951577259482, 362.76185536026657, 465.13747088421508, 571.18281734693892, 34.686060131195333, 120.36070035054841, 236.48560342044985, 378.05111293211166, 555.29757485422749, 765.42883865056228,  1015.365437422949, 1275.8005025510206, 31.008164431486883,  169.9679652367069, 344.03474002498956, 560.11282445647657, 836.48563393586016, 1168.0958219408578, 1565.5934039616827, 1980.4181877551023, 29.111674693877553, 223.27484935943363, 472.57878237192847,  795.7818930179094, 1175.6237736173819, 1541.1590488463144, 1943.5790696826327, 2358.9462570262394, 28.783633469387752, 279.83908324642516, 619.60790301957513,  1078.649848989588, 1565.7843751553519, 1888.1574883465219, 2169.9130701674308, 2450.3669604664724, 28.455592244897957, 336.40331713341664, 766.63702366722202, 1361.5178049612664, 1955.9449766933221, 2235.1559278467303, 2396.2470706522286, 2541.7876639067053, 28.127551020408163, 392.96755102040822, 913.66614431486892, 1644.3857609329452, 2346.1055782312924, 2582.1543673469391, 2622.5810711370268, 2633.2083673469388, 53.735102040816329, 40.592748299319723, 46.133615160349862, 58.163556851311952, 72.280184645286695, 85.457813411078718, 100.37657376093297, 115.92224489795919, 40.280739795918372, 63.127442950159661, 110.94077505206162, 165.19745674302376, 227.30377448563098, 297.13873097042898, 378.09813556559777, 462.08359023323618,  36.09679846938775,  103.6278374600861, 198.06740243648477, 308.87415564348191, 443.82413957378867, 602.75549106622248, 794.74167538629752, 995.65453097667637, 31.912857142857145, 144.12823197001251, 285.19402982090799, 452.55085454394003,  660.3445046619463, 908.37225116201591, 1211.3852152069971, 1529.2254717201167, 29.754170670553933, 187.14265184839655,  385.2952358788005, 631.76578637234491, 943.91654054171863, 1305.7580015160349, 1665.1022951680138, 2020.4748324198254,  29.37863078717201, 232.37055895821186,  496.8199798958766, 842.27054757184521, 1286.5245825833676, 1783.9422556268219, 2151.4609740381511,   2474.46192909621, 29.003090903790088, 277.59846606802722, 608.34472391295299, 1052.7753087713454, 1629.1326246250173,  2262.126509737609, 2637.8196529082879, 2928.4490257725947, 28.627551020408163, 322.82637317784258, 719.86946793002915, 1263.2800699708457, 1971.7406666666666, 2740.3107638483962, 3124.1783317784261, 3382.4361224489794, 54.803367346938785, 42.551428571428573, 46.797383381924206, 57.572142857142865, 71.539344995140908, 85.460104956268225, 100.61364635568516, 116.17724489795918, 42.693486005830906, 55.870048124392618, 94.186250687213658, 136.75799064139943, 184.63798636817992, 237.87313567055395, 298.93549172469818, 362.13140255102047, 38.033330174927109, 87.307938872691949, 162.15424984381508, 247.05119118075808, 346.97676925239477, 462.82645914723031, 601.36159565493563, 746.05213392857161, 33.373174344023326, 118.74582962099126,  230.1222490004165, 357.34439172011662, 509.31555213660971, 687.77978262390661, 903.78769958517296, 1129.9728653061225, 30.866328221574342, 152.16257690851037, 305.73657149104542, 489.11916252811346, 714.15369324309313, 980.85748641649332, 1274.2146453297796, 1573.2922929154522, 30.255375276967932, 187.32161913758156,  388.0831407122032, 639.80749805081223, 956.41062413161171, 1333.9156692728029,  1704.513299996335, 2068.9096238483962, 29.644422332361518, 222.48066136665278,  470.4297099333611, 790.49583357351116, 1198.6675550201303, 1686.9738521291129, 2134.8119546628905, 2564.5269547813414, 29.033469387755105, 257.63970359572403, 552.77627915451899,  941.1841690962101, 1440.9244859086491, 2040.0320349854228, 2565.1106093294461, 3060.1442857142861, 55.989591836734704, 45.084888241010695, 48.303724489795918, 57.560102040816339, 70.634364431486887, 84.610454810495639, 99.736604081632677, 115.18387755102043, 45.556975072886303,  49.60227519644593, 78.588905414410675, 110.30281877551025, 145.16576165798975, 183.33150930653898, 226.51136111141199, 271.06563177842565, 40.604647594752187, 71.693961432736359, 127.73951301541024,  189.3183402332362, 259.01403350426904, 337.82029577780099, 430.34326663317381, 526.61954008746363, 35.652320116618078, 93.785647669026815, 176.89012061640983, 268.33386169096218, 372.86230535054835, 492.30908224906295, 634.17517215493558,  782.1734483965015, 32.815774723032078, 117.78430685436626, 230.16663228862978, 358.34467954185766, 509.46220383867848, 684.55699630820516, 894.32220464506509, 1113.5187297959187, 31.842081107871721,  143.4619706855477, 287.07581851311954, 458.03636556434833, 666.09389504373178, 910.05013845231167, 1204.0521952249896, 1511.5949355102041, 30.868387492711371, 169.13963451672913,  343.9850047376093,   557.728051586839, 822.72558624878536, 1135.5432805964183, 1513.7821858049149, 1909.6711412244899,  29.89469387755102, 194.81729834791062, 400.89419096209917, 657.41973760932967, 979.35727745383872, 1361.0364227405248, 1823.5121763848401,  2307.747346938776, 57.334489795918373, 48.238981535471332, 50.777470845481048, 58.320390670553934, 69.809903790087461, 82.942221574344018, 97.654282798833833, 112.77306122448979, 48.897101749271144, 44.416583269471055, 64.618026504581422, 85.936928983756772, 109.01430224350965, 133.66214291128694, 160.91838596334861,  188.8981380466472, 43.893996355685132, 56.874853161876992, 95.006877954498108, 135.43520588296542, 179.44765501180061, 227.24285694502282, 281.04819719491877, 336.51170517492704, 38.890890962099128, 69.333123054282922, 125.39572940441479, 184.93348278217405, 249.88100778009158, 320.82357097875877, 401.17800842648887, 484.12527230320688, 35.776611953352777, 84.039488686380665, 158.22839931278631, 238.59834481632652, 328.36715979453004, 428.34693042065805, 543.46404571295295, 662.71498524781327, 34.325360349854229, 100.72520237012353, 193.21274239900035, 295.93169926530607, 413.94344245453271, 548.14616741357759, 705.25765204964591, 868.57781329446027, 32.874108746355681, 117.41091605386642, 228.19708548521447, 353.26505371428573, 499.51972511453545, 667.94540440649712, 867.05125838633887, 1074.4406413411075, 31.422857142857143, 134.09662973760931, 263.18142857142851, 410.59840816326533, 585.09600777453829, 787.74464139941654, 1028.8448647230321, 1280.3034693877548, 59.545714285714283, 53.636513119533532, 56.402507288629735, 62.446746355685136, 71.835520894071919, 83.065924198250727, 96.259005830903803, 109.96795918367349, 53.523951603498539, 40.934415568513117, 55.714955232194924, 72.914401411911712, 92.093786979036523, 112.61551650562268, 134.70237188796338, 157.20137346938779, 49.927170553935859, 44.844994023323622, 66.980597107975854, 91.361255549770945, 117.97211751006529, 146.45969159725115, 177.54288480945445, 209.32898469387757,  46.33038950437318, 48.755572478134127, 78.246238983756783, 109.80810968763019, 143.85044804109404, 180.30386668887968, 220.38339773094555, 261.45659591836744, 43.671299708454811, 54.491415689018481, 92.188539046230773, 132.48950584089971, 176.22964549354447, 223.57835120366525, 276.41757687830091,  330.7628246064142, 41.837805247813407, 61.834323125364449, 108.48751728446484, 158.89922741191177, 214.33256638345139, 275.15580264889644,  344.0681908051647, 415.19405994169114,  40.00431078717201, 69.177230561710417, 124.78649552269891, 185.30894898292388, 252.43548727335843, 326.73325409412757, 411.71880473202862, 499.62529527696813, 38.170816326530613, 76.520137998056384, 141.08547376093301, 211.71867055393597, 290.53840816326544, 378.31070553935876, 479.36941865889247, 584.05653061224518, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 62.158571428571435, 62.474680272108841, 65.263622448979589, 70.114571428571438, 77.096571428571423, 85.942285714285717, 97.254457142857149, 109.22857142857143, 58.262857142857143, 45.921564625850344, 52.308367346938766, 65.045142857142849, 77.828945578231298,             91.372, 106.60894285714286, 122.31857142857143, 44.968769387755103, 81.667181117589905, 149.25939165451894, 229.46841250728866, 324.56137429543242, 434.06000527696801, 519.98789020699724, 597.47068163265317, 41.390045918367349, 139.49805512633625, 278.04769670189501, 453.85388287172015, 677.71850182215735, 946.90760473760929, 1109.6637079956267, 1235.5958061224492, 37.811322448979595,  197.3289291350826, 406.83600174927108,  678.2393532361516, 1030.8756293488823, 1459.7552041982503, 1699.3395257842567, 1873.7209306122447, 35.913513877551026, 261.01980089795921, 566.42278024781342, 972.77621042565625, 1344.9703793955296, 1775.8852377376093, 2011.4623847731782, 2179.8683204081635, 35.495675918367347, 329.87013937414963, 753.12624125364414, 1329.0782355218657, 1624.6724479416907, 1918.8142265189504,  2079.212240188921, 2193.7241183673468, 35.077837959183668,  398.7204778503401,  939.8297022594752, 1685.3802606180761, 1904.3745164878519, 2061.7432153002915, 2146.9620956046647,  2207.579916326531, 34.659999999999997, 467.57081632653058,  1126.533163265306, 2041.6822857142861, 2184.0765850340135,  2204.672204081633, 2214.7119510204084, 2221.4357142857143, 59.250816326530618, 47.128581146744416, 52.486370262390672,  64.53625655976677, 78.221042759961136, 91.508804664723044, 106.61789912536443, 122.33448979591837, 46.594284110787171, 74.103621668749128,  130.7816618231987, 197.20189674302378, 275.04440522490631, 363.49796574760524, 465.58315291128707, 571.30487755102035, 42.618510568513116, 122.65845569901431, 237.68908620106205, 378.80234543940037,  555.9471837758573, 765.99438642232406, 1015.2798367555188, 1274.8845408163265, 38.642737026239075, 171.21328972927947,  344.5965105789254, 560.40279413577684, 836.84996232680828,  1168.490807097043, 1564.9765205997501, 1978.4642040816325, 36.504910553935865, 223.95387163293077, 472.82002618284048,   795.880128433153, 1177.0806018125782, 1541.7735672786343, 1944.2556298882137, 2360.1483955102044, 35.985314518950432,  280.3798178881022, 619.81140958350682, 1078.7936580224907, 1569.5467773016796, 1889.3350020924615,  2173.489664045981, 2458.4179779591836, 35.465718483965013, 336.80576414327362, 766.80279298417315, 1361.7071876118287, 1962.0129527907811, 2236.8964369062887, 2402.7236982037489, 2556.6875604081633, 34.946122448979594,  393.2317103984451,  913.7941763848396, 1644.6207172011666, 2354.4791282798833, 2584.4578717201166, 2631.9577323615158, 2654.9571428571426, 60.257346938775513, 48.615551020408162, 52.918950437317783, 64.014752186588922, 78.448797862001925, 91.900151603498529, 106.96968338192421, 122.64102040816326, 48.465252623906707, 66.897183684575879, 113.17199648583922, 166.66955366097463, 228.40143019644592, 298.05657811745107, 378.86752377176185, 462.70201137026237,  44.09741290087463, 106.25452511453561, 199.54193129164929, 309.81156245314457, 444.49987424857693, 603.41805562265722, 795.35004755518548, 996.17883600583082, 39.729573177842568, 145.61186654449534, 285.91186609745932, 452.95357124531449, 660.59831830070789, 908.77953312786326, 1211.8325713386089, 1529.6556606413994, 37.380458425655981, 188.00277742329587, 385.61219319658477, 631.93261384423158, 944.09728998750529, 1306.1641338558936, 1665.9893122109127, 2021.9097616909623, 36.808740991253643, 233.06461092877964, 497.04933628071626, 842.46456666389008, 1286.9394207302512, 1784.5709910837149, 2153.3165212601416, 2477.8690928279875, 36.237023556851312, 278.12644443426348, 608.48647936484792, 1052.9965194835486, 1629.7815514729973, 2262.9778483115365, 2640.6437303093708, 2933.8284239650147, 35.665306122448975, 323.18827793974737, 719.92362244897959, 1263.5284723032073, 1972.6236822157434,  2741.384705539358, 3127.9709393586008,  3389.787755102041, 61.325612244897961,  50.59113216715258, 54.021231778425658, 63.819169096209912,  77.43923955296404, 91.368839650145773, 106.68598221574345, 122.44173469387756, 50.753077988338191, 60.142305490073589, 96.994053412640568, 138.69980048937944, 186.09616485769817, 239.03959010412328, 299.92605376384847, 362.97087303207002, 46.086638119533532, 90.367618723101486, 164.03767859095169, 248.31595218138278, 347.89711251908926,  463.5815588088297, 602.05845435131198, 746.71301494169097, 41.420198250728866, 120.59293195612938,  231.0813037692628, 357.93210387338604,  509.6980601804803, 688.12352751353592, 904.19085493877537, 1130.4551568513118, 38.791289416909621, 153.32815000000002, 306.20374195960017, 489.37279846563945,  714.3394570284604, 981.06729901374433, 1274.5388618029158, 1573.7532998833822, 37.956335801749269, 188.27322730806605, 388.43921475635148, 640.02902988088294, 956.69995904262112,  1334.235807700125, 1704.9476681212827, 2069.4878461807575, 37.121382186588917, 223.21830461613217, 470.67468755310279, 790.68526129612678, 1199.0604610567816, 1687.4043163865053, 2135.3564744396504, 2565.2223924781338, 36.286428571428566, 258.16338192419823, 552.91016034985432, 941.34149271137051, 1441.4209630709424, 2040.5728250728862, 2565.7652807580175,   3060.95693877551, 62.580714285714294, 53.062978620019436, 55.941967930029158, 64.342813411078723, 76.733697764820221, 90.376813411078729, 105.45567172011664, 120.92336734693879, 53.418368440233237, 54.627930309593232, 82.141676134943765,  112.9480761037068, 147.23920874739693, 185.00273921907541,  227.9410196228655, 272.29399329446068, 48.669044278425659, 75.412394953491599, 130.18380627863391, 191.08248773948361, 260.36504014906984, 338.88505102561436, 431.28950164150365, 527.49252842565602,  43.91972011661808, 96.196859597389974, 178.22593642232403, 269.21689937526031, 373.49087155074272, 492.76736283215331, 634.63798366014157,  782.6910635568513, 41.065402682215748,  119.4345052172706, 230.90353145876719, 358.76728455060402, 509.74262913105656, 684.74303085047904, 894.60638631270331, 1113.9592523615161,  39.87955416909621, 144.83206729418299, 287.66244662848806, 458.36892536609753, 666.35573347966124, 910.25820035318611, 1204.4261042862142, 1512.2033791253646, 38.693705655976679, 170.22962937109537, 344.42136179820909, 557.97056618159115, 822.96883782826603, 1135.7733698558932, 1514.2458222597256, 1910.4475058892126, 37.507857142857148, 195.62719144800781, 401.18027696793007, 657.57220699708478, 979.58194217687071, 1361.2885393586007, 1824.0655402332366, 2308.6916326530618, 64.054897959183677, 56.073442176870756, 58.770911078717191, 65.722443148688043, 76.583597667638486, 89.116425655976684, 103.39863965014577, 118.14285714285714, 56.483193002915456, 50.214164310703879, 69.740226317159511, 89.513152065805912, 111.98036015896153, 136.14387730945438, 163.05492111911704, 190.72493513119531, 51.914769679300292, 61.358056993613772, 98.565858847875859, 137.87742531236984, 181.42378915903089, 228.87201972094954, 282.43067194814654, 337.67138301749264, 47.346346355685128, 72.501949676523665, 127.39149137859221, 186.24169855893376, 250.86721815910028, 321.60016213244467, 401.80642277717607, 484.61783090378992, 44.353391836734701, 86.397183517423287, 159.37886189087877, 239.29811926863809, 328.83018545883652, 428.69632148104949, 543.74940546855464, 662.94414040816321, 42.747567346938773,  102.7148505600444, 194.15000170762175, 296.48576658725528, 414.29422632791875, 548.44287904373152, 705.56182949604317, 868.89901877550983, 41.141742857142859, 119.03251760266554, 228.92114152436477, 353.67341390587251, 499.75826719700109, 668.18943660641366, 867.37425352353159, 1074.8538971428568, 39.535918367346937, 135.35018464528665, 263.69228134110779, 410.86106122448973, 585.22230806608343, 787.93599416909592, 1029.1866775510202, 1280.8087755102038, 66.255918367346951, 61.240672497570458, 63.702871720116619, 69.019994169096208, 78.722499514091339, 89.730740524781339, 102.64307113702625, 116.10510204081635,  60.58546997084548, 47.646256202970989, 61.682279732403174, 78.013611103706793,  96.19651496598641, 116.16448374427323,  138.1818657463557, 160.72900145772593, 57.465986880466474, 50.636563275718458, 71.759139512182443, 95.263734714702224, 121.06093610301266, 149.07635525822576, 180.02758182944609, 211.76333454810498, 54.346503790087468, 53.626870348465928, 81.835999291961684, 112.51385832569765, 145.92535724003892, 181.98822677217828, 221.87329791253654, 262.79766763848403, 51.959135860058311, 58.707403920033343, 94.991782767596874, 134.45399703623497, 177.69062205608782, 224.71000969429414, 277.35145383673489, 331.52987451895063, 50.216362682215745, 65.628288716368189, 110.85842126197421, 160.52348490462313, 215.53177200888524, 276.06898216743036, 344.83236514285733, 415.84427008746371,  48.47358950437318, 72.549173512703064, 126.72505975635156,  186.5929727730113, 253.37292196168272, 327.42795464056655, 412.31327644897988,  500.1586656559769, 46.730816326530615, 79.470058309037938, 142.59169825072891, 212.66246064139952,  291.2140719144802, 378.78692711370275, 479.79418775510237, 584.47306122449015, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 68.592857142857142, 68.908285714285725, 71.698571428571441, 76.548857142857145, 83.530857142857144, 92.372163265306114, 103.68302857142858, 115.65714285714284, 64.712857142857132, 53.557442176870751, 58.816683673469385, 70.900775510204085, 84.260476190476197, 97.388285714285715, 112.58954693877551, 128.45999999999998, 52.748182653061221, 85.355145024295439, 151.37753268950436, 230.89720982507293, 325.65589296890187, 434.95087268221573, 520.74165700874642, 598.10073367346945, 48.979854591836727, 142.05208967444122, 279.45563310860052, 454.76730240524785, 678.42769340379004, 947.53809628279896,  1110.234687303207, 1236.1060280612246, 45.211526530612247, 198.74903432458694, 407.53373352769665, 678.63739498542282, 1031.1994938386783, 1460.1253198833817, 1699.7277175976676, 1874.1113224489795, 43.128146530612241, 261.82915098542276, 566.74532498542271, 972.94391532361544, 1343.6101268299321, 1775.8477002682216, 2011.8575777714286, 2180.1909906122451, 42.528288163265302, 330.52936596307086, 753.36856699708437,  1329.266569399417, 1620.4845516734692,  1918.239364804665, 2079.7815348408158, 2194.0249461224489, 41.928429795918362, 399.22958094071907, 939.99180900874626, 1685.5892234752189, 1897.3589765170066, 2060.6310293411079, 2147.7054919102038, 2207.8589016326532, 41.328571428571422, 467.92979591836735, 1126.6150510204079, 2041.9118775510206, 2174.2334013605441, 2203.0226938775513, 2215.6294489795919, 2221.6928571428571, 65.676326530612243,  54.80744995140914, 59.220189504373181, 70.460827988338195, 84.561978620019431,  97.95702040816326, 113.25078017492712, 129.20999999999998, 54.344540087463557, 78.140149619602937, 133.23961245835071, 198.88613172844651, 276.32858550395667, 364.54594777592672, 466.48218972136618, 572.07231253644318, 50.255220481049562,  125.4922311779814, 239.33458502967511, 379.89978129945865, 556.74639936311257, 766.69030754893799, 1015.9311964150356, 1275.4968250728864, 46.165900874635568, 172.84431273635988,  345.4295576009996, 560.91343087047062, 837.16421322226836, 1168.8346673219492, 1565.3802031087048, 1978.9213376093294, 43.905033119533527, 224.92933695737889, 473.21781435027071, 796.10436262723886, 1177.2455106330697, 1542.6588818225741, 1944.2080931187011, 2358.5864576676386, 43.254035685131186, 281.18150614133003, 620.10604047063725, 1078.9959754735528, 1569.8578447952239, 1891.5484712636403, 2172.8116982307374, 2453.1958697376094, 42.603038250728858, 337.43367532528111, 766.99426659100368, 1361.8875883198671, 1962.4701789573785, 2240.4380607047065, 2401.4153033427738, 2547.8052818075803,  41.95204081632653, 393.68584450923231, 913.88249271137033, 1644.7792011661811, 2355.0825131195334, 2589.3276501457726, 2630.0189084548106, 2642.4146938775511, 66.684081632653061, 56.303830903790086, 59.940612244897956,  70.17620991253645, 84.324340136054403,  98.12841399416908, 113.52489679300294, 129.47367346938776, 56.158813848396498, 71.301270960710823, 116.06153235630985, 168.69136406497293, 229.95792505692071, 299.30983891295296, 379.95501488921292, 463.65766793002916, 51.775805758017491, 109.40699507496873, 201.49621349177426,  311.1482455643482, 445.50000577710671, 604.19835434194079, 796.07460174198263, 996.89266545189503, 47.392797667638483, 147.51271918922669, 286.93089462723862, 453.60512706372356, 661.04208649729276, 909.08686977092862, 1212.1941885947524, 1530.1276629737608,  44.94255055393586, 189.20649423295851, 386.13282219283633, 632.23018844981266, 944.29088821740947, 1306.5787578109125, 1666.1911586747194, 2021.7377268221576, 44.194013294460632,   234.059388121616, 497.45619633069543, 842.69971941357767, 1287.1524385978064, 1785.6037860624738, 2153.5374756568099, 2476.6988927113698, 43.445476034985418, 278.91228201027349, 608.77957046855465, 1053.1692503773427, 1630.0139889782033, 2264.6288143140355, 2640.8837926389006, 2931.6600586005829, 42.696938775510198, 323.76517589893103, 720.10294460641387, 1263.6387813411079, 1972.8755393586005, 2743.6538425655976, 3128.2301096209912, 3386.6212244897961, 67.839591836734684, 58.264557823129252,  61.36025510204081, 70.366043731778433, 83.474786200194359, 97.176367346938775, 112.49735655976679, 128.31193877551021, 58.337543950437308, 65.012728424267664, 100.46169121199499, 141.25476425239486, 188.09759507219215, 240.66786322573927, 301.33180256580602, 364.19080204081638, 53.788442237609331, 93.956994613355533, 166.42428774468971,  250.0277683152853, 349.21375492329582, 464.62589003019571, 602.95960326790919,  747.5045459183674,  49.23934052478134, 122.90126080244343, 232.38688427738438, 358.80077237817579, 510.32991477439958,  688.5839168346522, 904.58740397001247, 1130.8182897959182, 46.572032128279886, 154.89752658808834, 306.92795213036231, 489.79564979508552, 714.63175545411639, 981.29252155351946, 1274.7599770912957, 1573.9865580174928, 45.561558833819234, 189.58094192364291, 389.02198144939598, 640.35589578175757, 956.95666236623606, 1334.5329443981673,  1705.283165777093, 2069.8540658892125, 44.551085539358596, 224.26435725919757,  471.1160107684297, 790.91614176842984, 1199.2815692783561, 1687.7733672428153, 2135.8063544628903, 2565.7215737609331,  43.54061224489795, 258.94777259475217, 553.21004008746354, 941.47638775510222, 1441.6064761904761, 2041.0137900874633, 2566.3295431486886, 3061.5890816326532, 69.200306122448978, 60.671739552964048,  63.52545189504373, 71.280912536443154, 83.119922254616142,  96.22374344023325, 111.05739154518952,  126.3741836734694, 60.839812172011662, 59.956959480771907, 86.246452092877973, 116.21294913369431, 149.90872362973764,  187.2400733777593, 229.89253012411501,  273.9956650145773, 56.347050109329444, 79.527236561155092, 133.09589970324865, 193.32011869012916, 262.16258492468415,  340.3645023531862, 432.55602631820079,   528.571863702624, 51.854288046647234, 99.097513641538256, 179.94534731361932, 270.42728824656399,  374.4164462196307, 493.48893132861309, 635.21952251228663, 783.14806239067059, 49.047269154518943, 121.58750256337639, 231.96392633381924, 359.42442482215756, 510.21098856462595, 685.08667374843822, 894.86602422124145, 1114.1441776967933, 47.724471953352769, 146.64816756573651, 288.53369505830898,  458.8901471603499,  666.7320508474246, 910.55845305289461, 1204.6835212981259, 1512.4246218658895,  46.40167475218658, 171.70883256809662, 345.10346378279883, 558.35586949854246, 823.25311313022348,  1136.030232357351, 1514.5010183750107, 1910.7050660349855, 45.078877551020412, 196.76949757045679, 401.67323250728862,  657.8215918367348, 979.77417541302248, 1361.5020116618075, 1824.3185154518956, 2308.9855102040819, 70.790204081632652, 63.582870748299314, 66.486377551020396, 72.930787172011648, 83.412020408163258, 95.479591836734684, 109.40275043731779, 123.83367346938775, 63.684021282798824, 56.032394302374016,  74.93519507496876, 93.650451928363182,  115.5540697174788, 139.25062147022072, 165.81738906955434, 193.17544708454807, 59.510408892128275, 66.128685242260161, 102.34482181903373, 140.79374589754266, 183.88524464632786, 230.97259345064549, 284.27565736359838, 339.28806814868796, 55.336796501457719, 76.224976182146321, 129.75444856309866, 187.93703986672216, 252.21641957517693, 322.69456543107026, 402.73392565764254, 485.40068921282784, 52.494588921282798, 89.396448133833118, 160.94080566222405, 240.32216424822997,  329.5788563884492, 429.25225181341102, 544.19912061474383, 663.31013183673463, 50.824623906705533, 105.27547983784532, 195.45240542690539, 297.32248656393165, 414.89291656865174, 548.87212317201147, 705.92088857309432, 869.21525795918342, 49.154658892128282, 121.15451154185753, 229.96400519158675, 354.32280887963344, 500.20697674885452, 668.49199453061192, 867.64265653144503, 1075.1203840816322,  47.48469387755101, 137.03354324586974, 264.47560495626817, 411.32313119533524, 585.52103692905735, 788.11186588921259, 1029.3644244897957, 1281.0255102040815, 72.905510204081637, 68.563803692905736, 71.077120991253636, 75.914781341107869, 85.589916423712339, 96.362594752186595, 109.06196034985425, 122.33775510204082, 67.393582507288613, 54.540620347077606, 67.736975036443141, 83.472147576010002, 100.95758496806886, 120.27914075385254,  141.8996332428155, 164.11114679300292, 64.640740524781336, 56.773340733027908, 76.831404427842585, 99.628891930445661, 124.77043066257116, 152.23204951062058, 182.83833605997503, 214.28023287172016, 61.887898542274044, 59.006061118978209, 85.925833819241987, 115.78563628488132, 148.58327635707346, 184.18495826738859,  223.7770388771346, 264.44931895043737,  59.72750903790088, 63.514794050812185, 98.378039814660596, 136.99715897542697, 179.67971314396786, 226.31265821907547, 278.70245821824255, 332.66144845481074, 58.088747521865884, 70.027454537276157, 113.78661837984177, 162.65918858808837, 217.18902709759828, 277.39880810329038, 345.94255776259911, 416.75967311953377, 56.449986005830894, 76.540115023740128, 129.19519694502293, 188.32121820074977, 254.69834105122874, 328.48495798750537, 413.18265730695566, 500.85789778425681, 54.811224489795919,   83.0527755102041, 144.60377551020412,  213.9832478134112, 292.20765500485919,  379.5711078717203, 480.42275685131227, 584.95612244897984, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 75.034285714285716,  75.34114285714287, 78.131428571428572, 82.981714285714276, 89.963714285714289, 98.805020408163259, 110.11403673469388, 122.08571428571427, 71.116428571428571,  61.10750340136056, 65.428954081632668, 76.696102040816342, 90.798411564625852, 104.04695918367348, 119.28357346938778, 135.22428571428571, 60.362241326530615, 89.317583751214769,  153.7344131377551, 232.49880004373185, 326.88496795432462, 435.96384209912532, 521.61425893002922, 598.84723979591854,  56.47130229591837, 144.82169280855203, 281.04043463010203, 455.81168626093302, 679.18658752429531,  948.2211609329446, 1110.8727819314872, 1236.6947831632651, 52.580363265306126, 200.32580186588919,  408.3464561224489, 679.12457247813416, 1031.4882070942663, 1460.4784797667635, 1700.1313049329447, 1874.5423265306119, 50.359772244897968, 262.76034708260454, 567.13951286443148, 973.16545342857171, 1342.2924016394561, 1775.8486139650147, 2008.8059679900875, 2172.4742783673464, 49.609848163265312, 331.29683229543241,  753.6555000728863, 1329.4791798367348, 1616.5599639047616,  1917.868361690962, 2070.4390541702624, 2171.1257093877548, 48.859924081632663, 399.83331750826045,   940.171487281341, 1685.7929062448982, 1890.8275261700678, 2059.8881094169096, 2132.0721403504372, 2169.7771404081632, 48.109999999999999, 468.36980272108843, 1126.6874744897959, 2042.1066326530618, 2165.0950884353742, 2201.9078571428572, 2193.7052265306124, 2168.4285714285716, 72.071530612244899, 62.379587949465503, 66.047248542274062, 76.487087463556861, 90.734504373177842, 104.29795626822157, 119.86950524781344, 136.14938775510203, 61.923336734693891, 82.437039496737469, 135.95604211786758, 200.77318869221165, 277.78148773948357, 365.73044251561856, 467.51764844752199, 572.98639672011666, 57.782372448979594, 128.53930984138555, 241.17165438879636, 381.14853066430663, 557.67450927998755, 767.44354806330705, 1016.6428506377554,  1276.199816144315, 53.641408163265311, 174.64158018603362, 346.38726665972513, 561.52387263640162, 837.56753082049136, 1169.1566536109956, 1565.7680528279886, 1979.4132355685131,  51.29387250728864, 226.04059331139806, 473.70429810495631, 796.38402884048344, 1177.4181354766072, 1543.6723155027073, 1943.9937318179932, 2356.4179608746358, 50.525370787172001, 282.10315226266835, 620.48064358600584, 1079.2156343873387, 1570.0587184356518, 1894.2418410445646, 2171.7500013869221, 2446.2104773177848, 49.756869067055405, 338.16571121393861, 767.25698906705543, 1362.0472399341945, 1962.6993013946962, 2244.8113665864221, 2399.5062709558515, 2536.0029937609329,  48.98836734693878, 394.22827016520898, 914.03333454810502, 1644.8788454810501, 2355.3398843537416, 2595.3808921282798, 2627.2625405247813, 2625.7955102040823, 73.108877551020413, 63.866870748299327, 67.050681486880478, 76.465256559766772, 90.012107871720104, 103.97740233236151, 119.72960845481053, 136.04163265306124, 63.676325218658896, 75.948183138275738, 119.22094099854226, 170.94788845481051, 231.71317860648344, 300.72912876718038, 381.19249045418582, 464.75242514577258, 59.334055029154513, 112.77115832465637, 203.65228106778423, 312.65745029154522, 446.64875697365676, 605.09536119325287, 796.89424749947955, 997.68394697521887, 54.991784839650144, 149.59413351103709, 288.08362113702623, 454.36701212827995, 661.58433534083008, 909.46159361932519, 1212.5960045447732, 1530.6154688046649, 52.477782128279884, 190.56649828987923, 386.75012481257801, 632.58693203248663, 944.54386187505202, 1307.0270225297795, 1666.3319479824247, 2021.3904437900878,  51.57348740524781, 235.19221557517699, 497.95005526863793, 842.95259414410657, 1287.3954414443979, 1786.6501880033318, 2153.5552965169513, 2475.0484591253639, 50.669192682215737,  279.8179328604748,  609.1499857246979, 1053.3182562557267,  1630.247021013744,  2266.273353476884, 2640.7786450514786, 2928.7064744606409, 49.764897959183671, 324.44365014577261, 720.34991618075787, 1263.6839183673469, 1973.0986005830903, 2745.8965189504379, 3128.0019935860059, 3382.3644897959189, 74.370612244897956, 65.803327502429539, 68.732835276967933, 77.011565597667627, 89.534034499514092,  102.9295918367347, 118.12460670553939, 133.85030612244901, 65.754637500000001, 70.115757397265043, 104.18934396084964, 144.06404190753855, 190.32940926627796, 242.50130713869225, 302.92477395428995, 365.58242299562687, 61.359223852040813, 97.761924511488274, 169.01513843190335, 251.93054512703043, 350.70094811710396, 465.82326445231155, 603.99654023037294, 748.41451703717212, 56.963810204081639, 125.40809162571152, 233.84093290295706, 359.79704834652227,  511.0724869679301, 689.14522176593096, 905.06830650645588, 1131.2466110787173, 54.289249373177846, 156.64692293141749, 307.76790488442316, 490.30485743773448, 714.98446617576019, 981.55156518117462, 1275.0211279678056, 1574.2766264285719, 53.129822711370259, 191.04893485353324, 389.70805782382331, 640.74732672719699, 957.23186156129384, 1334.7836299458561, 1705.6206469736778,  2070.308227142857, 51.970396049562687, 225.45094677564904, 471.64821076322357, 791.18979601665978, 1199.4792569468277, 1688.0156947105372, 2136.2201659795505, 2566.3398278571431, 50.810969387755108, 259.85295869776479, 553.58836370262384, 941.63226530612258, 1441.7266523323615, 2041.2477594752188, 2566.8196849854235, 3062.3714285714291, 75.839234693877557, 68.159758503401363, 71.049289358600575,  78.20687900874637, 89.589525267249755, 102.19322886297378, 116.65451122448984, 131.64469387755105, 68.109530357142859, 65.464108386436223, 90.689632043940023, 119.72577230112455, 152.82222842843262, 189.71187807788422, 232.06940307965439, 275.91410557580184, 63.883496811224489, 83.839223197799541,  136.2708313072678, 195.75434195647651, 264.14824619602945, 342.02074521553527, 433.99063702337577, 529.81164622813424, 59.657463265306134, 102.21433800916287, 181.85203057059559, 271.78291161182841,  375.4742639636263, 494.32961235318623, 635.91187096709723, 783.70918688046652, 56.917552434402346, 123.95421768013331, 233.16818112609334, 360.19591574885476, 510.76567583458291, 685.50844932944619,  895.1930957108292, 1114.3835426093297, 55.486106384839658, 148.65662276829102, 289.53370056122446, 459.51286219991675, 667.15754854796603,  910.9105614577262, 1204.9772345555189, 1512.6564603790089, 54.054660335276971, 173.35902785644871, 345.89921999635567, 558.82980865097886, 823.54942126134949, 1136.3126735860058, 1514.7613734002086, 1910.9293781486883, 52.623214285714297, 198.06143294460642, 402.26473943148687, 658.14675510204097, 979.94129397473284, 1361.7147857142859, 1824.5455122448986, 2309.2022959183678, 77.530306122448977, 71.003030126336242, 74.023695335276969, 79.926501457725948, 90.216408163265299, 101.89638775510203, 115.49157172011661, 129.63571428571427, 70.753756924198242, 61.966804796612529, 80.353395212932128, 98.013037594752177, 119.36155244828542, 142.60704160974589, 168.83705797147022, 195.88813163265306, 66.953900327988336, 71.077067246286262,  106.3342094765202, 143.90548917638483, 186.54329739171175, 233.26956297896703, 286.31559120314444, 341.09845153061218, 63.154043731778422, 80.187329695960003, 132.31502374010827, 189.79794075801746, 253.72504233513806, 323.93208434818814,  403.7941244348188, 486.30877142857128, 60.468104577259474, 92.651954608357613, 162.68475192211577, 241.48909065805918, 330.46072488032763, 429.92019317784252, 544.73238800832985, 663.72998431486872, 58.762920058309035,  108.0699460913508,  196.9187230987088, 298.28573585172842, 415.60822345439385, 549.40179836734671, 706.32589803082033, 869.51145212827976, 57.057735539358603,   123.487937574344, 231.15269427530185, 355.08238104539771,  500.7557220284603, 668.88340355685102, 867.91940805331103, 1075.2929199416906, 55.352551020408157, 138.90592905733723,   265.386665451895, 411.87902623906706, 585.90322060252663, 788.36500874635533, 1029.5129180758017, 1281.0743877551017, 79.534897959183681, 75.776611273080661, 78.324952623906697, 82.949839650145776, 92.448035957240037, 102.97697084548105, 115.48754839650147, 128.59897959183672, 74.118098760932938, 61.482118946966544,  73.75662990680965, 88.960537865472745, 106.03761888588087, 124.63544445231153, 145.84836037463558, 167.76216698250732, 71.690979409620979, 63.023278382271279, 81.975647158735953,  104.1109291493128, 128.76367548070249, 155.61517139733445, 185.86418990160354, 217.03057707725949, 69.263860058309035, 64.564437817576021, 90.194664410662227, 119.26132043315289,  151.4897320755241,  186.5948983423574, 225.88001942857156, 266.29898717201178, 67.311766676384849, 68.536716697764859, 101.99149626822161, 139.77439929612666, 181.88766535860066, 228.11040336109963, 280.23454139358614, 333.96380341107897, 65.777912478134112, 74.649487840621987, 116.93843434402336, 165.00908543940034, 219.04034396112735, 278.90219116201592, 347.21364265889224, 417.82583492711393, 64.244058279883376, 80.762258983479128,  131.8853724198251, 190.24377158267401,  256.1930225636541, 329.69397896293219, 414.19274392419851, 501.68786644314889, 62.710204081632661, 86.875030126336284, 146.83231049562687, 215.47845772594764, 293.34570116618085, 380.48576676384857, 481.17184518950467, 585.54989795918402,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142,  81.47571428571429, 81.776142857142872, 84.566428571428574, 89.411000000000001, 96.393000000000015, 105.23816326530611, 116.54554897959186, 128.51357142857142, 77.474285714285713,  68.44531972789116, 72.204897959183668, 82.788734693877558, 96.723258503401368, 110.79787755102041, 126.38690612244899, 142.46285714285713, 67.667463265306125, 93.807839348882425, 156.59217142857142, 234.52539655976679, 328.47526393100094, 437.27670935860056, 522.77833685131202, 599.89384795918363, 63.765147959183679, 148.03160818756075, 283.00231122448974, 457.18151413994167, 680.22830866132153, 949.03724249271113, 1111.6265208454811, 1237.4319566326528, 59.862832653061233, 202.25537702623905, 409.41245102040807, 679.83763172011663, 1031.9813533916422, 1460.7977756268219, 1700.4747048396498, 1874.9700653061223, 57.584645306122454, 263.98708951992234, 567.72157048104964, 973.50615118367386, 1344.8584785170067, 1776.6259605481048,  2013.458392709038, 2183.1201065306122, 56.736430204081628, 332.32921160738579, 754.11631909620985, 1329.6979239183675, 1623.5070899863942, 1919.9444498892126, 2083.5556930985422, 2201.2586424489796, 55.888215102040817,  400.6713336948493, 940.51106771137029,  1685.889696653061,  1902.155701455782, 2063.2629392303206, 2153.6529934880464, 2219.3971783673469, 55.040000000000006, 469.01345578231292, 1126.9058163265306, 2042.0814693877553, 2180.8043129251701, 2206.5814285714282,  2223.750293877551, 2237.5357142857142, 78.501020408163271, 69.729063168124412, 73.070954810495635, 82.784017492711385, 96.394104956268237,  110.2111137026239, 125.94339533527699,  142.2895918367347, 69.189890816326539, 87.218506250173562,  139.2143690285298, 203.11597741774267, 279.63502181174516, 367.26084326114125, 468.86071678134124, 574.17517536443165, 65.098982142857153, 132.01717924475915, 243.42657058256978, 382.74393987921701, 558.91083474593916, 768.43128653685972, 1017.4969926530616, 1276.9501206268224, 61.008073469387767, 176.81585223934471, 347.63877213660976, 562.37190234069146, 838.18664768013321,  1169.601729812578, 1566.1332685247814, 1979.7250658892131, 58.633778950437325, 227.48386624517568, 474.41526255102048, 796.81778539608513, 1177.7114368127172, 1543.8841030154103, 1944.2479074242406, 2356.7501642565603, 57.770886647230327, 283.31957328557547, 621.05860311953347, 1079.5284030920452,  1570.282762863807, 1894.4927295127027, 2172.2258488756356, 2446.9670482798833, 56.907994344023336,  339.1552803259753, 767.70194368804653, 1362.2390207880053, 1962.8540889148965, 2245.1013560099959, 2400.2037903270302, 2537.1839323032073, 56.045102040816332, 394.99098736637518, 914.34528425655981, 1644.9496384839654, 2355.4254149659864, 2595.7099825072887, 2628.1817317784257, 2627.4008163265307, 79.643061224489799, 71.244176870748305, 74.275182215743442, 83.025626822157435, 95.978670553935871, 109.51602332361514, 124.94520291545192, 140.95489795918365, 70.884165451895043, 81.046130712203265, 122.94533385568512, 173.70611451895047, 233.91301293072331, 302.55978641815909, 382.79237429404424, 466.15409839650147, 66.664108965014577, 116.56238200749688, 206.25204783163264, 314.54171836734696, 448.13774433569347, 606.30778095585174, 797.92182852977953, 998.54898505830909, 62.444052478134118, 152.07863330279048, 289.55876180758014, 455.37732221574345, 662.36247574066351, 910.05577549354427, 1213.0512827655145, 1530.9438717201165, 59.922142099125381, 192.29664734194094,  387.6190901811745, 633.12402713536039, 944.92529642399006, 1307.3114491811746, 1666.5360986708877, 2021.5162611661813, 58.895373644314866, 236.65435381591004, 498.66931813827563, 843.36931060724692, 1287.6567573827572, 1786.8964646630568, 2153.7911001829239, 2475.2658067638476, 57.868605189504379, 281.01206028987923, 609.71954609537681, 1053.6145940791337, 1630.3882183415242, 2266.4814801449393, 2641.0461016949603, 2929.0153523615159, 56.841836734693878, 325.36976676384842, 720.76977405247794, 1263.8598775510206, 1973.1196793002912, 2746.0664956268224, 3128.3011032069971, 3382.7648979591836, 80.997551020408167, 73.136744412050547, 76.103913994169105, 83.814043731778426, 95.859140913508256, 108.77187755102042, 123.60906501457728, 139.01520408163265, 72.869984693877541, 75.648548958073036, 108.43719645980842, 147.41752441066225, 193.06088923712343, 244.79826807788422, 304.94924074552279, 367.37385036443152, 68.675931122448986, 102.00039773878939, 172.03843934818821, 254.25492510412332, 352.56774447278912, 467.37497254789673, 605.34028366305722, 749.57521118804664, 64.481877551020418, 128.35224651950574, 235.63968223656804, 361.09232579758435,  512.0745997084548, 689.95167701790922, 905.73132658059149, 1131.7765720116618, 61.848565976676397, 158.81312604303767,  308.8940714827155, 491.04388042399012, 715.51386082465638, 981.95552094877144, 1275.3440890780512, 1574.5272406122451, 60.589418134110787, 192.89182390559489, 390.64762564348183, 641.34641688129955, 957.63364628349302, 1335.0868652778011, 1705.9034889286133, 2070.5888746938772, 59.330270291545197,  226.9705217681522, 472.40117980424816, 791.64895333860898, 1199.7534317423294, 1688.2182096068304, 2136.4628887791755, 2566.6505087755099, 58.071122448979594, 261.04921963070944, 554.15473396501443, 941.95148979591852,  1441.873217201166, 2041.3495539358601, 2567.0222886297379, 3062.7121428571431, 82.537346938775514, 75.486156462585043, 78.299832361516053, 84.776157434402336, 96.129660835762877, 108.37130903790089, 122.50582448979594,  137.2090816326531, 75.102321428571429, 71.352195876023885, 96.477119200333192, 123.75539423573514, 156.25067502776625, 192.69318818408993, 234.74493505310295, 278.31844147230328, 71.145905612244903,  88.57989334825767, 140.39090472719698, 198.61867109537695, 266.54059115125642, 344.07840390462314, 435.81158364483559, 531.41785513848401, 67.189489795918377, 105.80759082049147,  184.3046902540608, 273.48194795501877, 376.83050727474665, 495.46361962515618, 636.87823223656824, 784.51726880466481,  64.54698638483967, 126.79769490045817, 234.71217373177848, 361.25402155268654,  511.5562338123006, 686.13363199416926,  895.7051578639738, 1114.7924177259479, 63.061324256559772, 151.10043022518397,  290.8370686880466, 460.39171504872979, 667.79660177703738, 911.39215602332376, 1205.3874317732616, 1513.0173601166182, 61.575662128279887, 175.40316554990977, 346.96196364431489,  559.5294085447731, 824.03696974177444, 1136.6506800524783, 1515.0697056825495,  1911.242302507289, 60.090000000000003,  199.7059008746356, 403.08685860058307, 658.66710204081653, 980.27733770651139, 1361.9092040816327, 1824.7519795918374, 2309.4672448979595, 84.257346938775498, 78.357035957240043, 80.931392128279882, 86.017586005830893,  96.86079591836733, 108.37914285714285, 121.71222740524782, 135.62714285714281, 77.590118221574357, 68.157757470498396, 87.654743164306538, 102.87241664723032, 123.64687196098848, 146.49584510620573, 172.40522453102872, 199.15333163265302, 74.111722667638475, 76.404518818547814, 111.62326841680547, 147.46283724489794, 189.61647861828402, 236.01259443981667, 288.79394121511865, 343.32613520408154,  70.63332711370262, 84.651280166597246, 135.59179366930442, 192.05325784256559, 255.58608527557951, 325.52934377342763, 405.18265789920861, 487.49893877551006, 68.106335568513117, 96.447911319727879, 164.85059914827153, 243.00784962432323, 331.64730254671656, 430.88739835568504,  545.5275787958351, 664.36136029154511,  66.41701282798833,  111.3700445403304, 198.76726046230729, 299.56581034235734, 416.59373295793404, 550.19301581341085, 706.96486924581404, 870.00553271136994, 64.727690087463557, 126.29217776093292, 232.68392177634311,  356.1237710603915, 501.54016336915163, 669.49863327113678, 868.40215969579333, 1075.6497051311951, 63.038367346938777, 141.21431098153545, 266.60058309037896, 412.68173177842567, 586.48659378036905, 788.80425072886271, 1029.8394501457724, 1281.2938775510202, 86.110408163265305, 82.820602526724983, 85.276107871720114, 89.844233236151609, 99.223706511175891,  109.5586472303207,  121.9212390670554, 134.90265306122447, 80.695893731778426,  68.42347623351381, 79.765203623490208, 94.423290637234501, 112.02711257600998, 129.57276566014161, 150.41864057434404, 172.31799125364435, 78.521977769679296, 69.470886738164666, 87.300518325697638, 108.77703617242818, 133.49534004754963, 159.53380981882552, 189.43835691690967, 220.48550291545195,  76.34806180758018, 70.518297242815493, 94.835833027905053, 123.13078170762185, 154.96356751908928, 189.49485397750942, 228.45807325947533, 268.65301457725957, 74.568589212828002, 74.042213613216745, 106.07007205539361, 143.03097594668895,  184.5614880583091, 230.35188345356113, 282.18602880583109, 335.66091451895062, 73.136406413994166,  79.74658296598642, 120.56104900874639, 167.81457094210748, 221.31724082410113, 280.80234212744699, 348.86393310087476, 419.25693620991268, 71.704223615160359, 85.450952318756094, 135.05202596209915, 192.59816593752612,  258.0729935898932, 331.25280080133291, 415.54183739591861, 502.85295790087486, 70.272040816326538, 91.155321671525769, 149.54300291545195,  217.3817609329447, 294.82874635568521, 381.70325947521883, 482.21974169096234, 586.44897959183697, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 87.917142857142863, 88.211142857142875, 91.001428571428576, 95.846000000000004,            102.828, 111.66820408163265, 122.97134693877553, 134.93571428571428, 83.903571428571439, 75.695653061224505, 78.991045918367348, 89.022081632653055, 102.55390476190476, 116.74389795918367, 132.48099387755104, 148.68000000000001, 74.841157142857142, 98.512903620019458, 159.71877328717201, 236.79151107871724, 330.28267158163266, 438.78504549562683, 524.11675836588927, 601.09353775510203, 70.969178571428571, 151.44579865160352,  285.1811609876093, 458.73780685131192, 681.44096823979578, 950.01011359329436, 1112.4814087572886, 1238.2008647959183, 67.097200000000015, 204.37869368318755, 410.64354868804662, 680.68410262390671, 1032.5992648979588, 1461.2351816909618, 1700.8460591486878, 1875.3081918367347, 64.790626530612258, 265.39258289990289, 568.43716552478145, 973.95164071720137, 1346.1999458950436, 1777.0899578309036,  2013.357338876968, 2182.6933979591836,  63.86232244897959, 333.52142714868802, 754.69698790087455, 1330.0143523148688, 1626.7328210728861, 1920.9299038600579,  2082.992026598251, 2199.7729795918367, 62.934018367346937, 401.65027139747326, 940.95681027696787, 1686.0770639125367, 1907.2656962507285, 2064.7698498892123, 2152.6267143195332, 2216.8525612244898, 62.005714285714291, 469.77911564625845, 1127.2166326530612, 2042.1397755102043, 2187.7985714285714, 2208.6097959183671, 2222.2614020408164, 2233.9321428571429, 85.011632653061241,  76.99275121477163, 79.909475218658898, 89.217580174927122, 102.23030660835764, 115.91581632653062, 131.53056967930033, 147.72928571428574, 76.324125291545201, 92.202041249479407, 142.73929148531863, 205.70945496043322, 281.72045146917958, 369.00782678467306,  470.4042903869223, 575.54805131195349, 72.318484511661808, 135.69924768846317, 245.89904954315909, 384.53226397334447,  560.3292306339373,  769.5898447209496, 1018.4962903613082, 1277.8157230320703, 68.312843731778443, 179.19645412744688, 349.05880760099961, 563.35507298625578, 838.93800979869491,  1170.171862657226, 1566.5882903356937, 1980.0833947521867, 65.930658717201183, 229.12403057420522, 475.26509916493137, 797.35643433402765, 1178.0933266980428, 1544.1626983082051, 1944.5030067489386, 2357.0041594752192, 64.977854110787177, 284.71326106955433, 621.76283578925438, 1079.9400504889632, 1570.5571701796475, 1894.7411837117866, 2172.5842135371931, 2447.4714804664723,   64.0250495043732, 340.30249156490351, 768.26057241357762, 1362.5236666438987, 1963.0210136612523, 2245.3196691153689, 2400.6654203254479, 2537.9388014577262, 63.072244897959195, 395.89172206025273, 914.75830903790097, 1645.1072827988341, 2355.4848571428574, 2595.8981545189504, 2628.7466271137027, 2628.4061224489797, 86.246122448979605, 78.522187560738587, 81.211217201166178, 89.568897959183687, 102.08576044703597, 115.19540233236151, 130.29715160349855, 146.00132653061223, 77.963850072886302, 86.334861428571429, 126.91843549302374, 176.72941174094132, 236.35794609919481, 304.62371403998338, 384.61568060058312, 467.76930772594756, 73.890023870262382, 120.55919088921283, 209.06446625754893, 316.64277679092049, 449.82183168211168, 607.70406514993761, 799.12613048104981, 999.58404045189513, 69.816197667638491,  154.7835203498542, 291.21049702207409, 456.55614184089967, 663.28571726502844, 910.78441625989171, 1213.6365803615161, 1531.3987731778425, 67.308610087463563, 194.24309129612669, 388.64022810287383, 633.79940307538538, 945.42147648035552, 1307.6769563823407, 1666.8242789461897, 2021.7412869970849, 66.180025772594746, 238.31205892044983, 499.52657578092453, 843.90997328446485, 1288.0197582107455, 1787.1669456543104, 2154.0655703178677, 2475.5693613994163, 65.051441457725957, 282.38102654477302, 610.41292345897534, 1054.0205434935444, 1630.6180399411355, 2266.6569349262804, 2641.3068616895457,  2929.397435801749, 63.922857142857154, 326.44999416909627, 721.29927113702615, 1264.1311137026241, 1973.2163216715257, 2746.1469241982509, 3128.5481530612242, 3383.2255102040817, 87.676275510204093, 80.411519436345969,  83.15217930029155, 90.144000000000005, 102.16813824101069,  114.7391778425656, 129.24439577259477, 144.35229591836736, 79.866295809037908, 81.381207334791071, 112.91312630284257, 151.02554463765102, 196.04620681729835,  247.3426873885881, 307.21482367794675, 369.40053866618075, 75.881728225218666,  106.4618704020894, 175.27208570973031, 256.79380141087051, 354.64417348153552, 469.11958448563098, 606.87130378774486, 750.92204965379017, 71.897160641399424, 131.54253346938776, 237.63104511661805, 362.56205818408995, 513.24214014577262, 690.89648158267391,  906.5277838975428, 1132.4435606413995, 69.321232944606422, 161.22435991545194, 310.19639106674299, 491.93029715243665, 716.18060848931009, 982.46086048146628, 1275.7552561241571,  1574.862078411079, 67.985549854227401, 194.95730554130222, 391.74800524468969, 642.07726807246991, 958.15891132847423,  1335.470005141191, 1706.2368845978344, 2070.8976951311952, 66.649866763848408, 228.69025116715258,  473.2996194226364, 792.22423899250316, 1200.1372141676384, 1688.4791498009163,  2136.718513071512, 2566.9333118513118, 65.314183673469387, 262.42319679300294, 554.85123360058299, 942.37120991253664, 1442.1155170068027, 2041.4882944606416,   2567.20014154519, 3062.9689285714289, 89.240051020408174, 82.767894557823141,  85.40344569970847, 91.135377551020412, 102.64835447035958, 114.62123323615162,  128.4601166180758, 142.90357142857144, 81.987805648688052, 77.443926325142314, 102.13334478602667, 128.01942767076221,  159.9276182283424, 195.92681224281552, 237.67391453748445, 280.97734679300299, 78.291978771865899, 93.552537081771504, 144.52187487635362, 201.69417678831741, 269.14893723231643, 346.34047752498964, 437.83697131039162, 533.23357981049571, 74.596151895043747, 109.66114783840068, 186.91040496668055, 275.36892590587263,  378.3702562362904,  496.7541428071637, 638.00002808329873, 785.48981282798843, 72.058390233236167, 129.90509701249482, 236.46563447834239, 362.48059144439833,  512.5058066966543, 686.88806597917551, 896.33393332490687, 1115.3125485276971, 70.540253352769696, 153.79002726003054, 292.33082325281134, 461.42286271220337, 668.57722697820361, 911.99390501624327,  1205.884676686964, 1513.4291309912537, 69.022116472303225, 177.67495750756632, 348.19601202728029, 560.36513398000852, 824.64864725975303, 1137.0997440533113, 1515.4354200490218, 1911.5457134548108, 67.503979591836739, 201.55988775510207, 404.06120080174929, 659.30740524781356, 980.72006754130234, 1362.2055830903792, 1824.9861634110794, 2309.6622959183678, 90.924183673469386, 85.644307094266281, 88.026916909620994,  92.71155685131194, 103.60274684159378, 114.87360349854228, 127.97449241982508, 141.68744897959181, 84.335300218658887, 74.525166280716363, 94.433721035506039, 108.16973317367763,   128.160165442871, 150.62835169512701,  176.2276158771345, 202.68193250728859, 81.151640488338174, 81.945652658614463, 116.70858327910243, 151.36832502082467, 192.90384126232121, 238.97190726780502,  291.4963882309454, 345.78760532069964, 77.967980758017504,  89.36613903651255, 138.98344552269884, 194.56691686797166, 257.64751708177141, 327.31546284048301,  406.7651605847563, 488.89327813411069, 75.595643673469397, 100.50975230431763, 167.26746189400248, 244.73344367180346, 333.01790078994856, 432.02143162182415, 546.49203087988337, 665.17021309037887, 73.937640000000002,  114.9314127768985, 200.84227051436898, 301.03492649062889, 417.74462692766883, 551.13377560016636, 707.75316577609306, 870.65296518950402, 72.279636326530621, 129.35307324947937, 234.41707913473545, 357.33640930945438, 502.47135306538928,  670.2461195785088, 869.01430067230308, 1076.1357172886296, 70.621632653061226,  143.7747337220602, 267.99188775510197, 413.63789212827987, 587.19807920310961, 789.35846355685112, 1030.2754355685129, 1281.6184693877547, 92.656326530612247, 89.789088435374154, 92.210189504373176, 96.755944606414005, 105.97121137026238, 116.12039650145772, 128.35422128279882, 141.22418367346938, 87.216933600583104, 75.371044949326674, 85.854285552894638, 100.28586499375263, 118.06296159551576, 135.18667849645979, 155.82607073760937, 177.66545801749274, 85.265884293002927, 76.018229310703873, 92.766461539462739,  113.7766418679717, 138.35370327728032, 163.95455205122869, 193.60977894314874, 224.51014358600588, 83.314834985422749, 76.665413672081087, 99.678637526030855, 127.26741874219078, 158.64444495904485, 192.72242560599756, 231.39348714868817, 271.35482915451905, 81.691243965014593, 79.782375448840781, 110.38329926905459, 146.52308022324041, 187.47577613216723, 232.83102243398594, 284.37239691278648, 337.59061032069985, 80.355965364431484, 85.073866107455245, 124.42707560599754, 170.85446553269475, 223.82671551020414, 282.92462039816752, 350.73000406430674, 420.89938647230338, 79.020686763848403, 90.365356766069709, 138.47085194294047, 195.18585084214919, 260.17765488824114,  333.0182183623491,   417.087611215827,  504.2081626239069, 77.685408163265308, 95.656847424684173, 152.51462827988342, 219.51723615160358, 296.52859426627799, 383.11181632653074, 483.44521836734725, 587.51693877551043, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 94.358571428571423, 94.643625850340143, 97.433239795918382, 102.28100000000001, 109.26005782312924, 118.09742857142857, 129.39868571428573, 141.36142857142858, 90.390000000000001, 82.876000000000005, 85.785357142857137, 95.367999999999995, 108.30919047619048, 122.04599999999999, 137.76768571428573, 154.08000000000001,  81.90962857142857, 103.38981482993199,          163.06045, 239.24924000000001, 332.26376857142861, 440.44975673469389, 525.59465475510206, 602.41569285714286, 78.101357142857125, 155.02340918367349,  287.5335892857143, 460.44327142857145, 682.79037857142862, 951.10841632653057, 1113.4172158673471, 1238.9951607142857, 74.293085714285723, 206.65700353741497, 412.00672857142854, 681.63730285714291, 1033.3169885714285,  1461.767075918367, 1701.2397769795919, 1875.5746285714285, 71.981494285714291, 266.94107703401363, 569.25959051020413, 974.48096367346966, 1346.5617256734695, 1777.3032756734692, 2009.4535022040818, 2173.4086600000001, 70.987662857142851, 334.84151167346937, 755.37353653061211, 1330.4089281632655, 1626.9814361632652, 1921.0428504489794, 2071.4841157551018,         2172.99244, 69.993831428571426, 402.74194631292517, 941.48748255102021, 1686.3368926530616, 1907.4011466530612, 2064.7824252244895, 2133.5147293061223, 2172.5762199999999,                 69, 470.64238095238096, 1127.6014285714284, 2042.2648571428576, 2187.8208571428572, 2208.5219999999999, 2195.5453428571427, 2172.1599999999999, 91.587142857142851, 84.187809523809534, 86.599846938775514, 95.760448979591843, 108.20778911564625, 121.45375510204082, 136.72837142857145, 152.60857142857145, 83.352504081632674, 97.347230845481064, 146.47749037900877, 208.50348355685134, 283.99139759475224, 370.92807653061226, 472.10826820991269, 577.06820510204091, 79.460301020408153,  139.5446753644315, 248.54557871720115, 386.47491997084558, 561.89328285957231,  770.8850586734693,  1019.611712623907, 1278.7735637755104, 75.568097959183675, 181.74211988338192, 350.61366705539359, 564.44635638483965, 839.79516812439249, 1170.8420408163265, 1567.1151570379011, 1980.4789224489796, 73.193117142857147, 230.92170801943638, 476.22603351311955, 797.97899717784276, 1178.5460874227408, 1544.4947398250729,  1944.758845048397, 2357.1955881632657, 72.153982857142864, 286.24876226239064, 622.56808696792996, 1080.4308008396504, 1570.8718678056366, 1894.9876904956268, 2172.8485919370264, 2447.7742016326533, 71.114848571428581, 341.57581650534496, 768.91014042274048,  1362.882604501458, 1963.1976481885322, 2245.4806411661812, 2400.9383388256565, 2538.3528151020409, 70.075714285714298, 396.90287074829934, 915.25219387755101, 1645.3344081632658, 2355.5234285714287, 2595.9735918367351, 2629.0280857142861, 2628.9314285714286,  92.90428571428572,   85.7207619047619, 87.916479591836719, 96.098489795918368, 108.30527210884352, 120.98738775510203, 135.75818367346938, 151.15428571428572, 84.941010204081635, 91.776218658892134, 131.09050415451895, 179.96476588921286,  238.9989583430515, 306.87425763848404, 386.61772488046654, 469.55534591836732, 81.032627551020397, 124.72046793002914, 212.04700601311953, 318.91726749271146,   451.661999016035, 609.24744088921284, 800.47180916909633, 1000.7551096938776, 77.124244897959187, 157.66471720116618,  293.0035078717201, 457.86976909621001, 664.32503968901847, 911.62062413994158, 1214.3258934577261, 1531.9548734693876, 74.648764489795909, 196.36257117201171, 389.78310406705543, 634.58540368513138, 946.00945294266273, 1308.1073280233238, 1667.1796828909623, 2022.0456795918369, 73.434890612244885, 240.12621751603496, 500.49422924198245, 844.54979293294468, 1288.4641069504373, 1787.4567901107871, 2154.3709736279884, 2475.9418816326529, 72.221016734693876, 283.88986386005831, 611.20535441690947, 1054.5141821807581, 1630.9187609582116, 2266.8062521982501, 2641.5622643650149, 2929.8380836734691, 71.007142857142853, 327.65351020408167, 721.91647959183661, 1264.4785714285715, 1973.3734149659863, 2746.1557142857141,  3128.753555102041, 3383.7342857142858, 94.396428571428572, 87.639380952380947, 89.942193877551006, 96.095938775510206, 108.46424829931973, 120.80648979591835,  135.0004244897959, 149.82714285714286, 86.767378061224491, 87.273759164237134, 117.57151802113704, 154.83719504373181, 199.23459448493682, 250.08507339650151, 309.67329952332369,  371.6154357142857,  82.99879719387755, 111.10174261418854, 178.67400842747813, 259.50427478134117, 356.88830861273084, 471.01851949708453, 608.55214526603504, 752.41780357142864, 79.230216326530609, 134.92972606413997, 239.77649883381923, 364.17135451895047, 514.54202274052477, 691.95196559766759, 907.43099100874645, 1133.2201714285713, 76.724580204081633, 163.83161839650148, 311.63963303935867, 492.93462887463573, 716.95723856656957,  983.0473070262392, 1276.2369879189507, 1575.2642951020412, 75.330910612244892, 197.20086124392611, 392.97703427113703, 642.91352809329464, 958.78296063168114, 1335.9170686297375, 1706.6107232384838,  2071.229053877551, 73.937241020408166,  230.5701040913508, 474.31443550291539, 792.89242731195338, 1200.6086826967928, 1688.7868302332358, 2136.9844585580176, 2567.1938126530613, 72.543571428571425, 263.93934693877554,  555.6518367346938, 942.87132653061235, 1442.4344047619047, 2041.6565918367346, 2567.3581938775515, 3063.1585714285716, 95.946428571428584, 90.013904761904769, 92.389515306122448, 97.326551020408175, 109.14989455782313, 120.92863265306121,  134.4967918367347, 148.70214285714286, 88.787444387755102, 83.698571141885338, 107.68456111516038, 132.47099030612247, 163.80335870991254, 199.36228746355687, 240.80565203061232, 283.83990765306123, 85.344983418367349, 98.710759681729854, 148.66156240889217, 204.93862372448984, 271.93008421404278, 348.76608309037908,  440.0259118112246, 535.21691709183676, 81.902522448979596, 113.72294822157437, 189.63856370262391, 277.40625714285716, 380.05680971817299, 498.16987871720119, 639.24617159183686, 786.59392653061229, 79.475369999999998, 133.22363903790091,  238.3866697376094,  363.8419326064141, 513.58259150631693, 687.74590102040838, 897.05607943207042, 1115.9216838775515, 77.942151428571435, 156.67625595724004, 293.97688696793006, 462.57584622740535, 669.47110975704561, 912.69177755102044, 1206.4515597574346, 1513.8815987755104, 76.408932857142872, 180.12887287657921, 349.56710419825077, 561.30975984839665, 825.35962800777452, 1137.6376540816327, 1515.8470400827991, 1911.8415136734698, 74.875714285714281, 203.58148979591837, 405.15732142857144, 660.04367346938784, 981.24814625850343, 1362.5835306122449,  1825.242520408164, 2309.8014285714289, 97.542857142857144, 92.878190476190468, 95.272704081632654, 99.887836734693877, 110.42274829931972, 121.37742857142857,  134.2700448979592, 147.80285714285714, 91.007538775510199, 81.033739999999995,  100.7948028425656, 113.81739967930028, 132.85583810009717, 154.95582075801747,  180.2533870524781, 206.42125408163264, 88.097234693877539, 87.657732312925162, 121.63090287900873,  155.5523245626822, 196.36254904033038, 242.10424518950433, 294.37811284985418, 348.43610459183668, 85.186930612244907,  94.28172462585033, 142.46700291545187,  197.2872494460641, 259.86925998056358, 329.25266962099113, 408.50283864723019, 490.45095510204067, 82.965813469387754,   104.784300707483, 169.88513705539359, 246.62450578425657, 334.53571549465494, 433.28892735860046, 547.59189200116623, 666.12104734693878, 81.351494693877555, 118.70179684353741, 203.09845871720111, 302.65527596501454,  419.0278284703594, 552.19416925947496, 708.66092255860042, 871.42307918367317, 79.737175918367356, 132.61929297959182, 236.31178037900867, 358.68604614577254, 503.51994144606397, 671.09941116034963, 869.72995311603472, 1076.7251110204079, 78.122857142857143, 146.53678911564623, 269.52510204081631, 414.71681632653053, 588.01205442176854, 790.00465306122419, 1030.7989836734691, 1282.0271428571425, 99.178571428571445, 96.697170068027219, 99.130612244897947, 103.68151020408163, 112.69618367346938, 122.66620408163264, 134.78663673469387,             147.56,   93.6925693877551,  82.32358280855199,  92.00777405247814, 106.46829661807581, 124.13589487852283,  141.3418646355685, 161.90322087172012,  183.6462387755102, 91.940117346938763, 82.645359256559772, 98.345262390670555, 119.04304650145774, 143.31342543731779, 168.77697733236155, 198.25900497813419, 228.99055612244905,  90.18766530612244, 82.967135704567553,   104.682750728863, 131.61779638483969, 162.49095599611275, 196.21209002915458, 234.61478908454822, 274.33487346938784, 88.706564489795923, 85.710269220602527, 114.88424762390676, 150.20360660058316, 190.58243650534507, 235.50028852478144,  286.7466695755104, 339.70637387755119, 87.464376326530612, 90.585299661807568, 128.48783174927115, 174.08188739358604, 226.52225245481054, 285.22460051311964, 352.76869944489817,  422.7109159183675, 86.222188163265301, 95.460330103012637,  142.0914158746356, 197.96016818658899, 262.46206840427612, 334.94891250145787, 418.79072931428595, 505.71545795918382,  84.97999999999999, 100.33536054421771, 155.69500000000005, 221.83844897959193, 298.40188435374159,  384.6732244897961, 484.81275918367379, 588.72000000000025,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999,              100.8, 101.07409523809524,           103.8625, 108.71599999999999, 115.68976190476189, 124.52600000000001, 135.82725714285715, 147.78999999999999
};

BOOST_AUTO_TEST_SUITE_END() // Integration tests//


