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
#include <memory>
#include <map>
#include <sstream>
#include <limits>
#include <vector>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/wells.h>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>

#include <opm/autodiff/VFPHelpersLegacy.hpp>
#include <opm/autodiff/VFPProdPropertiesLegacy.hpp>



const double max_d_tol = 1.0e-10;
const double sad_tol = 1.0e-8;


struct ConversionFixture {
    typedef Opm::VFPProdPropertiesLegacy::ADB ADB;

    ConversionFixture() :
            num_wells(5),
            aqua(ADB::null()),
            liquid(ADB::null()),
            vapour(ADB::null())
    {
        ADB::V aqua_v(num_wells);
        ADB::V liquid_v(num_wells);
        ADB::V vapour_v(num_wells);

        for (int i=0; i<num_wells; ++i) {
            aqua_v[i] = 300+num_wells*15;
            liquid_v[i] = 500+num_wells*15;
            vapour_v[i] = 700+num_wells*15;
        }

        aqua = ADB::constant(aqua_v);
        liquid = ADB::constant(liquid_v);
        vapour = ADB::constant(vapour_v);
    }

    ~ConversionFixture() {

    }

    int num_wells;

    ADB aqua;
    ADB liquid;
    ADB vapour;
};





BOOST_FIXTURE_TEST_SUITE( ConversionTests, ConversionFixture )


BOOST_AUTO_TEST_CASE(getFlo)
{
    //Compute reference solutions
    std::vector<double> ref_flo_oil(num_wells);
    std::vector<double> ref_flo_liq(num_wells);
    std::vector<double> ref_flo_gas(num_wells);
    for (int i=0; i<num_wells; ++i) {
        ref_flo_oil[i] = liquid.value()[i];
        ref_flo_liq[i] = aqua.value()[i] + liquid.value()[i];
        ref_flo_gas[i] = vapour.value()[i];
    }

    {
        ADB flo = Opm::detail::getFlo(aqua, liquid, vapour, Opm::VFPProdTable::FLO_OIL);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_flo_oil.begin(), ref_flo_oil.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getFlo(aqua, liquid, vapour, Opm::VFPProdTable::FLO_LIQ);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_flo_liq.begin(), ref_flo_liq.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getFlo(aqua, liquid, vapour, Opm::VFPProdTable::FLO_GAS);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_flo_gas.begin(), ref_flo_gas.end(), computed, computed+num_wells);
    }
}


BOOST_AUTO_TEST_CASE(getWFR)
{
    //Compute reference solutions
    std::vector<double> ref_wfr_wor(num_wells);
    std::vector<double> ref_wfr_wct(num_wells);
    std::vector<double> ref_wfr_wgr(num_wells);
    for (int i=0; i<num_wells; ++i) {
        ref_wfr_wor[i] = aqua.value()[i] / liquid.value()[i];
        ref_wfr_wct[i] = aqua.value()[i] / (aqua.value()[i] + liquid.value()[i]);
        ref_wfr_wgr[i] = aqua.value()[i] / vapour.value()[i];
    }

    {
        ADB flo = Opm::detail::getWFR(aqua, liquid, vapour, Opm::VFPProdTable::WFR_WOR);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_wfr_wor.begin(), ref_wfr_wor.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getWFR(aqua, liquid, vapour, Opm::VFPProdTable::WFR_WCT);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_wfr_wct.begin(), ref_wfr_wct.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getWFR(aqua, liquid, vapour, Opm::VFPProdTable::WFR_WGR);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_wfr_wgr.begin(), ref_wfr_wgr.end(), computed, computed+num_wells);
    }
}


BOOST_AUTO_TEST_CASE(getGFR)
{
    //Compute reference solutions
    std::vector<double> ref_gfr_gor(num_wells);
    std::vector<double> ref_gfr_glr(num_wells);
    std::vector<double> ref_gfr_ogr(num_wells);
    for (int i=0; i<num_wells; ++i) {
        ref_gfr_gor[i] = vapour.value()[i] / liquid.value()[i];
        ref_gfr_glr[i] = vapour.value()[i] / (liquid.value()[i] + aqua.value()[i]);
        ref_gfr_ogr[i] = liquid.value()[i] / vapour.value()[i];
    }

    {
        ADB flo = Opm::detail::getGFR(aqua, liquid, vapour, Opm::VFPProdTable::GFR_GOR);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_gfr_gor.begin(), ref_gfr_gor.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getGFR(aqua, liquid, vapour, Opm::VFPProdTable::GFR_GLR);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_gfr_glr.begin(), ref_gfr_glr.end(), computed, computed+num_wells);
    }

    {
        ADB flo = Opm::detail::getGFR(aqua, liquid, vapour, Opm::VFPProdTable::GFR_OGR);
        const double* computed = &flo.value()[0];
        BOOST_CHECK_EQUAL_COLLECTIONS(ref_gfr_ogr.begin(), ref_gfr_ogr.end(), computed, computed+num_wells);
    }
}



BOOST_AUTO_TEST_SUITE_END() // unit tests
























/**
 * Test fixture to set up axis etc.
 * All of our axes go from 0 to 1, but with a varying number of
 * values data is given at
 */
struct TrivialFixture {
    typedef Opm::VFPProdPropertiesLegacy::ADB ADB;
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
            size{{ nx, ny, nz, nu, nv }},
            data(size)
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
                            data[i][j][k][l][m] = value;
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
                            data[i][j][k][l][m] = x + 2*y + 3*z + 4*u + 5*v;
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
                            data[i][j][k][l][m] = randx / max_val;
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
                                          Opm::VFPProdTable::FLO_OIL,
                                          Opm::VFPProdTable::WFR_WOR,
                                          Opm::VFPProdTable::GFR_GOR,
                                          Opm::VFPProdTable::ALQ_UNDEF,
                                          flo_axis,
                                          thp_axis,
                                          wfr_axis,
                                          gfr_axis,
                                          alq_axis,
                                          data));

        properties.reset(new Opm::VFPProdPropertiesLegacy(table.get()));
    }



    /**
     * Helper function to simplify creating minimal ADB objects
     */
    inline ADB createConstantScalarADB(double value) {
        ADB::V v = ADB::V::Constant(1, value);
        ADB adb = ADB::constant(std::move(v));
        return adb;
    }

    std::shared_ptr<Opm::VFPProdPropertiesLegacy> properties;
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
    Opm::VFPProdTable::extents size;
    Opm::VFPProdTable::array_type data;
};





//Set F to be our test suite fixture for our "trivial" tests
BOOST_FIXTURE_TEST_SUITE( TrivialTests, TrivialFixture )


BOOST_AUTO_TEST_CASE(GetTable)
{
    fillDataRandom();
    initProperties();

    //Create wells
    const int nphases = 3;
    const int nwells  = 1;
    const int nperfs  = 1;
    std::shared_ptr<Wells> wells(create_wells(nphases, nwells, nperfs),
                                destroy_wells);
    const int cells[] = {5};
    add_well(INJECTOR, 100, 1, NULL, cells, NULL, 0, NULL, true, wells.get());

    //Create interpolation points
    double aqua_d   = -0.15;
    double liquid_d = -0.25;
    double vapour_d = -0.35;
    double thp_d    =  0.45;
    double alq_d    =  0.55;

    ADB aqua_adb   = createConstantScalarADB(aqua_d);
    ADB liquid_adb = createConstantScalarADB(liquid_d);
    ADB vapour_adb = createConstantScalarADB(vapour_d);
    ADB thp_adb    = createConstantScalarADB(thp_d);
    ADB alq_adb    = createConstantScalarADB(alq_d);

    ADB::V qs_adb_v(3);
    qs_adb_v << aqua_adb.value(), liquid_adb.value(), vapour_adb.value();
    ADB qs_adb = ADB::constant(qs_adb_v);

    //Check that our reference has not changed
    Opm::detail::VFPEvaluation ref = Opm::detail::bhp(table.get(), aqua_d, liquid_d, vapour_d, thp_d, alq_d);
    BOOST_CHECK_CLOSE(ref.value, 1.0923565702101556,  max_d_tol);
    BOOST_CHECK_CLOSE(ref.dthp,  0.13174065498177251, max_d_tol);
    BOOST_CHECK_CLOSE(ref.dwfr, -1.2298177745501071,  max_d_tol);
    BOOST_CHECK_CLOSE(ref.dgfr,  0.82988935779290274, max_d_tol);
    BOOST_CHECK_CLOSE(ref.dalq,  1.8148520254931713,  max_d_tol);
    BOOST_CHECK_CLOSE(ref.dflo,  9.0944843574181924,  max_d_tol);

    //Check that different versions of the prod_bph function work
    ADB a = properties->bhp(table_ids, aqua_adb, liquid_adb, vapour_adb, thp_adb, alq_adb);
    double b =properties->bhp(table_ids[0], aqua_d, liquid_d, vapour_d, thp_d, alq_d);
    ADB c = properties->bhp(table_ids, *wells, qs_adb, thp_adb, alq_adb);

    //Check that results are actually equal reference
    BOOST_CHECK_EQUAL(a.value()[0], ref.value);
    BOOST_CHECK_EQUAL(b,            ref.value);
    BOOST_CHECK_EQUAL(c.value()[0], ref.value);

    //Table 2 does not exist.
    std::vector<int> table_ids_wrong(1, 2);
    BOOST_CHECK_THROW(properties->bhp(table_ids_wrong, *wells, qs_adb, thp_adb, alq_adb), std::invalid_argument);
}


/**
 * Test that we can generate some dummy data representing an ND plane,
 * interpolate using ADBs as input, and compare against the analytic solution
 */
BOOST_AUTO_TEST_CASE(ExtrapolatePlaneADB)
{
    fillDataPlane();
    initProperties();

    //Check linear extrapolation (i.e., using values of x, y, etc. outside our interpolant domain)
    double sum = 0.0;
    double reference_sum = 0.0;
    double sad = 0.0; // Sum absolute difference
    double max_d = 0.0; // Maximum difference
    int n=1;
    int o=5;
    for (int i=0; i<=n+o; ++i) {
        const double x = i / static_cast<double>(n);
        for (int j=1; j<=n+o; ++j) {
            const double aqua = -j / static_cast<double>(n);
            for (int k=1; k<=n+o; ++k) {
                const double vapour = -k / static_cast<double>(n);
                for (int l=0; l<=n+o; ++l) {
                    const double u = l / static_cast<double>(n);
                    for (int m=1; m<=n+o; ++m) {
                        const double liquid = -m / static_cast<double>(n);

                        //Temporary variables used to represent independent wells
                        const int num_wells = 5;
                        ADB::V adb_v_x(num_wells);
                        ADB::V adb_v_aqua(num_wells);
                        ADB::V adb_v_vapour(num_wells);
                        ADB::V adb_v_u(num_wells);
                        ADB::V adb_v_liquid(num_wells);
                        table_ids.resize(num_wells);

                        for (unsigned int w=0; w<num_wells; ++w) {
                            table_ids[w] = 1;
                            adb_v_x[w] = x*(w+1);
                            adb_v_aqua[w] = aqua*(w+1);
                            adb_v_vapour[w] = vapour*(w+1);
                            adb_v_u[w] = u*(w+1);
                            adb_v_liquid[w] = liquid*(w+1);
                        }

                        ADB adb_x = ADB::constant(adb_v_x);
                        ADB adb_aqua = ADB::constant(adb_v_aqua);
                        ADB adb_vapour = ADB::constant(adb_v_vapour);
                        ADB adb_u = ADB::constant(adb_v_u);
                        ADB adb_liquid = ADB::constant(adb_v_liquid);

                        ADB bhp = properties->bhp(table_ids, adb_aqua, adb_liquid, adb_vapour, adb_x, adb_u);
                        ADB::V bhp_val = bhp.value();

                        double value = 0.0;
                        double reference = 0.0;
                        for (int w=0; w < num_wells; ++w) {
                            //Find values that should be in table
                            double v = Opm::detail::getFlo(aqua*(w+1), liquid*(w+1), vapour*(w+1), table->getFloType());
                            double y = Opm::detail::getWFR(aqua*(w+1), liquid*(w+1), vapour*(w+1), table->getWFRType());
                            double z = Opm::detail::getGFR(aqua*(w+1), liquid*(w+1), vapour*(w+1), table->getGFRType());

                            reference = x*(w+1) + 2*y + 3*z + 4*u*(w+1) - 5*v;
                            value = bhp_val[w];

                            sum += value;
                            reference_sum += reference;

                            double abs_diff = std::abs(value - reference);

                            sad += std::abs(abs_diff);
                            max_d = std::max(max_d, abs_diff);
                        }
                    }
                }
            }
        }
    }

    BOOST_CHECK_CLOSE(sum, reference_sum, 0.0001);
    BOOST_CHECK_SMALL(max_d, max_d_tol);
    BOOST_CHECK_SMALL(sad, sad_tol);
}



/**
 * Test that we can generate some dummy data representing an ND plane,
 * interpolate using well flow rates and ADBs as input, and compare against the analytic solution
 */
BOOST_AUTO_TEST_CASE(InterpolateADBAndQs)
{
    fillDataPlane();
    initProperties();

    //Create wells
    const int nphases = 3;
    const int nwells  = 5;
    const int nperfs  = 1;
    std::shared_ptr<Wells> wells(create_wells(nphases, nwells, nperfs),
                                destroy_wells);
    int cells = 1;
    for (int i=0; i<nwells; ++i) {
        //Just give the cells a set of different indices
        cells *= 2;

        std::stringstream ss;
        ss << "WELL_" << i;
        const bool ok = add_well(INJECTOR, 0.0, 1, NULL, &cells,
                                 NULL, 0, ss.str().c_str(), true, wells.get());
        BOOST_REQUIRE(ok);
    }

    //Create some artificial flow values for our wells between 0 and 1
    ADB::V qs_v(nphases*nwells);
    for (int j=0; j<nphases; ++j) {
        for (int i=0; i<nwells; ++i) {
            qs_v[j*nwells+i] = -(j*nwells+i) / static_cast<double>(nwells*nphases-1.0);
        }
    }
    ADB qs = ADB::constant(qs_v);

    //Create the THP for each well
    ADB::V thp_v(nwells);
    for (int i=0; i<nwells; ++i) {
        thp_v[i] = (i) / static_cast<double>(nwells-1.0);
    }
    ADB thp = ADB::constant(thp_v);

    //Create the ALQ for each well
    ADB::V alq_v(nwells);
    for (int i=0; i<nwells; ++i) {
        alq_v[i] = 0.0;
    }
    ADB alq = ADB::constant(alq_v);

    //Set which VFP table to use for each well
    table_ids.resize(nwells);
    for (int i=0; i<nwells; ++i) {
        table_ids[i] = 1;
    }

    //Call the bhp function
    ADB::V bhp = properties->bhp(table_ids, *wells, qs, thp, alq).value();

    //Calculate reference
    //First, find the three phases
    std::vector<double> water(nwells);
    std::vector<double> oil(nwells);
    std::vector<double> gas(nwells);
    for (int i=0; i<nwells; ++i) {
        water[i] = qs_v[i];
        oil[i] = qs_v[nwells+i];
        gas[i] = qs_v[2*nwells+i];
    }

    //Compute reference value
    std::vector<double> reference(nwells);
    for (int i=0; i<nwells; ++i) {
        double flo = oil[i];
        double wor = water[i]/oil[i];
        double gor = gas[i]/oil[i];
        reference[i] = thp_v[i] + 2*wor + 3*gor + 4*alq_v[i] - 5*flo;
    }

    //Check that interpolation matches
    BOOST_REQUIRE_EQUAL(bhp.size(), nwells);
    double sad = 0.0;
    double max_d = 0.0;
    for (int i=0; i<nwells; ++i) {
        double value = bhp[i];
        double ref = reference[i];
        double abs_diff = std::abs(value-ref);
        sad += abs_diff;
        max_d = std::max(abs_diff, max_d);
    }

    BOOST_CHECK_SMALL(max_d, max_d_tol);
    BOOST_CHECK_SMALL(sad, sad_tol);
}


BOOST_AUTO_TEST_SUITE_END() // Trivial tests
