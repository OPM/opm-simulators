#include <config.h>

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing
#define BOOST_TEST_MODULE BlackoilStateTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <memory>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>

#include "opm/core/grid/GridManager.hpp"
#include "opm/core/simulator/BlackoilState.hpp"

using namespace Opm;
using namespace std;



BOOST_AUTO_TEST_CASE(EqualsDifferentDeckReturnFalse) {

    const string filename1 = "testBlackoilState1.DATA";
    const string filename2 = "testBlackoilState2.DATA";
    const EclipseGridParser deck1 (filename1);
    const EclipseGridParser deck2 (filename2);

    GridManager gridManager1(deck1);
    const UnstructuredGrid* grid1 = gridManager1.c_grid();
    GridManager gridManager2(deck2);
    const UnstructuredGrid* grid2 = gridManager2.c_grid();
    
    BlackoilState state1;
    state1.init(*grid1, 3);
    
    BlackoilState state2;
    state2.init(*grid2, 3);

    BOOST_CHECK_EQUAL( false , state1.equals(state2) );
}



BOOST_AUTO_TEST_CASE(EqualsDifferentNumPhasesReturnFalse) {

    const string filename = "testBlackoilState1.DATA";
    const EclipseGridParser deck (filename);

    GridManager gridManager(deck);
    const UnstructuredGrid* grid = gridManager.c_grid();
    
    BlackoilState state1;
    state1.init(*grid, 3);
    
    BlackoilState state2;
    state2.init(*grid, 2);

    BOOST_CHECK_EQUAL( false , state1.equals(state2) );
}



BOOST_AUTO_TEST_CASE(EqualsNumericalDifferenceReturnFalse) {

    const string filename = "testBlackoilState1.DATA";
    const EclipseGridParser deck (filename);

    GridManager gridManager(deck);
    const UnstructuredGrid* grid = gridManager.c_grid();
    
    BlackoilState state1;
    state1.init(*grid, 3);
    
    BlackoilState state2;
    state2.init(*grid, 3);

    BOOST_CHECK_EQUAL( true , state1.equals(state2) );
    {
        std::vector<double>& p1 = state1.pressure();
        std::vector<double>& p2 = state2.pressure();
        p1[0] = p1[0] * 2 + 1;
        
        BOOST_CHECK_EQUAL( false , state1.equals(state2) );
        p1[0] = p2[0];
        BOOST_CHECK_EQUAL( true , state1.equals(state2) );
    }
    {
        std::vector<double>& gor1 = state1.gasoilratio();
        std::vector<double>& gor2 = state2.gasoilratio();
        gor1[0] = gor1[0] * 2 + 1;
        
        BOOST_CHECK_EQUAL( false , state1.equals(state2) );
        gor1[0] = gor2[0];
        BOOST_CHECK_EQUAL( true , state1.equals(state2) );
    }
    {
        std::vector<double>& p1 = state1.facepressure();
        std::vector<double>& p2 = state2.facepressure();
        p1[0] = p1[0] * 2 + 1;
        
        BOOST_CHECK_EQUAL( false , state1.equals(state2) );
        p1[0] = p2[0];
        BOOST_CHECK_EQUAL( true , state1.equals(state2) );
    }

    {
        std::vector<double>& f1 = state1.faceflux();
        std::vector<double>& f2 = state2.faceflux();
        if (f1.size() > 0 ) {
            f1[0] = f1[0] * 2 + 1;
            
            BOOST_CHECK_EQUAL( false , state1.equals(state2) );
            f1[0] = f2[0];
            BOOST_CHECK_EQUAL( true , state1.equals(state2) );
        }
    }
    {
        std::vector<double>& sv1 = state1.surfacevol();
        std::vector<double>& sv2 = state2.surfacevol();
        if (sv1.size() > 0) {
            sv1[0] = sv1[0] * 2 + 1;
            
            BOOST_CHECK_EQUAL( false , state1.equals(state2) );
            sv1[0] = sv2[0];
            BOOST_CHECK_EQUAL( true , state1.equals(state2) );
        }
    }
    {
        std::vector<double>& sat1 = state1.saturation();
        std::vector<double>& sat2 = state2.saturation();
        sat1[0] = sat1[0] * 2 + 1;
        
        BOOST_CHECK_EQUAL( false , state1.equals(state2) );
        sat1[0] = sat2[0];
        BOOST_CHECK_EQUAL( true , state1.equals(state2) );
    }
}
