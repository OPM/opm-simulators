/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
*/

#include "config.h"

/* --- Boost.Test boilerplate --- */
#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE UnitsTest
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

/* --- our own headers --- */

#include <opm/core/simulator/initStateEquil.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/cart_grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/core/pressure/msmfem/partition.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE (GwsegStandard)
{
    // This is the basic (no eps and hysteris) version of 
    // the Gwseg model.
    
    //std::cout << "==================================== GwsegStandard ====================================" << std::endl;

    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncStandard.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    for (int i=0; i<n; ++i) {
      cells[i] = 0;
      s[i*np+wpos] = i*0.1;
      s[i*np+opos] = 1.0-s[i*np+wpos];
      s[i*np+gpos] = 0.0;
    }
    
    props.relperm(n, s, cells, kr, dkrds);
    
    double krw[11] = {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7};
    double kro[11] = {1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0};
    double DkrwDsw[11] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    double DkroDsw[11] = {-2.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0};
    double DkroDsg[11] = {-5.0, -5.0, -3.0, -2.0, -0.66666666666666741, -0.75, -0.8, 
                          -0.83333333333333237, 0.14285714285714296, 0.0, 0.0};
    
    const double reltol = 1.0e-6;
    for (int i=0; i<n; ++i) {
      BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[i], reltol);
      BOOST_CHECK_CLOSE(kr[i*np+opos], kro[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[i], reltol);
    }

/*    
    std::cout << std::setw(12) << "sw";
    std::cout << std::setw(12) << "so";
    std::cout << std::setw(12) << "sg";
    std::cout << std::setw(12) << "krw";
    std::cout << std::setw(12) << "kro";
    std::cout << std::setw(12) << "krg";
    std::cout << std::setw(12) << "DkrwDsw";
    std::cout << std::setw(12) << "DkroDsw";
    std::cout << std::setw(12) << "DkrgDsw";
    std::cout << std::setw(12) << "DkrwDso";
    std::cout << std::setw(12) << "DkroDso";
    std::cout << std::setw(12) << "DkrgDso";
    std::cout << std::setw(12) << "DkrwDsg";
    std::cout << std::setw(12) << "DkroDsg";
    std::cout << std::setw(12) << "DkrgDsg";
    std::cout << std::endl;
    for (int i=0; i<n; ++i) {
      std::cout << std::setw(12) << s[i*np+wpos] << std::setw(12) << s[i*np+opos] << std::setw(12) << s[i*np+gpos]
                << std::setw(12) << kr[i*np+wpos] << std::setw(12) << kr[i*np+opos] << std::setw(12) << kr[i*np+gpos];
      for (int j=0; j<np*np; ++j) {
        std::cout << std::setw(12) << dkrds[i*np*np+j];
      }
      std::cout << std::endl;
    }
*/
}

BOOST_AUTO_TEST_CASE (GwsegEPSBase)
{
    // This is the eps (but no hysteris) version of the Gwseg model.
    // However, only default scaling parameters, i.e no scaling.
    
    //std::cout << "==================================== GwsegEPSBase ====================================" << std::endl;
    
    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncEPSBase.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    for (int i=0; i<n; ++i) {
      cells[i] = 0;
      s[i*np+wpos] = i*0.1;
      s[i*np+opos] = 1.0-s[i*np+wpos];
      s[i*np+gpos] = 0.0;
    }
    
    props.relperm(n, s, cells, kr, dkrds);
 
    double krw[11] = {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7};
    double kro[11] = {1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0};
    double DkrwDsw[11] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    double DkroDsw[11] = {-2.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0};
    double DkroDsg[11] = {-5.0, -5.0, -3.0, -2.0,-0.66666666666666741, -0.75, -0.8, -0.83333333333333237, 0.14285714285714296, 0.0, 0.0};
   
    const double reltol = 1.0e-6;
    for (int i=0; i<n; ++i) {
      BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[i], reltol);
      BOOST_CHECK_CLOSE(kr[i*np+opos], kro[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[i], reltol);
    }
 
/*
    std::cout << std::setw(12) << "sw";
    std::cout << std::setw(12) << "so";
    std::cout << std::setw(12) << "sg";
    std::cout << std::setw(12) << "krw";
    std::cout << std::setw(12) << "kro";
    std::cout << std::setw(12) << "krg";
    std::cout << std::setw(12) << "DkrwDsw";
    std::cout << std::setw(12) << "DkroDsw";
    std::cout << std::setw(12) << "DkrgDsw";
    std::cout << std::setw(12) << "DkrwDso";
    std::cout << std::setw(12) << "DkroDso";
    std::cout << std::setw(12) << "DkrgDso";
    std::cout << std::setw(12) << "DkrwDsg";
    std::cout << std::setw(12) << "DkroDsg";
    std::cout << std::setw(12) << "DkrgDsg";
    std::cout << std::endl;
    for (int i=0; i<n; ++i) {
      std::cout << std::setw(12) << s[i*np+wpos] << std::setw(12) << s[i*np+opos] << std::setw(12) << s[i*np+gpos]
                << std::setw(12) << kr[i*np+wpos] << std::setw(12) << kr[i*np+opos] << std::setw(12) << kr[i*np+gpos];
      for (int j=0; j<np*np; ++j) {
        std::cout << std::setw(12) << dkrds[i*np*np+j];
      }
      std::cout << std::endl;
    }
*/
}


BOOST_AUTO_TEST_CASE (GwsegEPS_A)    
{
    // This is the eps (but no hysteris) version of the Gwseg model.
    // Scaling parameters from keyword SWL etc.
    
    //std::cout << "==================================== GwsegEPS_A ====================================" << std::endl;
    
    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncEPS_A.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    const int ncell = 8;
 
    double krw[ncell][n] = {{0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7},
                            {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7}};
    double kro[ncell][n] = {{1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0}};
    double DkrwDsw[ncell][n] = {{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0},
                                {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0}};
    double DkroDsw[ncell][n] = {{0.0, 0.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0}};
    double DkroDsg[ncell][n] = {{-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0}}; 
    
    for (int icell=0; icell<ncell; ++icell) {
      for (int i=0; i<n; ++i) {
        cells[i] = icell;
        s[i*np+wpos] = i*0.1;
        s[i*np+opos] = 1.0-s[i*np+wpos];
        s[i*np+gpos] = 0.0;
      }
      
      props.relperm(n, s, cells, kr, dkrds);
      
      const double reltol = 1.0e-3;
      for (int i=0; i<n; ++i) {
        BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[icell][i], reltol);
        BOOST_CHECK_CLOSE(kr[i*np+opos], kro[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[icell][i], reltol);
      }
      
/*      
      std::cout << std::setw(12) << "sw: ";
      for (int i=0; i<n; ++i) std::cout << s[i*np+wpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "so: ";
      for (int i=0; i<n; ++i) std::cout << s[i*np+opos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "sg: ";
      for (int i=0; i<n; ++i) std::cout << s[i*np+gpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "krw: ";
      for (int i=0; i<n; ++i) std::cout << kr[i*np+wpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "kro: ";
      for (int i=0; i<n; ++i) std::cout << kr[i*np+opos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "krg: ";
      for (int i=0; i<n; ++i) std::cout << kr[i*np+gpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrwDsw: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+wpos*np+wpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkroDsw: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+wpos*np+opos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrgDsw: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+wpos*np+gpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrwDso: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+opos*np+wpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkroDso: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+opos*np+opos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrgDso: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+opos*np+gpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrwDsg: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+gpos*np+wpos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkroDsg: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+gpos*np+opos] << ", ";
      std::cout << std::endl;
      std::cout << std::setw(12) << "DkrgDsg: ";
      for (int i=0; i<n; ++i) std::cout << dkrds[i*np*np+gpos*np+gpos] << ", ";
      std::cout << std::endl;
      std::cout << std::endl;
*/
    
    }

}

BOOST_AUTO_TEST_CASE (GwsegEPS_B)    
{
    // This is the eps (but no hysteris) version of the Gwseg model.
    // Scaling parameters from ENPTVD table.
    
    //std::cout << "==================================== GwsegEPS_B ====================================" << std::endl;
 /*   
    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncEPS_B.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    const int ncell = 8;
 
    double krw[ncell][n] = {{0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7},
                            {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7}};
    double kro[ncell][n] = {{1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0}};
    double DkrwDsw[ncell][n] = {{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0},
                                {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0}};
    double DkroDsw[ncell][n] = {{0.0, 0.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0}};
    double DkroDsg[ncell][n] = {{-2.32831e-10, -2.32831e-10, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-2.32831e-10, -2.32831e-10, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-2.32831e-10, -2.32831e-10, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-2.32831e-10, -2.32831e-10, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-2.32831e-10, -2.32831e-10, -2.32831e-10, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-2.32831e-10, -2.32831e-10, -2.32831e-10, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-2.32831e-10, -2.32831e-10, -2.32831e-10, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-2.32831e-10, -2.32831e-10, -2.32831e-10, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0}}; 
    
    for (int icell=0; icell<ncell; ++icell) {
      for (int i=0; i<n; ++i) {
        cells[i] = icell;
        s[i*np+wpos] = i*0.1;
        s[i*np+opos] = 1.0-s[i*np+wpos];
        s[i*np+gpos] = 0.0;
      }
      
      props.relperm(n, s, cells, kr, dkrds);
      
      const double reltol = 1.0e-3;
      for (int i=0; i<n; ++i) {
        BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[icell][i], reltol);
        BOOST_CHECK_CLOSE(kr[i*np+opos], kro[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[icell][i], reltol);
      }
  
    }
    */

}

BOOST_AUTO_TEST_CASE (GwsegEPS_C)    
{
    // This is the eps (but no hysteris) version of the Gwseg model.
    // Scaling parameters given the "norne-way", i.e EQUALS, COPY, ADD, MULTIPLY.
    
    //std::cout << "==================================== GwsegEPS_C ====================================" << std::endl;
    
    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncEPS_C.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    const int ncell = 8;
 
    double krw[ncell][n] = {{0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7},
                            {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7},
                            {0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.14, 0.28, 0.42, 0.56, 0.7, 0.7},
                            {0, 0, 0, 0, 0, 0.233333, 0.466667, 0.7, 0.7, 0.7, 0.7}};
    double kro[ncell][n] = {{1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0},
                            {1, 1, 1, 0.766667, 0.533333, 0.35, 0.233333, 0.116667, 0, 0, 0}};
    double DkrwDsw[ncell][n] = {{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0},
                                {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                {0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 1.4, 1.4, 1.4, 1.4, 0, 0},
                                {0, 0, 0, 0, 0, 2.33333, 2.33333, 0, 0, 0, 0}};
    double DkroDsw[ncell][n] = {{0.0, 0.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, -2, -2, -1, -1, -1, -1, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0},
                                {0, 0, 0, -2.33333, -2.33333, -1.16667, -1.16667, -1.16667, 0, 0, 0}};
    double DkroDsg[ncell][n] = {{-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3, -2, -0.666667, -0.75, -0.8, -0.833333, 0.142857, 0, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0},
                                {-3.0, -3.0, -3.0, -3.14286, -2.14286, -0.809524, -0.892857, -0.942857, 0.190476, 3.17207e-16, 0}}; 
    
    for (int icell=0; icell<ncell; ++icell) {
      for (int i=0; i<n; ++i) {
        cells[i] = icell;
        s[i*np+wpos] = i*0.1;
        s[i*np+opos] = 1.0-s[i*np+wpos];
        s[i*np+gpos] = 0.0;
      }
      
      props.relperm(n, s, cells, kr, dkrds);
      
      const double reltol = 1.0e-3;
      for (int i=0; i<n; ++i) {
        BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[icell][i], reltol);
        BOOST_CHECK_CLOSE(kr[i*np+opos], kro[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[icell][i], reltol);
        BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[icell][i], reltol);
      }
    
    }

}

BOOST_AUTO_TEST_CASE (GwsegEPS_D)
{
    // This is the eps and hysteris version of the Gwseg model.
    // However, only default scaling parameters, i.e no scaling.
    
    //std::cout << "==================================== GwsegEPS_D ====================================" << std::endl;

    Opm::parameter::ParameterGroup param;
    
    Opm::GridManager gm(1, 1, 10, 1.0, 1.0, 5.0);
    const UnstructuredGrid& grid = *(gm.c_grid());
    Opm::ParserPtr parser(new Opm::Parser() );
    Opm::DeckConstPtr deck = parser->parseFile("satfuncEPS_D.DATA");
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));
    Opm::BlackoilPropertiesFromDeck props(deck, eclipseState, grid, param, false);
    
    const int np = props.numPhases();
    const int wpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Aqua];
    const int opos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Liquid];
    const int gpos = props.phaseUsage().phase_pos[Opm::BlackoilPhases::Vapour];
    
    BOOST_REQUIRE(np == 3);
    BOOST_REQUIRE(wpos == 0);
    BOOST_REQUIRE(opos == 1);
    BOOST_REQUIRE(gpos == 2);
    
    const int n=11;
    double s[n*np];
    int cells[n];
    double kr[n*np];
    double dkrds[n*np*np];
    
    for (int i=0; i<n; ++i) {
      cells[i] = 0;
      s[i*np+wpos] = i*0.1;
      s[i*np+opos] = 1.0-s[i*np+wpos];
      s[i*np+gpos] = 0.0;
    }
    
    props.relperm(n, s, cells, kr, dkrds);
    
    double krw[11] = {0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7};
    double kro[11] = {1.0, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0};
    double DkrwDsw[11] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    double DkroDsw[11] = {-2.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0};
    double DkroDsg[11] = {-5.0, -5.0, -3.0, -2.0, -0.66666666666666741, -0.75, -0.8, 
                          -0.83333333333333237, 0.14285714285714296, 0.0, 0.0};
    
    const double reltol = 1.0e-6;
    for (int i=0; i<n; ++i) {
      BOOST_CHECK_CLOSE(kr[i*np+wpos], krw[i], reltol);
      BOOST_CHECK_CLOSE(kr[i*np+opos], kro[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+wpos], DkrwDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+opos], DkroDsw[i], reltol);
      BOOST_CHECK_CLOSE(dkrds[i*np*np+np*gpos+opos], DkroDsg[i], reltol);
    }

/*    
    std::cout << std::setw(12) << "sw";
    std::cout << std::setw(12) << "so";
    std::cout << std::setw(12) << "sg";
    std::cout << std::setw(12) << "krw";
    std::cout << std::setw(12) << "kro";
    std::cout << std::setw(12) << "krg";
    std::cout << std::setw(12) << "DkrwDsw";
    std::cout << std::setw(12) << "DkroDsw";
    std::cout << std::setw(12) << "DkrgDsw";
    std::cout << std::setw(12) << "DkrwDso";
    std::cout << std::setw(12) << "DkroDso";
    std::cout << std::setw(12) << "DkrgDso";
    std::cout << std::setw(12) << "DkrwDsg";
    std::cout << std::setw(12) << "DkroDsg";
    std::cout << std::setw(12) << "DkrgDsg";
    std::cout << std::endl;
    for (int i=0; i<n; ++i) {
      std::cout << std::setw(12) << s[i*np+wpos] << std::setw(12) << s[i*np+opos] << std::setw(12) << s[i*np+gpos]
                << std::setw(12) << kr[i*np+wpos] << std::setw(12) << kr[i*np+opos] << std::setw(12) << kr[i*np+gpos];
      for (int j=0; j<np*np; ++j) {
        std::cout << std::setw(12) << dkrds[i*np*np+j];
      }
      std::cout << std::endl;
    }
*/
}

BOOST_AUTO_TEST_SUITE_END()
