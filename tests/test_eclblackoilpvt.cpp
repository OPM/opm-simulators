// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \brief This is the unit test for the black oil PVT classes
 *
 * This test requires the presence of opm-parser.
 */
#include "config.h"

#if !HAVE_OPM_PARSER
#error "The test for the black oil PVT classes requires the opm-parser module"
#endif

#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>

#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

// values of strings based on the first SPE1 test case of opm-data.  note that in the
// real world it does not make much sense to specify a fluid phase using more than a
// single keyword, but for a unit test, this saves a lot of boiler-plate code.
static const char* deckString1 =
    "RUNSPEC\n"
    "\n"
    "DIMENS\n"
    "   10 10 3 /\n"
    "\n"
    "TABDIMS\n"
    " * 2 /\n"
    "\n"
    "OIL\n"
    "GAS\n"
    "WATER\n"
    "\n"
    "DISGAS\n"
    "\n"
    "METRIC\n"
    "\n"
    "GRID\n"
    "\n"
    "DX\n"
    "   	300*1000 /\n"
    "DY\n"
    "	300*1000 /\n"
    "DZ\n"
    "	100*20 100*30 100*50 /\n"
    "\n"
    "TOPS\n"
    "	100*1234 /\n"
    "\n"
    "PROPS\n"
    "\n"
    "DENSITY\n"
    "      	100.0 500.0 0.1 /\n"
    "      	200.0 600.0 0.2 /\n"
    "\n"
    "PVTW\n"
    " 	1.0  1.1 1e-6 1.1 2.0e-9 /\n"
    " 	2.0  1.2 1e-7 1.2 3.0e-9 /\n"
    "\n"
    "PVCDO\n"
    " 	1.0  1.1 1e-6 1.1 2.0e-9 /\n"
    " 	2.0  1.2 1e-7 1.2 3.0e-9 /\n"
    "\n"
    "PVDG\n"
    "1.0	1.0	10.0\n"
    "2.0	*	*\n"
    "3.0	1e-10	30.0 /\n"
    "\n"
    "4.0	1.0	40.0\n"
    "5.0	*	*\n"
    "6.0	1e-10	60.0 /\n"
    "\n"
    "PVTG\n"
    "\n"
    "-- PVT region 2 --\n"
    "-- PRESSURE       RV        BG     VISCOSITY\n"
    "     1.00     1.1e-3       1.1      0.01\n"
    "              1.0e-3       1.15     0.005 /\n"
    "\n"
    "    500.00    0.9e-3       1.2     0.02\n"
    "              0.8e-3       1.25    0.015 /\n"
    "/\n"
    "\n"
    "-- PVT region 2 --\n"
    "-- PRESSURE       RV        BG     VISCOSITY\n"
    "     2.00     2.1e-3       2.1      0.02\n"
    "              2.0e-3       2.15     0.015 /\n"
    "\n"
    "    502.00    1.2e-3       2.2     2.02\n"
    "              1.1e-3       2.25    2.015 /\n"
    "/\n"
    "\n";

int main()
{
    typedef double Scalar;

    Opm::Parser parser;
    Opm::ParseMode parseMode;

    const auto deck = parser.parseString(deckString1, parseMode);
    const auto eclState = std::make_shared<Opm::EclipseState>(deck, parseMode);
    const auto eclGrid = eclState->getEclipseGrid();

    const auto pvtwKeyword = deck->getKeyword("PVTW");
    size_t numPvtRegions = pvtwKeyword->size();

    if (numPvtRegions != 2)
        OPM_THROW(std::logic_error,
                  "The number of PVT regions of the test deck must be 2. (is "
                  << numPvtRegions << ")");

    //////////
    // constant compressibility water
    //////////
    Opm::ConstantCompressibilityWaterPvt<Scalar> constCompWaterPvt;
    constCompWaterPvt.initFromDeck(deck, eclState);

    // make sure that the values at the reference points are the ones specified in the
    // deck.
    Scalar refTmp, tmp;

    refTmp = 1.1e-3; // the deck value is given in cP, while the SI units use Pa s...
    tmp = constCompWaterPvt.viscosity(/*regionIdx=*/0,
                                      /*temperature=*/273.15 + 20.0,
                                      /*pressure=*/1e5);
    if (std::abs(tmp - refTmp)  > 1e-30)
        OPM_THROW(std::logic_error,
                  "The reference water viscosity at region 0 is supposed to be " << refTmp
                  << ". (is " << tmp << ")");

    refTmp = 1.2e-3;
    tmp = constCompWaterPvt.viscosity(/*regionIdx=*/1,
                                      /*temperature=*/273.15 + 20.0,
                                      /*pressure=*/2e5);
    if (std::abs(tmp - refTmp)  > 1e-30)
        OPM_THROW(std::logic_error,
                  "The reference water viscosity at region 1 is supposed to be " << refTmp
                  << ". (is " << tmp << ")");

    refTmp = 500/1.1;
    tmp = constCompWaterPvt.density(/*regionIdx=*/0,
                                    /*temperature=*/273.15 + 20.0,
                                    /*pressure=*/1e5);
    if (std::abs(tmp - refTmp)  > 5e-14)
        OPM_THROW(std::logic_error,
                  "The reference water density at region 0 is supposed to be " << refTmp
                  << ". (is " << tmp << ")");


    refTmp = 600/1.2;
    tmp = constCompWaterPvt.density(/*regionIdx=*/1,
                                    /*temperature=*/273.15 + 20.0,
                                    /*pressure=*/2e5);
    if (std::abs(tmp - refTmp)  > 5e-14)
        OPM_THROW(std::logic_error,
                  "The reference water density at region 1 is supposed to be " << refTmp
                  << ". (is " << tmp << ")");

    //////////
    // the gas and oil PVT classes.
    //
    // TODO: check the results
    //////////
    Opm::GasPvtMultiplexer<Scalar> gasPvt;
    Opm::OilPvtMultiplexer<Scalar> oilPvt;
    Opm::WaterPvtMultiplexer<Scalar> waterPvt;

    gasPvt.initFromDeck(deck, eclState);
    oilPvt.initFromDeck(deck, eclState);
    waterPvt.initFromDeck(deck, eclState);

    gasPvt.initEnd(&oilPvt);
    oilPvt.initEnd(&gasPvt);
    waterPvt.initEnd();

    // make sure that the BlackOil fluid system's initFromDeck() method compiles.
    typedef Opm::FluidSystems::BlackOil<Scalar> BlackOilFluidSystem;
    BlackOilFluidSystem::initFromDeck(deck, eclState);

    return 0;
}
