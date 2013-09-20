// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This test makes sure that the programming interface is
 *        observed by all fluid systems
 */
#include "config.h"

#include "checkfluidsystem.hh"

// include all fluid systems in opm-material
#include <opm/material/fluidsystems/1pfluidsystem.hh>
#include <opm/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <opm/material/fluidsystems/blackoilfluidsystem.hh>
#include <opm/material/fluidsystems/brineco2fluidsystem.hh>
#include <opm/material/fluidsystems/h2on2fluidsystem.hh>
#include <opm/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <opm/material/fluidsystems/h2oairfluidsystem.hh>
#include <opm/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <opm/material/fluidsystems/h2oairxylenefluidsystem.hh>

// include all fluid states
#include <opm/material/fluidstates/pressureoverlayfluidstate.hh>
#include <opm/material/fluidstates/saturationoverlayfluidstate.hh>
#include <opm/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <opm/material/fluidstates/compositionalfluidstate.hh>
#include <opm/material/fluidstates/nonequilibriumfluidstate.hh>
#include <opm/material/fluidstates/immisciblefluidstate.hh>

// include the tables for CO2 which are delivered with opm-material by
// default
#include <opm/material/statictabulated2dfunction.hh>
namespace Opm {
namespace FluidSystemsTest {
#include <opm/material/components/co2tables.inc>
} }

int main(int argc, char **argv)
{
    typedef double Scalar;
    typedef Opm::H2O<Scalar> H2O;
    typedef Opm::N2<Scalar> N2;

    typedef Opm::LiquidPhase<Scalar, H2O> Liquid;
    typedef Opm::GasPhase<Scalar, N2> Gas;

    // check all fluid states
    {
        typedef Opm::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;

        // CompositionalFluidState
        {   Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // NonEquilibriumFluidState
        {   Opm::NonEquilibriumFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // ImmiscibleFluidState
        {   Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        typedef Opm::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
        BaseFluidState baseFs;

        // TemperatureOverlayFluidState
        {   Opm::TemperatureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // PressureOverlayFluidState
        {   Opm::PressureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // SaturationOverlayFluidState
        {   Opm::SaturationOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }
    }

    // black-oil
    {   typedef Opm::FluidSystems::BlackOil<Scalar> FluidSystem;
        if (false) checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2
    {   typedef Opm::FluidSystems::BrineCO2<Scalar, Opm::FluidSystemsTest::CO2Tables> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   typedef Opm::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- liquid phase
    {   typedef Opm::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/true> FluidSystem;
         checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air
    {   typedef Opm::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Opm::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Opm::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::H2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Opm::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::H2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Opm::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Opm::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   typedef Opm::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   typedef Opm::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {  typedef Opm::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 1p
    {   typedef Opm::FluidSystems::OneP<Scalar, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Opm::FluidSystems::OneP<Scalar, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    return 0;
}
