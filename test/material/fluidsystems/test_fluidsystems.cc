// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Katherina Baber                                   *
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

#include <dune/common/mpihelper.hh>

// include all fluid systems in ewoms-stable
#include <ewoms/material/fluidsystems/1pfluidsystem.hh>
#include <ewoms/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>
#include <ewoms/material/fluidsystems/brineco2fluidsystem.hh>
#include <ewoms/material/fluidsystems/h2on2fluidsystem.hh>
#include <ewoms/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairfluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairxylenefluidsystem.hh>

// include all fluid states
#include <ewoms/material/fluidstates/pressureoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/saturationoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/fluidstates/nonequilibriumfluidstate.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>

// include the tables for CO2 which are delivered with eWoms by
// default
#include <ewoms/common/statictabulated2dfunction.hh>
namespace Ewoms {
namespace FluidSystemsTest {
#include <ewoms/material/components/co2tables.inc>
} }

int main(int argc, char **argv)
{
    typedef double Scalar;
    typedef Ewoms::H2O<Scalar> H2O;
    typedef Ewoms::N2<Scalar> N2;

    typedef Ewoms::LiquidPhase<Scalar, H2O> Liquid;
    typedef Ewoms::GasPhase<Scalar, N2> Gas;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // check all fluid states
    {
        typedef Ewoms::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;

        // CompositionalFluidState
        {   Ewoms::CompositionalFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // NonEquilibriumFluidState
        {   Ewoms::NonEquilibriumFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // ImmiscibleFluidState
        {   Ewoms::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        typedef Ewoms::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
        BaseFluidState baseFs;

        // TemperatureOverlayFluidState
        {   Ewoms::TemperatureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // PressureOverlayFluidState
        {   Ewoms::PressureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // SaturationOverlayFluidState
        {   Ewoms::SaturationOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }
    }

    // black-oil
    {   typedef Ewoms::FluidSystems::BlackOil<Scalar> FluidSystem;
        if (false) checkFluidSystem<Scalar, FluidSystem>(); }

    // Brine -- CO2
    {   typedef Ewoms::FluidSystems::BrineCO2<Scalar, Ewoms::FluidSystemsTest::CO2Tables> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   typedef Ewoms::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- liquid phase
    {   typedef Ewoms::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/true> FluidSystem;
         checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air
    {   typedef Ewoms::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Ewoms::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Ewoms::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::H2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Ewoms::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::H2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Ewoms::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Ewoms::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   typedef Ewoms::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   typedef Ewoms::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {  typedef Ewoms::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 1p
    {   typedef Ewoms::FluidSystems::OneP<Scalar, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Ewoms::FluidSystems::OneP<Scalar, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    return 0;
}
