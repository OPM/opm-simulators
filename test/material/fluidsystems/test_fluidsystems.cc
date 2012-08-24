// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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

// include the property system just to make sure that all fluid system
// type tag adapter behave nicely together
#include <dumux/common/propertysystem.hh>

#include "checkfluidsystem.hh"

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/blackoilfluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairxylenefluidsystem.hh>

// include all fluid states
#include <dumux/material/fluidstates/pressureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/saturationoverlayfluidstate.hh>
#include <dumux/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/fluidstates/nonequilibriumfluidstate.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

int main()
{
    typedef double Scalar;
    typedef Dumux::H2O<Scalar> H2O;
    typedef Dumux::N2<Scalar> N2;

    typedef Dumux::LiquidPhase<Scalar, H2O> Liquid;
    typedef Dumux::GasPhase<Scalar, N2> Gas;

    // check all fluid states
    {
        typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;

        // CompositionalFluidState
        {   Dumux::CompositionalFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // NonEquilibriumFluidState
        {   Dumux::NonEquilibriumFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        // ImmiscibleFluidState
        {   Dumux::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            checkFluidState<Scalar>(fs); }

        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> BaseFluidState;
        BaseFluidState baseFs;

        // TemperatureOverlayFluidState
        {   Dumux::TemperatureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // PressureOverlayFluidState
        {   Dumux::PressureOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }

        // SaturationOverlayFluidState
        {   Dumux::SaturationOverlayFluidState<Scalar, BaseFluidState> fs(baseFs);
            checkFluidState<Scalar>(fs); }
    }

    // black-oil
    {   typedef Dumux::FluidSystems::BlackOil<Scalar> FluidSystem;
        if (false) checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2
    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- N2 -- liquid phase
    {   typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, /*enableComplexRelations=*/true> FluidSystem;
         checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air
    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::SimpleH2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=false;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::H2O<Scalar> H2O;
        const bool enableComplexRelations=true;
        typedef Dumux::FluidSystems::H2OAir<Scalar, H2O, enableComplexRelations> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Mesitylene
    {   typedef Dumux::FluidSystems::H2OAirMesitylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // H2O -- Air -- Xylene
    {   typedef Dumux::FluidSystems::H2OAirXylene<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 2p-immiscible
    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Liquid, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {  typedef Dumux::FluidSystems::TwoPImmiscible<Scalar, Gas, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // 1p
    {   typedef Dumux::FluidSystems::OneP<Scalar, Liquid> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::OneP<Scalar, Gas> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    return 0;
}
