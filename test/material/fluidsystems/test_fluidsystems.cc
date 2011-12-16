// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
#include <test/boxmodels/1p2c/interstitialfluidtrailfluidsystem.hh>
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <appl/lecture/msm/1p2cvs2p/watercontaminantfluidsystem.hh>

int main()
{
    typedef double Scalar;
    typedef Dumux::H2O<Scalar> H2O;
    typedef Dumux::N2<Scalar> N2;

    typedef Dumux::LiquidPhase<Scalar, H2O> Liquid;
    typedef Dumux::GasPhase<Scalar, N2> Gas;

    // H2O -- N2
    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/false> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    {   typedef Dumux::FluidSystems::H2ON2<Scalar, /*enableComplexRelations=*/true> FluidSystem;
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

    // water -- contaminant
    {   typedef Dumux::FluidSystems::WaterContaminant<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    // interstitial fluid -- TRAIL
    {   typedef Dumux::FluidSystems::InterstitialFluidTrail<Scalar> FluidSystem;
        checkFluidSystem<Scalar, FluidSystem>(); }

    return 0;
}
