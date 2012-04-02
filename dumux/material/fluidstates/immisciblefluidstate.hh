// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_IMMISCIBLE_FLUID_STATE_HH
#define DUMUX_IMMISCIBLE_FLUID_STATE_HH

#include "modularfluidstate.hh"

#include <dumux/common/valgrind.hh>
#include <algorithm>

namespace Dumux
{
/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
template <class Scalar, class FluidSystem>
class ImmiscibleFluidState
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateImmiscibleCompositionModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateImmiscibleFugacityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitEnthalpyModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem> > >                            
{
};
} // end namespace

#endif
