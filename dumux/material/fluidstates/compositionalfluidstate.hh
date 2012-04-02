// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
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
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_COMPOSITIONAL_FLUID_STATE_HH
#define DUMUX_COMPOSITIONAL_FLUID_STATE_HH

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
class CompositionalFluidState
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitCompositionModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> >,
                               FluidStateExplicitEnthalpyModule<Scalar, FluidSystem, CompositionalFluidState<Scalar, FluidSystem> > >                            
{
};

} // end namepace Dumux

#endif
