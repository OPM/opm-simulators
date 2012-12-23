// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 *        multi-phase, multi-component fluid system _not_ assuming
 *        thermodynamic equilibrium.
 */
#ifndef EWOMS_NON_EQUILIBRIUM_FLUID_STATE_HH
#define EWOMS_NON_EQUILIBRIUM_FLUID_STATE_HH

#include "modularfluidstate.hh"

#include <ewoms/common/valgrind.hh>
#include <algorithm>

namespace Ewoms {

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system _not_ assuming
 *        thermodynamic equilibrium.
 */
template <class Scalar, class FluidSystem, bool storeEnthalpy=true>
class NonEquilibriumFluidState;


// specialization for the enthalpy enabled case
template <class Scalar, class FluidSystem>
class NonEquilibriumFluidState<Scalar, FluidSystem, true>
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitTemperatureModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitCompositionModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitEnthalpyModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, true> > >
{
};

// specialization for the enthalpy disabled case
template <class Scalar, class FluidSystem>
class NonEquilibriumFluidState<Scalar, FluidSystem, false>
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitTemperatureModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitCompositionModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> >,
                               FluidStateNullEnthalpyModule<Scalar, FluidSystem, NonEquilibriumFluidState<Scalar, FluidSystem, false> > >
{
};

} // end namespace

#endif
