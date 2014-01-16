/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
#ifndef OPM_IMMISCIBLE_FLUID_STATE_HPP
#define OPM_IMMISCIBLE_FLUID_STATE_HPP

#include "ModularFluidState.hpp"

#include <opm/material/Valgrind.hpp>

#include <algorithm>

#include <string.h>

namespace Opm {

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
template <class Scalar, class FluidSystem, bool storeEnthalpy=true>
class ImmiscibleFluidState;

// specialization for the enthalpy enabled case
template <class Scalar, class FluidSystem>
class ImmiscibleFluidState<Scalar, FluidSystem, true>
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateImmiscibleCompositionModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitEnthalpyModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> > >
{
public:
    ImmiscibleFluidState()
    {}

    ImmiscibleFluidState(const ImmiscibleFluidState &fs)
    { memcpy(this, &fs, sizeof(fs)); }
};

// specialization for the enthalpy disabled case
template <class Scalar, class FluidSystem>
class ImmiscibleFluidState<Scalar, FluidSystem, false>
    : public ModularFluidState<Scalar,
                               FluidSystem,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateImmiscibleCompositionModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateNullEnthalpyModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> > >
{
public:
    ImmiscibleFluidState()
    {}

    ImmiscibleFluidState(const ImmiscibleFluidState &fs)
    { memcpy(this, &fs, sizeof(fs)); }
};
} // namespace Opm

#endif
