// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::SimpleModularFluidState
 */
#ifndef OPM_SIMPLE_MODULAR_FLUID_STATE_HPP
#define OPM_SIMPLE_MODULAR_FLUID_STATE_HPP

#include "FluidStatePressureModules.hpp"
#include "FluidStateTemperatureModules.hpp"
#include "FluidStateCompositionModules.hpp"
#include "FluidStateFugacityModules.hpp"
#include "FluidStateSaturationModules.hpp"
#include "FluidStateDensityModules.hpp"
#include "FluidStateViscosityModules.hpp"
#include "FluidStateEnthalpyModules.hpp"
#include "ModularFluidState.hpp"

#include <type_traits>

namespace Opm {
// this macro is a small hack to prevent death-through verbosity
#define OPM_SMFS SimpleModularFluidState<ScalarT, \
                                         numPhasesV, \
                                         numComponentsV,    \
                                         FluidSystem,       \
                                         storePressure,     \
                                         storeTemperature,  \
                                         storeComposition,  \
                                         storeFugacity,     \
                                         storeSaturation,   \
                                         storeDensity,      \
                                         storeViscosity,    \
                                         storeEnthalpy>

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 *
 * This class uses simpler and slightly less flexible template parameters as
 * ModularFluidState. Except for this, it is identical.
 */
template <class ScalarT,
          unsigned numPhasesV,
          unsigned numComponentsV,
          class FluidSystem, // only needed if the compositional stuff enabled
          bool storePressure,
          bool storeTemperature,
          bool storeComposition,
          bool storeFugacity,
          bool storeSaturation,
          bool storeDensity,
          bool storeViscosity,
          bool storeEnthalpy>
class SimpleModularFluidState
    : public ModularFluidState<ScalarT, numPhasesV, numComponentsV,
                               typename std::conditional<storePressure,
                                                         FluidStateExplicitPressureModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullPressureModule<ScalarT> >::type,
                               typename std::conditional<storeTemperature,
                                                         FluidStateExplicitTemperatureModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullTemperatureModule<ScalarT> >::type,
                               typename std::conditional<storeComposition,
                                                         FluidStateExplicitCompositionModule<ScalarT, FluidSystem, OPM_SMFS>,
                                                         FluidStateNullCompositionModule<ScalarT> >::type,
                               typename std::conditional<storeFugacity,
                                                         FluidStateExplicitFugacityModule<ScalarT, numPhasesV, numComponentsV, OPM_SMFS>,
                                                         FluidStateNullFugacityModule<ScalarT> >::type,
                               typename std::conditional<storeSaturation,
                                                         FluidStateExplicitSaturationModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullSaturationModule<ScalarT> >::type,
                               typename std::conditional<storeDensity,
                                                         FluidStateExplicitDensityModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullDensityModule<ScalarT, numPhasesV, OPM_SMFS> >::type,
                               typename std::conditional<storeViscosity,
                                                         FluidStateExplicitViscosityModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullViscosityModule<ScalarT, numPhasesV, OPM_SMFS> >::type,
                               typename std::conditional<storeEnthalpy,
                                                         FluidStateExplicitEnthalpyModule<ScalarT, numPhasesV, OPM_SMFS>,
                                                         FluidStateNullEnthalpyModule<ScalarT, numPhasesV, OPM_SMFS> >::type
                               >
{};

// we don't need the death-prevention macro anymore
#undef OPM_SMFS

} // namespace Opm

#endif
