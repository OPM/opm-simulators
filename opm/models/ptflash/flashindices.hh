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
 *
 * \copydoc Opm::FlashIndices
 */
#ifndef OPM_PTFLASH_INDICES_HH
#define OPM_PTFLASH_INDICES_HH

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/multiphasebaseproperties.hh>

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief Defines the primary variable and equation indices for the
 *        compositional multi-phase model based on PT flash calculations.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
// TODO: The indices class should handle whether phase is active, not the FluidSystem
template <class TypeTag, int PVOffset>
class FlashIndices
    : public EnergyIndices<PVOffset + getPropValue<TypeTag, Properties::NumComponents>(),
                           getPropValue<TypeTag, Properties::EnableEnergy>()>
{
    static constexpr int numComponents = getPropValue<TypeTag, Properties::NumComponents>();
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableWater = getPropValue<TypeTag, Properties::EnableWater>() };
    using EnergyIndices = ::Opm::EnergyIndices<PVOffset + numComponents, enableEnergy>;

public:
    //! All phases active (note: immiscible/"dummy" water phase)
    static constexpr bool waterEnabled = enableWater;
    static constexpr bool gasEnabled = true;
    static constexpr bool oilEnabled = true;

    //! number of active phases
    static constexpr int numPhases = enableWater ? 3 : 2;

    //! number of equations/primary variables
    static constexpr int numEq = numComponents + EnergyIndices::numEq_ + (enableWater ? 1 : 0);

    // Primary variable indices

    //! Index of the pressure
    static constexpr int pressure0Idx = PVOffset;

    //! Index of the molefraction of the first component
    static constexpr int z0Idx = pressure0Idx + 1;

    //! Index of water saturation
    static constexpr int water0Idx = enableWater ? z0Idx + numComponents - 1 : -1000;

    // equation indices

    //! Index of the mass conservation equation for the first
    //! component.
    static constexpr int conti0EqIdx = PVOffset;

    //! MICP implementation is an extension of the black-oil model
    static constexpr bool enableMICP = false;
};

} // namespace Opm

#endif
