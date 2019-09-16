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
#ifndef EWOMS_FLASH_INDICES_HH
#define EWOMS_FLASH_INDICES_HH

#include "flashproperties.hh"
#include <opm/models/common/energymodule.hh>

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief Defines the primary variable and equation indices for the
 *        compositional multi-phase model based on flash calculations.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class FlashIndices
    : public EnergyIndices<PVOffset + GET_PROP_VALUE(TypeTag, NumComponents),
                           GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Opm::EnergyIndices<PVOffset + numComponents, enableEnergy> EnergyIndices;

public:
    //! number of equations/primary variables
    static const int numEq = numComponents + EnergyIndices::numEq_;

    // Primary variable indices

    //! Index of the total concentration of the first component in the
    //! pore space.
    static const int cTot0Idx = PVOffset;

    // equation indices

    //! Index of the mass conservation equation for the first
    //! component.
    static const int conti0EqIdx = PVOffset;
};

} // namespace Opm

#endif
