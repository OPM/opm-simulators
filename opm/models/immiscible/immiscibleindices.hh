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
 * \copydoc Opm::ImmiscibleIndices
 */
#ifndef EWOMS_IMMISCIBLE_INDICES_HH
#define EWOMS_IMMISCIBLE_INDICES_HH

#include "immiscibleproperties.hh"
#include <opm/models/common/energymodule.hh>

namespace Opm {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief The indices for the isothermal multi-phase model.
 */
template <class TypeTag, int PVOffset>
struct ImmiscibleIndices
    : public EnergyIndices<PVOffset + getPropValue<TypeTag, Properties::NumPhases>(),
                           getPropValue<TypeTag, Properties::EnableEnergy>()>
{
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    using EnergyIndices = Opm::EnergyIndices<PVOffset + numPhases, enableEnergy>;

public:
    // number of equations/primary variables
    static const int numEq = numPhases + EnergyIndices::numEq_;

    // Primary variable indices

    //! Index for wetting/non-wetting phase pressure
    //! (depending on formulation) in a solution vector
    static const int pressure0Idx = PVOffset + 0;
    //! Index of the saturation of the non-wetting/wetting phase
    static const int saturation0Idx = PVOffset + 1;

    // indices of the equations

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
};
} // namespace Opm

#endif
