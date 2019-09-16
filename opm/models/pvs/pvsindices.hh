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
 * \copydoc Opm::PvsIndices
 */
#ifndef EWOMS_PVS_INDICES_HH
#define EWOMS_PVS_INDICES_HH

#include "pvsproperties.hh"

#include <opm/models/common/energymodule.hh>

namespace Opm {
/*!
 * \ingroup PvsModel
 *
 * \brief The indices for the compositional multi-phase primary
 *        variable switching model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class PvsIndices
    : public EnergyIndices<PVOffset + GET_PROP_VALUE(TypeTag, NumComponents),
                           GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Opm::EnergyIndices<PVOffset + numComponents, enableEnergy> EnergyIndices;

public:
    //! Number of partial differential equations or primary variables, respectively
    static const int numEq = numComponents + EnergyIndices::numEq_;

    // Primary variable indices

    //! Index for the pressure of the first phase
    static const int pressure0Idx = PVOffset + 0;
    //! Index of the either the saturation or the mole
    //! fraction of the phase with the lowest index
    static const int switch0Idx = PVOffset + 1;

    // equation indices

    //! Index of the mass conservation equation for the first component
    static const int conti0EqIdx = PVOffset;
};

} // namespace Opm

#endif
