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
 * \copydoc Opm::NcpExtensiveQuantities
 */
#ifndef EWOMS_NCP_EXTENSIVE_QUANTITIES_HH
#define EWOMS_NCP_EXTENSIVE_QUANTITIES_HH

#include "ncpproperties.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>
#include <opm/models/common/multiphasebaseextensivequantities.hh>

namespace Opm {

/*!
 * \ingroup NcpModel
 * \ingroup ExtensiveQuantities
 *
 * \brief This template class represents the extensive quantities of the compositional
 *        NCP model.
 */
template <class TypeTag>
class NcpExtensiveQuantities
    : public MultiPhaseBaseExtensiveQuantities<TypeTag>
    , public EnergyExtensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>
    , public DiffusionExtensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>()>
{
    using ParentType = MultiPhaseBaseExtensiveQuantities<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    using DiffusionExtensiveQuantities = Opm::DiffusionExtensiveQuantities<TypeTag, enableDiffusion>;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    using EnergyExtensiveQuantities = Opm::EnergyExtensiveQuantities<TypeTag, enableEnergy>;

public:
    /*!
     * \copydoc MultiPhaseBaseExtensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);
        DiffusionExtensiveQuantities::update_(elemCtx, scvfIdx, timeIdx);
        EnergyExtensiveQuantities::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc MultiPhaseBaseExtensiveQuantities::updateBoundary
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context& context,
                        unsigned bfIdx,
                        unsigned timeIdx,
                        const FluidState& fluidState)
    {
        ParentType::updateBoundary(context, bfIdx, timeIdx, fluidState);
        DiffusionExtensiveQuantities::updateBoundary_(context, bfIdx, timeIdx, fluidState);
        EnergyExtensiveQuantities::updateBoundary_(context, bfIdx, timeIdx, fluidState);
    }
};

} // namespace Opm

#endif
