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
 * \copydoc Opm::BlackOilExtensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_EXTENSIVE_QUANTITIES_HH
#define EWOMS_BLACK_OIL_EXTENSIVE_QUANTITIES_HH

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilenergymodules.hh"

#include <opm/models/common/multiphasebaseextensivequantities.hh>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 * \ingroup ExtensiveQuantities
 *
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the black-oil model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class BlackOilExtensiveQuantities
    : public MultiPhaseBaseExtensiveQuantities<TypeTag>
    , public BlackOilSolventExtensiveQuantities<TypeTag>
    , public BlackOilPolymerExtensiveQuantities<TypeTag>
    , public BlackOilEnergyExtensiveQuantities<TypeTag>
{
    typedef MultiPhaseBaseExtensiveQuantities<TypeTag> MultiPhaseParent;

    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief Update the extensive quantities for a given sub-control-volume-face.
     *
     * \param elemCtx Reference to the current element context.
     * \param scvfIdx The local index of the sub-control-volume face for
     *                which the extensive quantities should be calculated.
     * \param timeIdx The index used by the time discretization.
     */
    void update(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        MultiPhaseParent::update(elemCtx, scvfIdx, timeIdx);

        asImp_().updateSolvent(elemCtx, scvfIdx, timeIdx);
        asImp_().updatePolymer(elemCtx, scvfIdx, timeIdx);
        asImp_().updateEnergy(elemCtx, scvfIdx, timeIdx);
    }

    template <class Context, class FluidState>
    void updateBoundary(const Context& ctx,
                        unsigned bfIdx,
                        unsigned timeIdx,
                        const FluidState& fluidState)
    {
        MultiPhaseParent::updateBoundary(ctx, bfIdx, timeIdx, fluidState);

        asImp_().updateEnergyBoundary(ctx, bfIdx, timeIdx, fluidState);
    }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Opm

#endif
