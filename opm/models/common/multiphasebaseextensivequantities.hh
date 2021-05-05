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
 * \copydoc Opm::MultiPhaseBaseExtensiveQuantities
 */
#ifndef EWOMS_MULTI_PHASE_BASE_EXTENSIVE_QUANTITIES_HH
#define EWOMS_MULTI_PHASE_BASE_EXTENSIVE_QUANTITIES_HH

#include "multiphasebaseproperties.hh"

#include <opm/models/common/quantitycallbacks.hh>
#include <opm/models/discretization/common/fvbaseextensivequantities.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Opm {
/*!
 * \ingroup Discretization
 *
 * \brief This class calculates the pressure potential gradients and
 *        the filter velocities for multi-phase flow in porous media
 */
template <class TypeTag>
class MultiPhaseBaseExtensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscExtensiveQuantities>
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxExtensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscExtensiveQuantities>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using FluxExtensiveQuantities = typename FluxModule::FluxExtensiveQuantities;

public:
    /*!
     * \brief Register all run-time parameters for the extensive quantities.
     */
    static void registerParameters()
    {
        FluxModule::registerParameters();
    }

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
        ParentType::update(elemCtx, scvfIdx, timeIdx);

        // compute the pressure potential gradients
        FluxExtensiveQuantities::calculateGradients_(elemCtx, scvfIdx, timeIdx);

        // Check whether the pressure potential is in the same direction as the face
        // normal or in the opposite one
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(upstreamScvIdx_[phaseIdx]);
                Valgrind::SetUndefined(downstreamScvIdx_[phaseIdx]);
                continue;
            }

            upstreamScvIdx_[phaseIdx] = FluxExtensiveQuantities::upstreamIndex_(phaseIdx);
            downstreamScvIdx_[phaseIdx] = FluxExtensiveQuantities::downstreamIndex_(phaseIdx);
        }

        FluxExtensiveQuantities::calculateFluxes_(elemCtx, scvfIdx, timeIdx);
    }


    /*!
     * \brief Update the extensive quantities for a given boundary face.
     *
     * \param context Reference to the current execution context.
     * \param bfIdx The local index of the boundary face for which
     *              the extensive quantities should be calculated.
     * \param timeIdx The index used by the time discretization.
     * \param fluidState The FluidState on the domain boundary.
     * \param paramCache The FluidSystem's parameter cache.
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context& context,
                        unsigned bfIdx,
                        unsigned timeIdx,
                        const FluidState& fluidState)
    {
        ParentType::updateBoundary(context, bfIdx, timeIdx, fluidState);

        FluxExtensiveQuantities::calculateBoundaryGradients_(context.elementContext(),
                                                             bfIdx,
                                                             timeIdx,
                                                             fluidState);
        FluxExtensiveQuantities::calculateBoundaryFluxes_(context.elementContext(),
                                                          bfIdx,
                                                          timeIdx);
    }

    /*!
     * \brief Return the local index of the upstream control volume for a given phase as
     *        a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIndex(unsigned phaseIdx) const
    { return upstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIndex(unsigned phaseIdx) const
    { return downstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the weight of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar upstreamWeight(unsigned phaseIdx OPM_UNUSED) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(unsigned phaseIdx) const
    { return 1.0 - upstreamWeight(phaseIdx); }

private:
    short upstreamScvIdx_[numPhases];
    short downstreamScvIdx_[numPhases];
};

} // namespace Opm

#endif
