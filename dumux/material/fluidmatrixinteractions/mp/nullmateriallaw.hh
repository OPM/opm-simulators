// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \file nullmateriallaw.hh
 *
 * Implements a dummy linear saturation-capillary pressure relation
 * which just disables capillary pressure.
 */
#ifndef DUMUX_NULL_MATERIAL_LAW_HH
#define DUMUX_NULL_MATERIAL_LAW_HH

#include "nullmateriallawparams.hh"

#include <dune/common/exceptions.hh>

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements a dummy linear saturation-capillary pressure
 *        relation which just disables capillary pressure.
 */
template <int numPhasesV, class ScalarT, class ParamsT = NullMaterialLawParams<numPhasesV, ScalarT> >
class NullMaterialLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = numPhasesV };

    /*!
     * \brief Returns constant 0 for all phases.
     *
     * \param values Container for the return values
     * \param params Parameters
     * \param state The fluid state
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            values[phaseIdx] = 0.0;
    }

    /*!
     * \brief The inverse of the capillary pressure
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &state)
    {
        DUNE_THROW(Dune::NotImplemented, "MpLinearMaterial::saturations()");
    }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            values[phaseIdx] = std::max(std::min(state.saturation(phaseIdx),1.0),0.0);
        }
    }
};
}

#endif
