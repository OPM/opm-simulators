// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \file Mplinearmaterial.hh
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 */
#ifndef MP_LINEAR_MATERIAL_HH
#define MP_LINEAR_MATERIAL_HH

#include "Mplinearmaterialparams.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 * \sa MpLinearMaterialParams
 */
template <int numPhasesV, class ScalarT, class ParamsT = MpLinearMaterialParams<numPhasesV, ScalarT> >
class MpLinearMaterial
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = numPhasesV };

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
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
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar S = state.saturation(phaseIdx);
            values[phaseIdx] =
                S*params.pcMaxSat(phaseIdx) +
                (1 - S)*params.pcMinSat(phaseIdx);
        }
    }

#if 0
#warning TODO
    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    template <class SatContainerT, class FluidState>
    static void saturations(SatContainerT &saturations,
                            const Params &params,
                            const FluidState &state)
    {
        for (int i = 0; i < numPhases; ++i) {
            saturations[i] =
                (pc[i] - params.pcMaxSat(i))
                /
                (params.pcMinSat(i) - params.pcMaxSat(i));
        }
    }

#endif

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            values[phaseIdx] = std::max(std::min(state.saturation(phaseIdx),1.0),0.0);
    };
};
}

#endif
