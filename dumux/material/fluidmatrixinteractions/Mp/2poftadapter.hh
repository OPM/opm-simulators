/*****************************************************************************
 *   Copyright (C) 2010 by Philipp Nuske                                     *
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
 * \file 2poftadapter.hh
 *
 * Makes the twophase capillary pressure-saturation relations
 * available under the M-phase API for material laws.
 *
 * Also use the temperature dependent version of the material laws.
 */
#ifndef DUMUX_MP_2P_OFT_ADAPTER_HH
#define DUMUX_MP_2P_OFT_ADAPTER_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Adapts the interface of the MpNc material law to the standard-\Dumux material law.
 *
 *        Also use the temperature dependent version of the material laws.
 */
template <int wPhaseIdx, class TwoPLaw>
class TwoPOfTAdapter
{
    enum { nPhaseIdx = (wPhaseIdx == 0)?1:0 };

public:
    typedef typename TwoPLaw::Params Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = 2 };

    /*!
     * \brief The capillary pressure-saturation curve.
     */
    template <class pcContainerT, class FluidState>
    static void capillaryPressures(pcContainerT &pc,
                   const Params &params, 
                   const FluidState &fluidState)
    {
        // non-wetting phase gets the capillary pressure added
        pc[nPhaseIdx] = 0;

        // wetting phase does not get anything added
        pc[wPhaseIdx] = - TwoPLaw::pC(params, fluidState.saturation(wPhaseIdx), fluidState.temperature(wPhaseIdx));
    }



    /*!
     * \brief The relative permeability of all phases.
     */
    template <class krContainerT, class FluidState>
    static void relativePermeabilities(krContainerT &kr,
                   const Params &params, 
                   const FluidState &fluidState)
    {
        kr[wPhaseIdx] = TwoPLaw::krw(params, fluidState.saturation(wPhaseIdx));
        kr[nPhaseIdx] = TwoPLaw::krn(params, fluidState.saturation(wPhaseIdx));
    };
};
}

#endif
