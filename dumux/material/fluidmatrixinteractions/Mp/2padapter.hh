/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \file 2padapter.hh
 *
 * Makes the twophase capillary pressure-saturation relations
 * available under the M-phase API for material laws
 */
#ifndef DUMUX_MP_2P_ADAPTER_HH
#define DUMUX_MP_2P_ADAPTER_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements a brookscorey saturation-capillary pressure relation
 *
 * Implements a brookscorey saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 * \sa MpBrookscoreyMaterialParams
 */
template <int wPhaseIdx, class TwoPLaw >
class TwoPAdapter
{
    enum { nPhaseIdx = (wPhaseIdx == 0)?1:0 };

public:
    typedef typename TwoPLaw::Params Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = 2 };

    /*!
     * \brief The capillary pressure-saturation curve.
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params, 
                                   const FluidState &state)
    {
        // non-wetting phase gets the capillary pressure added
        values[nPhaseIdx] = 0;

        // wetting phase does not get anything added
        values[wPhaseIdx] = - TwoPLaw::pC(params, state.saturation(wPhaseIdx)); 
    }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params, 
                                       const FluidState &state)
    {
        values[wPhaseIdx] = TwoPLaw::krw(params, state.saturation(wPhaseIdx));
        values[nPhaseIdx] = TwoPLaw::krn(params, state.saturation(wPhaseIdx));
    };
};
}

#endif
