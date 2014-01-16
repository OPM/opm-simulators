/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2012 by Bernd Flemisch

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
*/
/*!
 * \file
 * \copydoc Opm::ThreePAdapter
 */
#ifndef OPM_3P_ADAPTER_HPP
#define OPM_3P_ADAPTER_HPP

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/material/Valgrind.hpp>

#include <algorithm>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Makes the three-phase capillary pressure-saturation relations
 *        available under the M-phase API for material laws.
 */
template <int wPhaseIdx, int nPhaseIdx, int gPhaseIdx, class ThreePLaw>
class ThreePAdapter
{
public:
    typedef typename ThreePLaw::Params Params;
    typedef typename ThreePLaw::Scalar Scalar;
    enum { numPhases = 3 };

    /*!
     * \brief The capillary pressure-saturation curve.
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &fluidState)
    {
        Scalar p_cgw = ThreePLaw::pCGW(params, fluidState.saturation(wPhaseIdx));
        Scalar p_cnw = ThreePLaw::pCNW(params, fluidState.saturation(wPhaseIdx));
        Scalar p_cgn = ThreePLaw::pCGN(params,
                                       fluidState.saturation(wPhaseIdx)
                                       + fluidState.saturation(nPhaseIdx));
        Scalar p_cAlpha = ThreePLaw::pCAlpha(params,
                                             fluidState.saturation(nPhaseIdx));
        Scalar p_cnw1 = 0.0;
        Valgrind::CheckDefined(p_cgw);
        Valgrind::CheckDefined(p_cnw);
        Valgrind::CheckDefined(p_cgn);
        Valgrind::CheckDefined(p_cAlpha);

        values[gPhaseIdx] = 0;
        values[nPhaseIdx] = - (p_cAlpha*p_cgn + (1 - p_cAlpha)*(p_cgw - p_cnw1));
        values[wPhaseIdx] = values[nPhaseIdx] - (p_cAlpha*p_cnw + (1 - p_cAlpha)*p_cnw1);
    }

    static Scalar pCGW(const Params &params, Scalar Sw)
    { return ThreePLaw::pCGW(params, Sw); }

    /*!
     * \brief The inverse of the capillary pressure-saturation curve.
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &fluidState)
    {  OPM_THROW(std::runtime_error, "Not implemented: Inverse capillary pressure curves"); }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &fluidState)
    {
        Scalar Sw = fluidState.saturation(wPhaseIdx);
        Scalar Sn = fluidState.saturation(nPhaseIdx);
        Scalar Sg = fluidState.saturation(gPhaseIdx);

        values[wPhaseIdx] = ThreePLaw::krw(params, Sw, Sn, Sg);
        values[nPhaseIdx] = ThreePLaw::krn(params, Sw, Sn, Sg);
        values[gPhaseIdx] = ThreePLaw::krg(params, Sw, Sn, Sg);
    }
};
} // namespace Opm

#endif
