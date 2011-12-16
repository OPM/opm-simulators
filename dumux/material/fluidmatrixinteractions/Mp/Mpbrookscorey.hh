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
 * \file
 * Implements a Brooks-Corey saturation-capillary pressure relation
 * for M-phase fluid systems.
 */
#ifndef DUMUX_MP_BROOKSCOREY_MATERIAL_HH
#define DUMUX_MP_BROOKSCOREY_MATERIAL_HH

#include "Mpbrookscoreymaterialparams.hh"

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
template <int numPhasesV, class ScalarT,
          class ParamsT = MpBrookscoreyMaterialParams<numPhasesV, ScalarT> >
class MpBrookscoreyMaterial
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = numPhasesV };

    /*!
     * \brief The Brooks-Corey capillary pressure-saturation curve.
     *
     * The Brooks-Corey material law is given by the relation:
     * \f[
     p_i - p_{i - 1} = p_{e,i}\overline{S}_{i}^{-1/\alpha}
     \f]
    */
    template <class pcContainerT, class SatContainerT>
    static void pC(pcContainerT &pc,
                   const Params &params,
                   const SatContainerT &saturations,
                   Scalar temperature)
    {
        for (int i = 0; i < numPhases; ++i) {
            Scalar S = saturations[i];
            assert(0 <= S && S <= 1);
            pc[i] =
                params.entryPressure(i) *
                std::pow(S, -1.0/params.alpha(i));
        }
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param saturations The container of saturation values to be filled
     * \param params Parameters
     * \param pc Capillary pressures
     * \param temperature Temperature
     */
    template <class SatContainerT, class pcContainerT>
    static void S(SatContainerT &saturations,
                  const Params &params,
                  const pcContainerT &pc,
                  Scalar temperature)
    {
        int refPhaseIdx = -1;
        for (int i = 0; i < numPhases; ++i) {
            Scalar p_Ci = pc[i];
            assert(0 =< p_Ci);
            if (params.entryPressure(i) == 0) {
                assert(refPhaseIdx == -1);
                refPhaseIdx = i;
                continue;
            }

            Scalar tmp = pow(p_Ci/params.entryPressure(i), -params.alpha(i));
            saturations[i] = tmp;
        }
    }

#warning TODO
#if 0
    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     \f]
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        return - (params.maxPC() - params.entryPC());
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return - 1/(params.maxPC() - params.entryPC());
    }
#endif

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class krContainerT, class SatContainerT>
    static void kr(krContainerT &kr,
                   const Params &params,
                   const SatContainerT &saturations,
                   Scalar temperature)
    {
        for (int i = 0; i < numPhases; ++i)
            // TODO: probably incorrect!
            kr[i] = pow(saturations[i], 2.0/params.alpha(i) + 3);
    };
};
}

#endif
