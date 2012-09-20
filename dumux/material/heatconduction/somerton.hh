// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \file somerton.hh
 *
 * \brief Implements the Somerton law of heat conductivity in a
 * porous medium.
 *
 * See:
 *
 * W.H. Somerton, A.H. El-Shaarani and S.M. Mobarak: High
 * Temperature Behavior of Rocks Associated with Geothermal Type
 * Reservoirs, paper SPE-4897 presentet at SPE California Regional
 * Meeting 1974, 1974
 *
 * or
 *
 * H. Class: Theorie und numerische Modellierung nichtisothermer
 * Mehrphasenprozesse in NAPL kontaminierten poroesen Medien, PhD
 * thesis, Technical University of Braunschweig, 2000
 */
#ifndef DUMUX_SOMERTON_HH
#define DUMUX_SOMERTON_HH

#include "somertonparams.hh"

#include <dumux/common/spline.hh>
#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements the Somerton law of heat conductivity in a
 * porous medium.
 *
 * See:
 *
 * W.H. Somerton, A.H. El-Shaarani and S.M. Mobarak: High
 * Temperature Behavior of Rocks Associated with Geothermal Type
 * Reservoirs, paper SPE-4897 presentet at SPE California Regional
 * Meeting 1974, 1974
 *
 * or
 *
 * H. Class: Theorie und numerische Modellierung nichtisothermer
 * Mehrphasenprozesse in NAPL kontaminierten poroesen Medien, PhD
 * thesis, Technical University of Braunschweig, 2000
 */
template <class FluidSystem,
          class ScalarT,
          class ParamsT = SomertonParams<FluidSystem::numPhases, ScalarT> >
class Somerton
{
    enum { numPhases = FluidSystem::numPhases };

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective heat conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     *
     * For two phases, the Somerton law is given by:
     * \f[
     \lambda_{pm} =
     \lambda_{ful,g} +
     \sqrt{S_w}(\lambda_{ful,w} - \lambda_{vac}) +
     \sqrt{S_n}(\lambda_{ful,n} - \lambda_{vac})
     \f]
     *
     * where \f$\lambda_{vac}\f$ is the heat conductivity of the
     * porous medium at vacuum, \f$\lambda_{ful,\alpha}\f$ is the heat
     * conductivty of the porous medium if it is fully saturated by
     * phase \f$\alpha\f$ and \f$S_\alpha\f$ is the saturation of
     * phase \f$\alpha\f$.
     */
    template <class FluidState>
    static Scalar heatConductivity(const Params &params,
                                   const FluidState &fluidState)
    {
        Valgrind::CheckDefined(params.vacuumLambda());

        Scalar lambda = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Valgrind::CheckDefined(params.fullySaturatedLambda(phaseIdx));

            if (FluidSystem::isLiquid(phaseIdx)) {
                lambda +=
                    regularizedSqrt_(fluidState.saturation(phaseIdx))
                    * (params.fullySaturatedLambda(phaseIdx) - params.vacuumLambda());
            }
            else { // gas phase
                lambda += params.fullySaturatedLambda(phaseIdx) - params.vacuumLambda();
            }
        };

        lambda += params.vacuumLambda();
        return lambda;
    }

protected:
    static Scalar regularizedSqrt_(Scalar x)
    {
        const Scalar xMin = 1e-2;
        const Scalar sqrtXMin = std::sqrt(xMin);
        const Scalar fPrimeXMin = 1.0/(2*std::sqrt(xMin));
        const Scalar fPrime0 = 2*fPrimeXMin;
        typedef Dumux::Spline<Scalar, 2> Spline;
        static Spline sqrtRegSpline(0, xMin, // x0, x1
                                    0, sqrtXMin, // y0, y1
                                    fPrime0, fPrimeXMin); // m0, m1

        if (x > xMin)
            return std::sqrt(x);
        else if (x <= 0)
            return fPrime0 * x;
        else
            return sqrtRegSpline.eval(x);
    }
};
}

#endif
