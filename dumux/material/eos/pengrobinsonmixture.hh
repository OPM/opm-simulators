// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \file
 *
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
#ifndef DUMUX_PENG_ROBINSON_MIXTURE_HH
#define DUMUX_PENG_ROBINSON_MIXTURE_HH

#include <stdlib.h>
#include "pengrobinson.hh"

#include <dumux/material/constants.hh>

namespace Dumux
{

/*!
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
template <class Scalar, class StaticParameters>
class PengRobinsonMixture
{
    enum { numComponents = StaticParameters::numComponents };
    typedef Dumux::PengRobinson<Scalar> PengRobinson;

    // this class cannot be instantiated!
    PengRobinsonMixture() {};

    // the ideal gas constant
    static constexpr Scalar R = Dumux::Constants<Scalar>::R;

    // the u and w parameters as given by the Peng-Robinson EOS
    static constexpr Scalar u = 2.0;
    static constexpr Scalar w = -1.0;

public:
    /*!
     * \brief Computes molar volumes where the Peng-Robinson EOS is
     *        true.
     *
     * \return Number of solutions.
     */
    template <class MutableParams, class FluidState>
    static int computeMolarVolumes(Scalar *Vm,
                                   const MutableParams &params,
                                   int phaseIdx,
                                   const FluidState &fs)
    {
        return PengRobinson::computeMolarVolumes(Vm, params, phaseIdx, fs);
    }

    /*!
     * \brief Returns the fugacity coefficient of an individual
     *        component in the phase.
     *
     * The fugacity coefficient \f$\phi_i\f$ of a component $i$ is
     * defined as
     * \f[
     f_i = \phi_i x_i \;,
     \f]
     * where \f$f_i\f$ is the component's fugacity and \f$x_i\f$ is
     * the component's mole fraction.
     *
     * See:
     *
      * R. Reid, et al.: The Properties of Gases and Liquids,
      * 4th edition, McGraw-Hill, 1987, pp. 42-44, 143-145
      */
    template <class Params>
    static Scalar computeFugacityCoeff(const Params &params,
                                       int phaseIdx,
                                       int compIdx)
    {
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.

        Scalar Vm = params.molarVolume(phaseIdx);

        // Calculate b_i / b
        Scalar bi_b = params.bPure(phaseIdx, compIdx) / params.b(phaseIdx);

        // Calculate the compressibility factor
        Scalar RT = R*params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx); // molar volume in [bar]
        Scalar Z = p*Vm/RT; // compressibility factor

        // Calculate A^* and B^* (see: Reid, p. 42)
        Scalar Astar = params.a(phaseIdx)*p/(RT*RT);
        Scalar Bstar = params.b(phaseIdx)*p/(RT);


        // calculate delta_i (see: Reid, p. 145)
        Scalar deltai = 2*std::sqrt(params.aPure(phaseIdx, compIdx))/params.a(phaseIdx);
        Scalar tmp = 0;
        for (int j = 0; j < numComponents; ++j) {
            tmp +=
                params.moleFrac(phaseIdx, j) / params.sumMoleFrac(phaseIdx) *
                std::sqrt(params.aPure(phaseIdx, j))
                *(1.0 - StaticParameters::interactionCoefficient(compIdx, j));
        };
        deltai *= tmp;

        Scalar base =
            (2*Z + Bstar*(u + std::sqrt(u*u - 4*w))) /
            (2*Z + Bstar*(u - std::sqrt(u*u - 4*w)));
        Scalar expo =  Astar/(Bstar*std::sqrt(u*u - 4*w))*(bi_b - deltai);

        Scalar fugCoeff =
            std::exp(bi_b*(Z - 1))/std::max(1e-9, Z - Bstar) *
            std::pow(base, expo);

        ////////
        // limit the fugacity coefficient to a reasonable range:
        //
        // on one side, we want the mole fraction to be at
        // least 10^-3 if the fugacity is at the current pressure
        //
        fugCoeff = std::min(1e10, fugCoeff);
        //
        // on the other hand, if the mole fraction of the component is 100%, we want the
        // fugacity to be at least 10^-3 Pa
        //
        fugCoeff = std::max(1e-10, fugCoeff);
        ///////////
        if (!std::isfinite(fugCoeff)) {
            std::cout << "Non finite phi: " << fugCoeff << "\n";
        }

        return fugCoeff;
    };

};
} // end namepace Dumux

#endif
