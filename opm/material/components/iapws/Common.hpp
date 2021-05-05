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
 * \copydoc Opm::IAPWS::Common
 */
#ifndef OPM_IAPWS_COMMON_HPP
#define OPM_IAPWS_COMMON_HPP

#include <opm/material/Constants.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
namespace IAPWS {

/*!
 *
 *  \ingroup IAPWS
 *
 * \brief Implements relations which are common for all regions of the IAPWS '97
 *        formulation.
 *
 * \tparam Scalar The type used for scalar values
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class Common
{
public:
    //! The molar mass of water \f$\mathrm{[kg/mol]}\f$
    static const Scalar molarMass;

    //! Specific gas constant of water \f$\mathrm{[J/(kg*K)]}\f$
    static const Scalar Rs;

    //! Critical temperature of water \f$\mathrm{[K]}\f$
    static const Scalar criticalTemperature;

    //! Critical pressure of water \f$\mathrm{[Pa]}\f$
    static const Scalar criticalPressure;

    //! Density of water at the critical point \f$\mathrm{[kg/m^3]}\f$
    static const Scalar criticalDensity;

    //! Critical molar volume of water \f$\mathrm{[m^3/mol]}\f$
    static const Scalar criticalMolarVolume;

    //! The acentric factor of water \f$\mathrm{[-]}\f$
    static const Scalar acentricFactor;

    //! Triple temperature of water \f$\mathrm{[K]}\f$
    static const Scalar tripleTemperature;

    //! Triple pressure of water \f$\mathrm{[Pa]}\f$
    static const Scalar triplePressure;

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[(N/m^2)*s]}\f$of pure water.
     *
     * This relation is valid for all regions of the IAPWS '97
     * formulation.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param rho density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    template <class Evaluation>
    static Evaluation viscosity(const Evaluation& temperature, const Evaluation& rho)
    {
        Evaluation rhoBar = rho/322.0;
        Evaluation TBar = temperature/criticalTemperature;

        // muBar = muBar_1
        const Scalar Hij[6][7] = {
            { 5.20094e-1, 2.22531e-1,-2.81378e-1, 1.61913e-1,-3.25372e-2, 0, 0 },
            { 8.50895e-2, 9.99115e-1,-9.06851e-1, 2.57399e-1, 0, 0, 0 },
            { -1.08374, 1.88797 ,-7.72479e-1, 0, 0, 0, 0 },
            { -2.89555e-1, 1.26613 ,-4.89837e-1, 0, 6.98452e-2, 0,-4.35673e-3 },
            { 0, 0,-2.57040e-1, 0, 0, 8.72102e-3, 0 },
            { 0, 1.20573e-1, 0, 0, 0, 0,-5.93264e-4 }
        };

        Evaluation tmp, tmp2, tmp3 = 1;
        Evaluation muBar = 0;
        for (int i = 0; i <= 5; ++i) {
            tmp = 0;
            tmp2 = 1;
            for (int j = 0; j <= 6; ++j) {
                tmp += Hij[i][j]*tmp2;
                tmp2 *= (rhoBar - 1);
            };
            muBar += tmp3 * tmp;
            tmp3 *= 1.0/TBar - 1;
        };
        muBar *= rhoBar;
        muBar = exp(muBar);

        // muBar *= muBar_0
        muBar  *= 100*sqrt(TBar);
        const Scalar H[4] = {
            1.67752, 2.20462, 0.6366564, -0.241605
        };

        tmp = 0, tmp2 = 1;
        for (int i = 0; i < 4; ++i) {
            tmp += H[i]/tmp2;
            tmp2 *= TBar;
        };
        muBar /= tmp;

        return 1e-6*muBar;
    }

    /*!
    * \brief Thermal conductivity \f$\mathrm{[[W/(m^2 K/m)]}\f$ water (IAPWS) .
    *
    * Implementation taken from:
    * freesteam - IAPWS-IF97 steam tables library
    * Copyright (C) 2004-2009  John Pye
    *
    * Appendix B: Recommended Interpolating equation for Industrial Use
    * see http://www.iapws.org/relguide/thcond.pdf
    *
    * \param T absolute temperature in K
    * \param rho density of water in kg/m^3
    */
    template <class Evaluation>
    static Evaluation thermalConductivityIAPWS(const Evaluation& T, const Evaluation& rho)
    {
        static const Scalar thcond_tstar = 647.26 ;
        static const Scalar thcond_rhostar = 317.7 ;
        /*static const Scalar thcond_kstar = 1.0 ;*/

        static const Scalar thcond_b0 = -0.397070 ;
        static const Scalar thcond_b1 = 0.400302 ;
        static const Scalar thcond_b2 = 1.060000 ;
        static const Scalar thcond_B1 = -0.171587 ;
        static const Scalar thcond_B2 = 2.392190 ;

        static const Scalar thcond_c1 = 0.642857 ;
        static const Scalar thcond_c2 = -4.11717 ;
        static const Scalar thcond_c3 = -6.17937 ;
        static const Scalar thcond_c4 = 0.00308976 ;
        static const Scalar thcond_c5 = 0.0822994 ;
        static const Scalar thcond_c6 = 10.0932 ;

        static const Scalar thcond_d1 = 0.0701309 ;
        static const Scalar thcond_d2 = 0.0118520 ;
        static const Scalar thcond_d3 = 0.00169937 ;
        static const Scalar thcond_d4 = -1.0200 ;
        static const int thcond_a_count = 4;
        static const Scalar thcond_a[thcond_a_count] = {
            0.0102811
            ,0.0299621
            ,0.0156146
            ,-0.00422464
        };

        Evaluation Tbar = T / thcond_tstar;
        Evaluation rhobar = rho / thcond_rhostar;

        /* fast implementation... minimised calls to 'pow' routine... */
        Evaluation Troot = sqrt(Tbar);
        Evaluation Tpow = Troot;
        Evaluation lam = 0;

        for(int k = 0; k < thcond_a_count; ++k) {
            lam += thcond_a[k] * Tpow;
            Tpow *= Tbar;
        }

        lam +=
            thcond_b0 + thcond_b1
            * rhobar + thcond_b2
            * exp(thcond_B1 * ((rhobar + thcond_B2)*(rhobar + thcond_B2)));

        Evaluation DTbar = abs(Tbar - 1) + thcond_c4;
        Evaluation DTbarpow = pow(DTbar, 3./5);
        Evaluation Q = 2. + thcond_c5 / DTbarpow;

        Evaluation S;
        if(Tbar >= 1)
            S = 1. / DTbar;
        else
            S = thcond_c6 / DTbarpow;

        Evaluation rhobar18 = pow(rhobar, 1.8);
        Evaluation rhobarQ = pow(rhobar, Q);

        lam +=
            (thcond_d1 / pow(Tbar,10.0) + thcond_d2) * rhobar18 *
            exp(thcond_c1 * (1 - rhobar * rhobar18))
            + thcond_d3 * S * rhobarQ *
            exp((Q/(1+Q))*(1 - rhobar*rhobarQ))
            + thcond_d4 *
            exp(thcond_c2 * pow(Troot,3.0) + thcond_c3 / pow(rhobar,5.0));
        return /*thcond_kstar * */ lam;
    }
};

template <class Scalar>
const Scalar Common<Scalar>::molarMass = 18.01518e-3;
template <class Scalar>
const Scalar Common<Scalar>::Rs = Constants<Scalar>::R/molarMass;
template <class Scalar>
const Scalar Common<Scalar>::criticalTemperature = 647.096;
template <class Scalar>
const Scalar Common<Scalar>::criticalPressure = 22.064e6;
template <class Scalar>
const Scalar Common<Scalar>::criticalDensity = 322.0;
template <class Scalar>
const Scalar Common<Scalar>::criticalMolarVolume = molarMass/criticalDensity;
template <class Scalar>
const Scalar Common<Scalar>::acentricFactor = 0.344;
template <class Scalar>
const Scalar Common<Scalar>::tripleTemperature = 273.16;
template <class Scalar>
const Scalar Common<Scalar>::triplePressure = 611.657;

} // namespace IAPWS
} // namespace Opm

#endif
