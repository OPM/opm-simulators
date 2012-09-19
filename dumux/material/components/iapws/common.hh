// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
 *   Copyright (C) 2010 by Felix Bode                                        *
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
 * \brief Implements relations common for all regions of the IAPWS '97
 *        formulation.
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
#ifndef DUMUX_IAPWS_COMMON_HH
#define DUMUX_IAPWS_COMMON_HH

#include <dumux/material/constants.hh>

#include <cmath>

namespace Dumux {
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
    static constexpr Scalar molarMass = 18.01518e-3;

    //! Specific gas constant of water \f$\mathrm{[J/(kg*K)]}\f$
    static constexpr Scalar Rs = Dumux::Constants<Scalar>::R/molarMass;

    //! Critical temperature of water \f$\mathrm{[K]}\f$
    static constexpr Scalar criticalTemperature = 647.096;

    //! Critical pressure of water \f$\mathrm{[Pa]}\f$
    static constexpr Scalar criticalPressure = 22.064e6;

    //! Density of water at the critical point \f$\mathrm{[kg/m^3]}\f$
    static constexpr Scalar criticalDensity = 322.0;

    //! Critical molar volume of water \f$\mathrm{[m^3/mol]}\f$
    static constexpr Scalar criticalMolarVolume = molarMass/criticalDensity;

    //! The acentric factor of water \f$\mathrm{[-]}\f$
    static constexpr Scalar acentricFactor = 0.344;

    //! Triple temperature of water \f$\mathrm{[K]}\f$
    static constexpr Scalar tripleTemperature = 273.16;

    //! Triple pressure of water \f$\mathrm{[Pa]}\f$
    static constexpr Scalar triplePressure = 611.657;

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
    static Scalar viscosity(Scalar temperature, Scalar rho)
    {
        Scalar rhoBar = rho/322.0;
        Scalar TBar = temperature/criticalTemperature;

        // muBar = muBar_1
        const Scalar Hij[6][7] = {
            { 5.20094e-1, 2.22531e-1,-2.81378e-1, 1.61913e-1,-3.25372e-2, 0, 0 },
            { 8.50895e-2, 9.99115e-1,-9.06851e-1, 2.57399e-1, 0, 0, 0 },
            { -1.08374, 1.88797 ,-7.72479e-1, 0, 0, 0, 0 },
            { -2.89555e-1, 1.26613 ,-4.89837e-1, 0, 6.98452e-2, 0,-4.35673e-3 },
            { 0, 0,-2.57040e-1, 0, 0, 8.72102e-3, 0 },
            { 0, 1.20573e-1, 0, 0, 0, 0,-5.93264e-4 }
        };

        Scalar tmp, tmp2, tmp3 = 1;
        Scalar muBar = 0;
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
        muBar = std::exp(muBar);

        // muBar *= muBar_0
        muBar  *= 100*std::sqrt(TBar);
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
    static Scalar thermalConductivityIAPWS(Scalar T, Scalar rho)
    {
        static constexpr Scalar thcond_tstar = 647.26 ;
        static constexpr Scalar thcond_rhostar = 317.7 ;
        /*static constexpr Scalar thcond_kstar = 1.0 ;*/

        static constexpr Scalar thcond_b0 = -0.397070 ;
        static constexpr Scalar thcond_b1 = 0.400302 ;
        static constexpr Scalar thcond_b2 = 1.060000 ;
        static constexpr Scalar thcond_B1 = -0.171587 ;
        static constexpr Scalar thcond_B2 = 2.392190 ;

        static constexpr Scalar thcond_c1 = 0.642857 ;
        static constexpr Scalar thcond_c2 = -4.11717 ;
        static constexpr Scalar thcond_c3 = -6.17937 ;
        static constexpr Scalar thcond_c4 = 0.00308976 ;
        static constexpr Scalar thcond_c5 = 0.0822994 ;
        static constexpr Scalar thcond_c6 = 10.0932 ;

        static constexpr Scalar thcond_d1 = 0.0701309 ;
        static constexpr Scalar thcond_d2 = 0.0118520 ;
        static constexpr Scalar thcond_d3 = 0.00169937 ;
        static constexpr Scalar thcond_d4 = -1.0200 ;
        static constexpr int thcond_a_count = 4;
        static constexpr Scalar thcond_a[thcond_a_count] = {
            0.0102811
            ,0.0299621
            ,0.0156146
            ,-0.00422464
        };

        Scalar Tbar = T / thcond_tstar;
        Scalar rhobar = rho / thcond_rhostar;

        /* fast implementation... minimised calls to 'pow' routine... */
        Scalar Troot = std::sqrt(Tbar);
        Scalar Tpow = Troot;
        Scalar lam = 0;

        for(int k = 0; k < thcond_a_count; ++k) {
            lam += thcond_a[k] * Tpow;
            Tpow *= Tbar;
        }

        lam += 
            thcond_b0 + thcond_b1
            * rhobar + thcond_b2
            * std::exp(thcond_B1 * ((rhobar + thcond_B2)*(rhobar + thcond_B2)));

        Scalar DTbar = std::abs(Tbar - 1) + thcond_c4;
        Scalar DTbarpow = std::pow(DTbar, 3./5);
        Scalar Q = 2. + thcond_c5 / DTbarpow;

        Scalar S;
        if(Tbar >= 1){
            S = 1. / DTbar;
        }else{
            S = thcond_c6 / DTbarpow;
        }

        Scalar rhobar18 = std::pow(rhobar, 1.8);
        Scalar rhobarQ = std::pow(rhobar, Q);

        lam +=
            (thcond_d1 / std::pow(Tbar,10) + thcond_d2) * rhobar18 *
            std::exp(thcond_c1 * (1 - rhobar * rhobar18))
            + thcond_d3 * S * rhobarQ *
            std::exp((Q/(1+Q))*(1 - rhobar*rhobarQ))
            + thcond_d4 *
            std::exp(thcond_c2 * std::pow(Troot,3) + thcond_c3 / std::pow(rhobar,5));
        return /*thcond_kstar * */ lam;
    }
};

} // end namepace IAPWS
} // end namepace Dune

#endif
