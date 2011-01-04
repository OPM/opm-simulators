// $Id$
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
 *
 * \ingroup IAPWS
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

#include <cmath>
#include <iostream>

namespace Dumux
{
namespace IAPWS
{

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
    static const Scalar molarMass = 18.01518e-3;

    //! Specific gas constant of water \f$\mathrm{[J/(kg*K)]}\f$
    static const Scalar R = 461.526;

    //! Critical temperature of water \f$\mathrm{[K]}\f$
    static const Scalar criticalTemperature = 647.096;

    //! Critical pressure of water \f$\mathrm{[Pa]}\f$
    static const Scalar criticalPressure = 22.064e6;

    //! Critical molar volume of water \f$\mathrm{[m^3/mol]}\f$
    static const Scalar criticalMolarVolume = molarMass/322.0;

    //! The acentric factor of water \f$\mathrm{[-]}\f$
    static const Scalar acentricFactor = 0.344;

    //! Density of water at the critical point \f$\mathrm{[kg/m^3]}\f$
    static const Scalar criticalDensity = 322;

    //! Triple temperature of water \f$\mathrm{[K]}\f$
    static const Scalar tripleTemperature = 273.16;

    //! Triple pressure of water \f$\mathrm{[Pa]}\f$
    static const Scalar triplePressure = 611.657;

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
            {-1.08374 , 1.88797 ,-7.72479e-1, 0, 0, 0, 0 },
            {-2.89555e-1, 1.26613 ,-4.89837e-1, 0, 6.98452e-2, 0,-4.35673e-3 },
            {          0, 0,-2.57040e-1, 0, 0, 8.72102e-3, 0 },
            {          0, 1.20573e-1, 0, 0, 0, 0,-5.93264e-4 }
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
};

} // end namepace IAPWS
} // end namepace Dune

#endif
