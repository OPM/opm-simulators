// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
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
 * \copydoc Ewoms::IAPWS::Region4
 */
#ifndef EWOMS_IAPWS_REGION4_HH
#define EWOMS_IAPWS_REGION4_HH

#include <cmath>

namespace Ewoms {
namespace IAPWS {

/*!
 * \ingroup IAPWS
 *
 * \brief Implements the equations for region 4 of the IAPWS '97 formulation.
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
class Region4
{
public:
    /*!
     * \brief Returns the saturation pressure in \f$\mathrm{[Pa]}\f$ of pure water at a given
     *        temperature.
     *
     *\param temperature temperature of component in \f$\mathrm{[K]}\f$
     *
     * The saturation pressure is often also called vapor pressure.
     */
    static Scalar saturationPressure(Scalar temperature)
    {
        static constexpr Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        Scalar sigma = temperature + n[8]/(temperature - n[9]);

        Scalar A = (sigma + n[0])*sigma + n[1];
        Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        Scalar tmp = 2*C/(std::sqrt(B*B - 4*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Returns the saturation temperature in \f$\mathrm{[K]}\f$ of pure water at a given
     *        pressure.
     *
     *\param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * The saturation pressure is often also called vapor pressure.
     */
    static Scalar vaporTemperature(Scalar pressure)
    {
        static const Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };
        Scalar beta = pow((pressure/1e6 /*from Pa to MPa*/), (1./4.));
        Scalar beta2 = pow(beta, 2.);
        Scalar E = beta2 + n[2] * beta + n[5];
        Scalar F = n[0]*beta2 + n[3]*beta + n[6];
        Scalar G = n[1]*beta2 + n[4]*beta + n[7];

        Scalar D = ( 2.*G)/(-F -std::sqrt(pow(F,2.) - 4.*E*G));

        Scalar temperature = (n[9] + D - std::sqrt(pow(n[9]+D , 2.) - 4.* (n[8] + n[9]*D)) ) * 0.5;

        return temperature;
    }
};

} // end namepace IAPWS
} // end namepace Dune

#endif
