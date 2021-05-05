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
 * \copydoc Opm::IAPWS::Region4
 */
#ifndef OPM_IAPWS_REGION4_HPP
#define OPM_IAPWS_REGION4_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
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
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     *
     * The saturation pressure is often also called vapor pressure.
     */
    template <class Evaluation>
    static Evaluation saturationPressure(const Evaluation& temperature)
    {
        static const Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        const Evaluation& sigma = temperature + n[8]/(temperature - n[9]);

        const Evaluation& A = (sigma + n[0])*sigma + n[1];
        const Evaluation& B = (n[2]*sigma + n[3])*sigma + n[4];
        const Evaluation& C = (n[5]*sigma + n[6])*sigma + n[7];

        Evaluation tmp = 2*C/(sqrt(B*B - 4*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Returns the saturation temperature in \f$\mathrm{[K]}\f$ of pure water at a given
     *        pressure.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * The saturation pressure is often also called vapor pressure.
     */
    template <class Evaluation>
    static Evaluation vaporTemperature(const Evaluation& pressure)
    {
        static const Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };
        const Evaluation& beta = pow((pressure/1e6 /*from Pa to MPa*/), (1./4.));
        const Evaluation& beta2 = pow(beta, 2.);
        const Evaluation& E = beta2 + n[2] * beta + n[5];
        const Evaluation& F = n[0]*beta2 + n[3]*beta + n[6];
        const Evaluation& G = n[1]*beta2 + n[4]*beta + n[7];

        const Evaluation& D = ( 2.*G)/(-F -sqrt(pow(F,2.) - 4.*E*G));

        const Evaluation& temperature = (n[9] + D - sqrt(pow(n[9]+D , 2.) - 4.* (n[8] + n[9]*D)) ) * 0.5;

        return temperature;
    }
};

} // namespace IAPWS
} // namespace Opm

#endif
