// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::BinaryCoeff::H2O_CO2
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CO2_HH
#define DUMUX_BINARY_COEFF_H2O_CO2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleco2.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and CO2.
 */
class H2O_CO2
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular CO2 in liquid water.
     *
     * See:
     *
     * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
     * Distribution Constant for Gases in H2O and D2O at High
     * Temperatures"
     * http://www.iapws.org/relguide/HenGuide.pdf
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        const Scalar E = 1672.9376;
        const Scalar F = 28.1751;
        const Scalar G = -112.4619;
        const Scalar H = 85.3807;

        return henryIAPWS(E, F, G, H, temperature);
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and CO2.
     *
     * To calculate the values, the \ref fullerMethod is used.
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::H2O<Scalar> H2O;
        typedef Dumux::SimpleCO2<Scalar> CO2;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  26.9 /* CO2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, CO2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular CO2 in liquid water.
     *
     * \todo
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Binary liquid diffusion coefficients of CO2 and CH4");
    }
};

}
} // end namepace

#endif
