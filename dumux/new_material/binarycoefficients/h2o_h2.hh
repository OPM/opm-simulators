/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Binary coefficients for water and hydrogen.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_H2_HH
#define DUMUX_BINARY_COEFF_H2O_H2_HH

#include "henryiapws.hh"
#include "fullermethod.hh"

#include <dumux/new_material/components/h2.hh>
#include <dumux/new_material/components/h2o.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and hydrogen.
 */
class H2O_H2
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular hydrogen in liquid water.
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
        const Scalar E = 2286.4159;
        const Scalar F = 11.3397;
        const Scalar G = -70.7279;
        const Scalar H = 63.0631;

        return henryIAPWS(E, F, G, H, temperature);
    };

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular hydrogen and water steam.
     *
     * \copybody fullerMethod()
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::H2<Scalar> H2;
        typedef Dumux::H2O<Scalar> H2O;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 6.12 /* H2 */,  13.1 /* H2O */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2::molarMass()*1e3, H2O::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    };

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular nitrogen in liquid water.
     *
     * The empirical equations for estimating the diffusion coefficient in
     * infinite solution which are presented in Reid, 1987 all show a
     * linear dependency on temperature. We thus simply scale the
     * experimentally obtained diffusion coefficient of Ferrell and
     * Himmelblau by the temperature.
     *
     * See:
     *
     * R. Reid et al.: "The properties of Gases and Liquids", 4th edition,
     * pp. 599, McGraw-Hill, 1987
     *
     * R. Ferrell, D. Himmelblau: "Diffusion Coeffients of Nitrogen and
     * Oxygen in Water", Journal of Chemical Engineering and Data,
     * Vol. 12, No. 1, pp. 111-115, 1967
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        const Scalar Texp = 273.15 + 25; // [K]
        const Scalar Dexp = 4.5e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    };
};

}
} // end namepace

#endif
