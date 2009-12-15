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
 * \brief Binary coefficients for water and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_O2_HH
#define DUMUX_BINARY_COEFF_H2O_O2_HH

#include "henryiapws.hh"
#include "fullermethod.hh"

#include <dumux/new_material/components/o2.hh>
#include <dumux/new_material/components/h2o.hh>

namespace Dune
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and oxygen.
 */
class H2O_O2
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular oxygen in liquid water.
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
        const Scalar E =  2305.0674;
        const Scalar F = -11.3240;
        const Scalar G =  25.3224;
        const Scalar H = -15.6449;
        
        return henryIAPWS(E, F, G, H, temperature);
    };

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and oxygen.
     *
     * \copybody fullerMethod()
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dune::H2O<Scalar> H2O;
        typedef Dune::O2<Scalar>  O2;
        
        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  16.3 /* O2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, O2::molarMass()*1e3 };
        
        return fullerMethod(M, SigmaNu, temperature, pressure);
    };

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular oxygen in liquid water.
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
        const Scalar Dexp = 2.2e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    };
};

}
} // end namepace

#endif
