/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
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
 * \brief Relations valid for an ideal gas.
 */
#ifndef DUMUX_IDEAL_GAS_HH
#define DUMUX_IDEAL_GAS_HH

namespace Dumux
{

/*!
 * \brief Relations valid for an ideal gas.
 */
template <class Scalar>
class IdealGas
{
public:
    //! The ideal gas constant [J/mol/K]
    static const Scalar R = 8.3144;
    
    /*!
     * \brief The density of the gas in [kg/m^3], depending on
     *        pressure, temperature and average molar mass of the gas.
     */
    static Scalar density(Scalar avgMolarMass, 
                          Scalar temperature,
                          Scalar pressure)
    { return pressure*avgMolarMass/(R*temperature); } 

    /*!
     * \brief The pressure of the gas in [N/m^2], depending on
     *        concentration and temperature.
     */
    static Scalar pressure(Scalar temperature,
                           Scalar concentration)
    { return R*temperature*concentration; } 

    /*!
     * \brief The molar concentration of the gas in [mol/m^3], depending on
     *        pressure and temperature.
     */
    static Scalar concentration(Scalar temperature,
                                Scalar pressure)
    { return pressure/(R*temperature); } 
};

} // end namepace

#endif
