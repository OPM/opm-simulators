// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
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
 * \brief Binary coefficients for water and xylene.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_XYLENE_HH
#define DUMUX_BINARY_COEFF_H2O_XYLENE_HH

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/xylene.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and xylene.
 */
class H2O_Xylene
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for xylene in liquid water.
     *
     * See:
     *
     *  Sanders1999 Henry collection
     */

    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        // after Sanders
        Scalar sanderH = 1.5e-1;    //[M/atm]
        //conversion to our Henry definition
        Scalar dumuxH = sanderH / 101.325; // has now [(mol/m^3)/Pa]
        dumuxH *= 18.02e-6;  //multiplied by molar volume of reference phase = water
        return 1.0/dumuxH; // [Pa]
    };

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and xylene.
     *
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::H2O<Scalar> H2O;
        typedef Dumux::Xylene<Scalar> Xylene;

        temperature = std::max(temperature, 1e-9); // regularization
        temperature = std::min(temperature, 500.0); // regularization
        pressure = std::max(pressure, 0.0); // regularization
        pressure = std::min(pressure, 1e8); // regularization

        const Scalar M_x = 1e3*Xylene::molarMass(); // [g/mol] molecular weight of xylene
        const Scalar M_w = 1e3*H2O::molarMass(); // [g/mol] molecular weight of water
        const Scalar Tb_x = 412.9;        // [K] boiling temperature of xylene
        const Scalar Tb_w = 373.15;       // [K] boiling temperature of water (at p_atm)
        const Scalar V_B_w = 18.0;                // [cm^3/mol] LeBas molal volume of water
        const Scalar  sigma_w = 1.18*std::pow(V_B_w, 0.333);     // charact. length of air
        const Scalar  T_scal_w = 1.15*Tb_w;     // [K] (molec. energy of attraction/Boltzmann constant)
        const Scalar V_B_x = 140.4;       // [cm^3/mol] LeBas molal volume of xylene
        const Scalar sigma_x = 1.18*std::pow(V_B_x, 0.333);     // charact. length of xylene
        const Scalar sigma_wx = 0.5*(sigma_w + sigma_x);
        const Scalar T_scal_x = 1.15*Tb_x;
        const Scalar T_scal_wx = std::sqrt(T_scal_w*T_scal_x);

        Scalar T_star = temperature/T_scal_wx;
        T_star = std::max(T_star, 1e-5); // regularization

        const Scalar Omega = 1.06036/std::pow(T_star, 0.1561) + 0.193/std::exp(T_star*0.47635)
            + 1.03587/std::exp(T_star*1.52996) + 1.76474/std::exp(T_star*3.89411);
        const Scalar  B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_w + 1.0/M_x);
        const Scalar Mr = (M_w + M_x)/(M_w*M_x);
        const Scalar D_wx = (B_*std::pow(temperature,1.6)*std::sqrt(Mr))
                           /(1e-5*pressure*std::pow(sigma_wx, 2.0)*Omega); // [cm^2/s]

        return D_wx*1e-4; // [m^2/s]
    };

    /*!
     * \brief Diffusion coefficent [m^2/s] for xylene in liquid water.
     *
     * \todo
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 1.e-9;  // This is just an order of magnitude. Please improve it!
    };
};

}
} // end namepace

#endif
