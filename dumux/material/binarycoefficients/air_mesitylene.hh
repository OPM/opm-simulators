// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief Binary coefficients for air and mesitylene.
 */
#ifndef DUMUX_BINARY_COEFF_AIR_MESITYLENE_HH
#define DUMUX_BINARY_COEFF_AIR_MESITYLENE_HH

#include <dumux/material/components/air.hh>
#include <dumux/material/components/mesitylene.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and mesitylene.
 */
class Air_Mesitylene
{
public:
    /*!
     *
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    { DUNE_THROW(Dune::NotImplemented,
                 "Henry coefficient of air in mesitylene");
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for air and mesitylene.
     * I used the method according to Wilke and Lee
     * see Handbook of chem. property's Estimation Methods
     * W.J. Lyman, W.F. Reehl, D.H. Rosenblatt
     *
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::Air<Scalar> Air;
        typedef Dumux::Mesitylene<Scalar> Mesitylene;

        temperature = std::max(temperature, 1e-9); // regularization
        temperature = std::min(temperature, 500.0); // regularization
        pressure = std::max(pressure, 0.0); // regularization
        pressure = std::min(pressure, 1e8); // regularization

        const Scalar M_m = 1e3*Mesitylene::molarMass(); // [g/mol] molecular weight of mesitylene
        const Scalar M_a = 1e3*Air::molarMass(); // [g/mol] molecular weight of air
        const Scalar Tb_m = 437.9;        // [K] boiling temperature of mesitylene
        const Scalar sigma_a = 3.711;     // charact. length of air
        const Scalar T_scal_a = 78.6;     // [K] (molec. energy of attraction/Boltzmann constant)
        const Scalar V_B_m = 162.6;       // [cm^3/mol] LeBas molal volume of mesitylene
        const Scalar sigma_m = 1.18*std::pow(V_B_m, 0.333);     // charact. length of mesitylene
        const Scalar sigma_am = 0.5*(sigma_a + sigma_m);
        const Scalar T_scal_m = 1.15*Tb_m;
        const Scalar T_scal_am = std::sqrt(T_scal_a*T_scal_m);

        Scalar T_star = temperature/T_scal_am;
        T_star = std::max(T_star, 1e-5); // regularization

        const Scalar Omega = 1.06036/std::pow(T_star, 0.1561) + 0.193/std::exp(T_star*0.47635)
            + 1.03587/std::exp(T_star*1.52996) + 1.76474/std::exp(T_star*3.89411);
        const Scalar B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_a + 1.0/M_m);
        const Scalar Mr = (M_a + M_m)/(M_a*M_m);
        const Scalar D_am = (B_*std::pow(temperature, 1.5) * std::sqrt(Mr))
                           /(1e-5*pressure*std::pow(sigma_am, 2.0) * Omega); // [cm^2/s]

        return 1e-4*D_am; // [m^2/s]
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular mesitylene in liquid water.
     *
     * \todo
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented,
                 "Binary liquid diffusion coefficients of air and mesitylene");
    }
};

}
} // end namepace

#endif
