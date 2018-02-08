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
 * \copydoc Opm::BinaryCoeff::Air_Xylene
 */
#ifndef OPM_BINARY_COEFF_AIR_XYLENE_HPP
#define OPM_BINARY_COEFF_AIR_XYLENE_HPP

#include <opm/material/components/Air.hpp>
#include <opm/material/components/Xylene.hpp>

namespace Opm {
namespace BinaryCoeff {

/*!
 * \brief Binary coefficients for water and xylene.
 */
class Air_Xylene
{
public:
    /*!
     *
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& /*temperature*/)
    { throw std::runtime_error("Not implemented: Henry coefficient of air in xylene"); }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for air and xylene.
     * method according to Wilke and Lee
     * see Handbook of chem. property's Estimation Methods
     * W.J. Lyman, W.F. Reehl, D.H. Rosenblatt
     *
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(Evaluation temperature, Evaluation pressure)
    {
        typedef Opm::Air<double> Air;
        typedef Opm::Xylene<double> Xylene;

        temperature = Opm::max(temperature, 1e-9); // regularization
        temperature = Opm::min(temperature, 500.0); // regularization
        pressure = Opm::max(pressure, 0.0); // regularization
        pressure = Opm::min(pressure, 1e8); // regularization

        const double M_x = 1e3*Xylene::molarMass(); // [g/mol] molecular weight of xylene
        const double M_a = 1e3*Air::molarMass(); // [g/mol] molecular weight of air
        const double Tb_x = 412.0;        // [K] boiling temperature of xylene
        const double sigma_a = 3.711;     // charact. length of air
        const double T_scal_a = 78.6;     // [K] (molec. energy of attraction/Boltzmann constant)
        const double V_B_x = 140.4;       // [cm^3/mol] LeBas molal volume of xylene
        const double sigma_x = 1.18*std::pow(V_B_x, 0.333);     // charact. length of xylene
        const double sigma_ax = 0.5*(sigma_a + sigma_x);
        const double T_scal_x = 1.15*Tb_x;
        const double T_scal_ax = std::sqrt(T_scal_a*T_scal_x);

        const Evaluation& T_star = Opm::max(temperature/T_scal_ax, 1e-5);

        const Evaluation& Omega = 1.06036/Opm::pow(T_star, 0.1561) + 0.193/Opm::exp(T_star*0.47635)
            + 1.03587/Opm::exp(T_star*1.52996) + 1.76474/Opm::exp(T_star*3.89411);
        const double B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_a + 1.0/M_x);
        const double Mr = (M_a + M_x)/(M_a*M_x);
        return 1e-4
            *(B_*Opm::pow(temperature,1.5)*std::sqrt(Mr))
            /(1e-5*pressure*std::pow(sigma_ax, 2.0)*Omega);
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular xylene in liquid water.
     *
     * \todo
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { throw std::runtime_error("Not implemented: Binary liquid diffusion coefficients of air and xylene"); }
};

} // namespace BinaryCoeff
} // namespace Opm

#endif
