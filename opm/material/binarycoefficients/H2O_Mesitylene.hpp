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
 * \copydoc Opm::BinaryCoeff::H2O_Mesitylene
 */
#ifndef OPM_BINARY_COEFF_H2O_MESITYLENE_HPP
#define OPM_BINARY_COEFF_H2O_MESITYLENE_HPP

#include <opm/material/components/H2O.hpp>
#include <opm/material/components/Mesitylene.hpp>

namespace Opm
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and mesitylene.
 */
class H2O_Mesitylene
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for mesitylene in liquid water.
     *
     * See:
     *
     *  Sanders1999 Henry collection
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& /*temperature*/)
    {
        // after Sanders
        double sanderH = 1.7e-1; // [M/atm]
        //conversion to our Henry definition
        double opmH = sanderH / 101.325; // has now [(mol/m^3)/Pa]
        opmH *= 18.02e-6; // multiplied by molar volume of reference phase = water
        return 1.0/opmH; // [Pa]
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and mesitylene.
     *
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(Evaluation temperature, Evaluation pressure)
    {
        typedef H2O<double> H2O;
        typedef Mesitylene<double> Mesitylene;

        temperature = max(temperature, 1e-9); // regularization
        temperature = min(temperature, 500.0); // regularization
        pressure = max(pressure, 0.0); // regularization
        pressure = min(pressure, 1e8); // regularization

        const double M_m = 1e3*Mesitylene::molarMass(); // [g/mol] molecular weight of mesitylene
        const double M_w = 1e3*H2O::molarMass(); // [g/mol] molecular weight of water
        const double Tb_m = 437.9;        // [K] boiling temperature of mesitylen
        const double Tb_w = 373.15;       // [K] boiling temperature of water (at p_atm)
        const double V_B_w = 18.0;                // [cm^3/mol] LeBas molal volume of water
        const double sigma_w = 1.18*std::pow(V_B_w, 0.333);     // charact. length of air
        const double T_scal_w = 1.15*Tb_w;     // [K] (molec. energy of attraction/Boltzmann constant)
        const double V_B_m = 162.6;       // [cm^3/mol] LeBas molal volume of mesitylen
        const double sigma_m = 1.18*std::pow(V_B_m, 0.333);     // charact. length of mesitylen
        const double sigma_wm = 0.5*(sigma_w + sigma_m);
        const double T_scal_m = 1.15*Tb_m;
        const double T_scal_wm = std::sqrt(T_scal_w*T_scal_m);

        const Evaluation T_star = max(temperature/T_scal_wm, 1e-5);

        const Evaluation& Omega = 1.06036/pow(T_star,0.1561) + 0.193/exp(T_star*0.47635)
            + 1.03587/exp(T_star*1.52996) + 1.76474/exp(T_star*3.89411);
        const double B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_w + 1.0/M_m);
        const double Mr = (M_w + M_m)/(M_w*M_m);
        const Evaluation D_wm = (B_*pow(temperature,1.6)*std::sqrt(Mr))
                           /(1e-5*pressure*std::pow(sigma_wm, 2.0)*Omega); // [cm^2/s]

        return D_wm*1e-4;   //  [m^2/s]
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for mesitylene in liquid water.
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        // This is just an order of magnitude estimate. Please improve!
        return 1e-9;
    }
};

} // namespace BinaryCoeff
} // namespace Opm

#endif
