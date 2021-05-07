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
 * \copydoc Opm::BinaryCoeff::H2O_Air
 */
#ifndef OPM_BINARY_COEFF_H2O_AIR_HPP
#define OPM_BINARY_COEFF_H2O_AIR_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and nitrogen.
 */
class H2O_Air
{
public:
    /*!
     * \brief Henry coefficent \f$\mathrm{[N/m^2]}\f$  for air in liquid water.
     *
     *
     * Henry coefficent See:
     * Stefan Finsterle, 1993
     * Inverse Modellierung zur Bestimmung hydrogeologischer Parameter eines Zweiphasensystems
     * (fitted to data from the book "Tchobanoglous & Schroeder: Water Quality:
     * Characteristics", Addison-Wesley, 1985)
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& temperature)
    { return 1.0/((0.8942+1.47*exp(-0.04394*(temperature-273.15)))*1e-10); }

    /*!
     * \brief Binary diffusion coefficent \f$\mathrm{[m^2/s]}\f$ for molecular water and air
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     * Vargaftik : Tables on the thermophysical properties of liquids and gases. John Wiley &      * Sons, New York, 1975.
     *
     * Walker, Sabey, Hampton: Studies of heat transfer and water migration in soils.
     * Dep. of Agricultural and Chemical Engineering, Colorado State University,
     * Fort Collins, 1981.
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        double Theta=1.8;
        double Daw=2.13e-5;  /* reference value */
        double pg0=1.e5;     /* reference pressure */
        double T0=273.15;    /* reference temperature */

        return Daw*(pg0/pressure)*pow((temperature/T0),Theta);
    }

    /*!
     * Lacking better data on water-air diffusion in liquids, we use at the
     * moment the diffusion coefficient of the air's main component nitrogen!!
     * \brief Diffusion coefficent \f$\mathrm{[m^2/s]}\f$ for molecular nitrogen in liquid water.
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
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        const double Texp = 273.15 + 25; // [K]
        const double Dexp = 2.01e-9; // [m^2/s]
        return Dexp/Texp*temperature;
    }
};

} // namespace BinaryCoeff
} // namespace Opm

#endif
