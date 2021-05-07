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
 * \copydoc Opm::BinaryCoeff::H2O_N2
 */
#ifndef OPM_BINARY_COEFF_H2O_N2_HPP
#define OPM_BINARY_COEFF_H2O_N2_HPP

#include "HenryIapws.hpp"
#include "FullerMethod.hpp"

#include <opm/material/components/N2.hpp>
#include <opm/material/components/H2O.hpp>

namespace Opm {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and nitrogen.
 */
class H2O_N2
{
public:
    /*!
     * \brief Henry coefficent \f$\mathrm{[N/m^2]}\f$  for molecular nitrogen in liquid water.
     *
     * \copydetails Opm::henryIAPWS
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& temperature)
    {
        const double E = 2388.8777;
        const double F = -14.9593;
        const double G = 42.0179;
        const double H = -29.4396;

        return henryIAPWS(E, F, G, H, temperature);
    }

    /*!
     * \brief Binary diffusion coefficent \f$\mathrm{[m^2/s]}\f$ for molecular water and nitrogen.
     *
     * Uses fullerMethod to determine the diffusion of water in nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        typedef H2O<double> H2O;
        typedef N2<double> N2;

        // atomic diffusion volumes
        const double SigmaNu[2] = { 13.1 /* H2O */,  18.5 /* N2 */ };
        // molar masses [g/mol]
        const double M[2] = { H2O::molarMass()*1e3, N2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
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

        return Dexp * temperature/Texp;
    }
};

} // namespace BinaryCoeff
} // namespace Opm

#endif
