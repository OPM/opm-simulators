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
 *
 * \copydoc Opm::BinaryCoeff::H2O_CO2
 */
#ifndef OPM_BINARY_COEFF_H2O_CO2_HPP
#define OPM_BINARY_COEFF_H2O_CO2_HPP

#include <opm/material/binarycoefficients/HenryIapws.hpp>
#include <opm/material/binarycoefficients/FullerMethod.hpp>

#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleCO2.hpp>

namespace Opm {
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
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation henry(const Evaluation& temperature)
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
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        typedef Opm::H2O<Scalar> H2O;
        typedef Opm::SimpleCO2<Scalar> CO2;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  26.9 /* CO2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, CO2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular CO2 in liquid water.
     */
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { throw std::runtime_error("Not implemented: Binary liquid diffusion coefficients of CO2 and CH4"); }
};

} // namespace BinaryCoeff
} // namespace Opm

#endif
