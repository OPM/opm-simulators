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
 * \copydoc Opm::BinaryCoeff::fullerMethod
 */
#ifndef OPM_FULLERMETHOD_HPP
#define OPM_FULLERMETHOD_HPP

#include <opm/material/common/Means.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Estimate binary diffusion coefficents \f$\mathrm{[m^2/s]}\f$ in gases according to
 *        the method by Fuller.
 *
 * \param M molar masses \f$\mathrm{[g/mol]}\f$
 * \param SigmaNu atomic diffusion volume
 * \param temperature The temperature \f$\mathrm{[K]}\f$
 * \param pressure phase pressure \f$\mathrm{[Pa]}\f$
 *
 * This function estimates the diffusion coefficents in binary gases
 * using to the method proposed by Fuller. This method and is only
 * valid at "low" pressures.
 *
 * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
 * edition, McGraw-Hill, 1987, pp. 587-588
 */
template <class Scalar, class Evaluation = Scalar>
inline Evaluation fullerMethod(const Scalar* M, // molar masses [g/mol]
                               const Scalar* SigmaNu, // atomic diffusion volume
                               const Evaluation& temperature, // [K]
                               const Evaluation& pressure) // [Pa]
{
    // "effective" molar mass in [g/m^3]
    Scalar Mab = harmonicMean(M[0], M[1]);

    // Fuller's method
    const Evaluation& tmp = std::pow(SigmaNu[0], 1./3) + std::pow(SigmaNu[1], 1./3);
    return 1e-4 * (143.0*pow(temperature, 1.75))/(pressure*std::sqrt(Mab)*tmp*tmp);
}

} // namespace BinaryCoeff
} // namespace Opm

#endif
