// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Various relations for molecular diffusion coefficients.
 */
#ifndef DUMUX_FULLERMETHOD_HH
#define DUMUX_FULLERMETHOD_HH

namespace Dumux
{
namespace BinaryCoeff
{
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
template <class Scalar>
inline Scalar fullerMethod(const Scalar *M, // molar masses [g/mol]
                           const Scalar *SigmaNu, // atomic diffusion volume
                           const Scalar temperature, // [K]
                           const Scalar pressure) // [Pa]
{
    // "effective" molar mass in [g/m^3]
    Scalar Mab = harmonicMean(M[0], M[1]);

    // Fuller's method
    Scalar tmp = pow(SigmaNu[0], 1./3) + pow(SigmaNu[1], 1./3);
    return 1e-4 * (143.0*pow(temperature, 1.75))/(pressure*sqrt(Mab)*tmp*tmp);
}

} // end namepace BinaryCoeff
} // end namepace Dumux

#endif
