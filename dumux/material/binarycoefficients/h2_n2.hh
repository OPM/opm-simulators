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
 *
 * \brief Binary coefficients for hydrogen and nitrogen.
 */
#ifndef DUMUX_BINARY_COEFF_H2_N2_HH
#define DUMUX_BINARY_COEFF_H2_N2_HH

#include "henryiapws.hh"
#include "fullermethod.hh"

#include <dumux/material/components/h2.hh>
#include <dumux/material/components/n2.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for hydrogen and nitrogen.
 */
class H2_N2
{
public:
    /*!
     * \brief Henry coefficent \f$\mathrm{[N/m^2]}\f$ for molecular nitrogen in liquid hydrogen.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        DUNE_THROW(Dune::NotImplemented, "henry coefficient for nitrogen in liquid hydrogen");
    };

    /*!
     * \brief Binary diffusion coefficent \f$\mathrm{[m^2/s]}\f$ for molecular hydrogen and nitrogen.
     *
     * This function estimates the diffusion coefficents in binary gases
     * using to the method proposed by Fuller. This method and is only
     * valid at "low" pressures.
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp. 587-588
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::H2<Scalar> H2;
        typedef Dumux::N2<Scalar> N2;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 6.12 /* H2 */,  18.5 /* N2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2::molarMass()*1e3, N2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    };

    /*!
     * \brief Diffusion coefficent \f$\mathrm{[m^2/s]}\f$ for molecular nitrogen in liquid hydrogen.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for liquid nitrogen and hydrogen");
    };
};

}
} // end namepace

#endif
