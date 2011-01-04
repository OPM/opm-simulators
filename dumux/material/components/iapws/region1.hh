// $Id$
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
 *\ingroup IAPWS
 *
 * \brief Implements the equations for region 1 of the IAPWS '97 formulation.
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
#ifndef DUMUX_IAPWS_REGION1_HH
#define DUMUX_IAPWS_REGION1_HH

#include <cmath>
#include <iostream>

namespace Dumux
{
namespace IAPWS
{
/*!
 * \ingroup IAPWS
 *
 * \brief Implements the equations for region 1 of the IAPWS '97 formulation.
 *
 * \tparam Scalar The type used for scalar values
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class Region1
{
public:
    /*!
     * \brief Returns true if IAPWS region 1 applies for a
     *        (temperature in \f$\mathrm{[K]}\f$, pressure in \f$\mathrm{[Pa]}\f$) pair.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static bool isValid(Scalar temperature, Scalar pressure)
    {
        return
            temperature <= 623.15 &&
            pressure <= 100e6;

        // actually this is:
        /*
        return
           (
           273.15 <= temperature &&
           temperature <= 623.15 &&
           pressure >= vaporPressure(temperature) &&
           pressure <= 100e6
           );
        */
    }

    /*!
     * \brief Returns the reduced temperature for IAPWS region 1.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar tau(Scalar temperature)
    { return 1386.0 / temperature; }

    /*!
     * \brief Returns the derivative of the reduced temperature to the
     *        temperature for IAPWS region 1 in \f$\mathrm{[1/K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar dtau_dT(Scalar temperature)
    { return - 1386.0 / (temperature*temperature); }

    /*!
     * \brief Returns the reduced pressure for IAPWS region 1.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar pi(Scalar pressure)
    { return pressure / 16.53e6; }

    /*!
     * \brief Returns the derivative of the reduced pressure to the
     *        pressure for IAPWS region 1 in \f$\mathrm{[1/Pa]}\f$.
     *
     * \param pressure temperature of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dpi_dp(Scalar pressure)
    { return 1.0 / 16.53e6; }

    /*!
     * \brief Returns the derivative of the pressure to the
     *        reduced pressure for IAPWS region 1 in \f$\mathrm{[Pa]}\f$.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dp_dpi(Scalar pressure)
    { return 16.53e6; }


    /*!
     * \brief The Gibbs free energy (dimensionless) for IAPWS region 1 (i.e. liquid)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gamma(Scalar temperature, Scalar pressure)
    {
        Scalar tau_ = tau(temperature);   /* reduced temperature */
        Scalar pi_ = pi(pressure);    /* reduced pressure */

        Scalar result = 0;
        for (int i = 0; i < 34; ++i) {
            result += n(i)*pow(7.1 - pi_, I(i))*pow(tau_ - 1.222, J(i));
        }

        return result;
    }


    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized temperature for IAPWS region 1 (i.e. liquid) (dimensionless).
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar dgamma_dtau(Scalar temperature, Scalar pressure)
    {
        Scalar tau_ = tau(temperature);   /* reduced temperature */
        Scalar pi_ = pi(pressure);    /* reduced pressure */

        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                n(i) *
                std::pow(7.1 - pi_, I(i)) *
                std::pow(tau_ - 1.222,  J(i)-1) *
                J(i);
        }

        return result;
    }

    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized pressure for IAPWS region 1 (i.e. liquid) dimensionless).
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar dgamma_dpi(Scalar temperature, Scalar pressure)
    {
        Scalar tau_ = tau(temperature);   /* reduced temperature */
        Scalar pi_ = pi(pressure);    /* reduced pressure */

        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                -n(i) *
                I(i) *
                std::pow(7.1 - pi_, I(i) - 1) *
                std::pow(tau_ - 1.222, J(i));
        }

        return result;
    }

    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized pressure and to the normalized temperature
     *        for IAPWS region 1 (i.e. liquid water) (dimensionless).
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar ddgamma_dtaudpi(Scalar temperature, Scalar pressure)
    {
        Scalar tau_ = tau(temperature);   /* reduced temperature */
        Scalar pi_ = pi(pressure);    /* reduced pressure */

        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                -n(i) *
                I(i) *
                J(i) *
                std::pow(7.1 - pi_, I(i) - 1) *
                std::pow(tau_ - 1.222, J(i) - 1);
        }

        return result;
    }

    /*!
     * \brief The second partial derivative of the Gibbs free energy
     *        to the normalized pressure for IAPWS region 1
     *        (i.e. liquid water) (dimensionless).
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar ddgamma_ddpi(Scalar temperature, Scalar pressure)
    {
        Scalar tau_ = tau(temperature);   /* reduced temperature */
        Scalar pi_ = pi(pressure);    /* reduced pressure */

        Scalar result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                n(i) *
                I(i) *
                (I(i) - 1) *
                std::pow(7.1 - pi_, I(i) - 2) *
                std::pow(tau_ - 1.222, J(i));
        }

        return result;
    }

private:
    static Scalar n(int i)
    {
        static const Scalar n[34] = {
            0.14632971213167, -0.84548187169114, -0.37563603672040e1,
            0.33855169168385e1, -0.95791963387872, 0.15772038513228,
           -0.16616417199501e-1, 0.81214629983568e-3, 0.28319080123804e-3,
           -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,
           -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,
           -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5,
           -0.72694996297594e-15,-0.31679644845054e-4, -0.28270797985312e-5,
           -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6,
           -0.14341729937924e-12,-0.40516996860117e-6, -0.12734301741641e-8,
           -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,
            0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,
           -0.93537087292458e-25
        };
        return n[i];
    }

    static short int I(int i)
    {
        static const short int I[34] = {
            0, 0, 0,
            0, 0, 0,
            0, 0, 1,
            1, 1, 1,
            1, 1, 2,
            2, 2, 2,
            2, 3, 3,
            3, 4, 4,
            4, 5, 8,
            8, 21, 23,
            29, 30, 31,
            32
        };
        return I[i];
    }

    static short int J(int i)
    {
        static const short int J[34] = {
             -2, -1, 0,
              1, 2, 3,
              4, 5, -9,
             -7, -1, 0,
              1, 3, -3,
              0, 1, 3,
             17, -4, 0,
              6, -5, -2,
             10, -8, -11,
             -6, -29, -31,
            -38, -39, -40,
            -41
        };
        return J[i];
    }

};

} // end namepace IAPWS
} // end namepace Dune

#endif
