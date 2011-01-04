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
 * \ingroup IAPWS
 *
 * \brief Implements the equations for region 2 of the IAPWS '97 formulation.
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
#ifndef DUMUX_IAPWS_REGION2_HH
#define DUMUX_IAPWS_REGION2_HH

#include <cmath>
#include <iostream>

namespace Dumux
{
namespace IAPWS
{
/*!
 *
 * \ingroup IAPWS
 *
 * \brief Implements the equations for region 2 of the IAPWS '97 formulation.
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
class Region2
{
public:
    /*!
     * \brief Returns true if IAPWS region 2 applies for a
     *        (temperature, pressure) pair.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static bool isValid(Scalar temperature, Scalar pressure)
    {
        return
            temperature <= 623.15 && pressure <= 100e6;

        // actually this is:
        /*
        return
           (273.15 <= temperature && temperature <= 623.15 && pressure <= vaporPressure(temperature)) ||
           (623.15 < temperature && temperature <= 863.15 && pressure <= auxPressure(temperature)) ||
           (863.15 < temperature && temperature <= 1073.15 && pressure < 100e6);
        */
    }

    /*!
     * \brief Returns the reduced temperature (dimensionless) for IAPWS region 2.
     *
     * \param temperature temperature of component
     */
    static Scalar tau(Scalar temperature)
    { return 540.0 / temperature; }

    /*!
     * \brief Returns the derivative of the reduced temperature to the
     *        temperature for IAPWS region 2.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar dtau_dT(Scalar temperature)
    { return - 540.0 / (temperature*temperature); }

    /*!
     * \brief Returns the reduced pressure (dimensionless) for IAPWS region 2.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar pi(Scalar pressure)
    { return pressure / 1e6; }

    /*!
     * \brief Returns the derivative of the reduced pressure to the
     *        pressure for IAPWS region 2 in \f$\mathrm{[1/Pa]}\f$.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dpi_dp(Scalar pressure)
    { return 1.0 / 1e6; }

    /*!
     * \brief Returns the derivative of the pressure to the
     *        reduced pressure for IAPWS region 2 (dimensionless).
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar dp_dpi(Scalar pressure)
    { return 1e6; }

    /*!
     * \brief The Gibbs free energy for IAPWS region 2 (i.e. sub-critical
     *        steam) (dimensionless).
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
        Scalar tau = tau(temperature);   /* reduced temperature */
        Scalar pi = pi(pressure);    /* reduced pressure */

        Scalar result;

        // ideal gas part
        result = ln(pi);
        for (int i = 0; i < 9; ++i)
            result += n_g(i)*std::pow(tau, J_g(i));

        // residual part
        for (int i = 0; i < 43; ++i)
            result +=
                n_r(i)*
                std::pow(pi, I_r(i))*
                std::pow(tau - 0.5, J_r(i));
        return result;
    }

    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized temperature for IAPWS region 2 (i.e. sub-critical
     *        steam) dimensionless).
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

        // ideal gas part
        Scalar result = 0;
        for (int i = 0; i < 9; i++) {
            result +=
                n_g(i) *
                J_g(i) *
                std::pow(tau_, J_g(i) - 1);
        }

        // residual part
        for (int i = 0; i < 43; i++) {
            result +=
                n_r(i) *
                std::pow(pi_,  I_r(i)) *
                J_r(i) *
                std::pow(tau_ - 0.5, J_r(i) - 1);
        }

        return result;
    }

    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized pressure for IAPWS region 2 (i.e. sub-critical
     *        steam) (dimensionless).
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

        // ideal gas part
        Scalar result = 1/pi_;

        // residual part
        for (int i = 0; i < 43; i++) {
            result +=
                n_r(i) *
                I_r(i) *
                std::pow(pi_, I_r(i) - 1) *
                std::pow(tau_ - 0.5, J_r(i));
        }

        return result;
    }

    /*!
     * \brief The partial derivative of the Gibbs free energy to the
     *        normalized pressure and to the normalized temperature
     *        for IAPWS region 2 (i.e. sub-critical steam)  (dimensionless).
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

        // ideal gas part
        Scalar result = 0;

        // residual part
        for (int i = 0; i < 43; i++) {
            result +=
                n_r(i) *
                I_r(i) *
                J_r(i) *
                std::pow(pi_, I_r(i) - 1) *
                std::pow(tau_ - 0.5, J_r(i) - 1);
        }

        return result;
    }

    /*!
     * \brief The second partial derivative of the Gibbs free energy
     *        to the normalized pressure for IAPWS region 2
     *        (i.e. sub-critical steam) (dimensionless).
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

        // ideal gas part
        Scalar result = -1/(pi_*pi_);

        // residual part
        for (int i = 0; i < 43; i++) {
            result +=
                n_r(i) *
                I_r(i) *
                (I_r(i) - 1) *
                std::pow(pi_, I_r(i) - 2) *
                std::pow(tau_ - 0.5, J_r(i));
        }

        return result;
    }


private:
    static Scalar n_g(int i)
    {
        static const Scalar n[9] = {
            -0.96927686500217e1, 0.10086655968018e2, -0.56087911283020e-2,
            0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1,
            -0.43839511319450e1, -0.28408632460772, 0.21268463753307e-1
        };
        return n[i];
    }

    static Scalar n_r(int i)
    {
        static const Scalar n[43] = {
            -0.17731742473213e-2, -0.17834862292358e-1, -0.45996013696365e-1,
            -0.57581259083432e-1, -0.50325278727930e-1, -0.33032641670203e-4,
            -0.18948987516315e-3, -0.39392777243355e-2, -0.43797295650573e-1,
            -0.26674547914087e-4, 0.20481737692309e-7, 0.43870667284435e-6,
            -0.32277677238570e-4, -0.15033924542148e-2, -0.40668253562649e-1,
            -0.78847309559367e-9, 0.12790717852285e-7, 0.48225372718507e-6,
             0.22922076337661e-5, -0.16714766451061e-10, -0.21171472321355e-2,
            -0.23895741934104e2, -0.59059564324270e-17, -0.12621808899101e-5,
            -0.38946842435739e-1, 0.11256211360459e-10, -0.82311340897998e1,
             0.19809712802088e-7, 0.10406965210174e-18, -0.10234747095929e-12,
            -0.10018179379511e-8, -0.80882908646985e-10, 0.10693031879409,
            -0.33662250574171, 0.89185845355421e-24, 0.30629316876232e-12,
            -0.42002467698208e-5, -0.59056029685639e-25, 0.37826947613457e-5,
            -0.12768608934681e-14, 0.73087610595061e-28, 0.55414715350778e-16,
            -0.94369707241210e-6
        };
        return n[i];
    }

    static Scalar I_r(int i)
    {
        static const short int I[43] = {
            1, 1, 1,
            1, 1, 2,
            2, 2, 2,
            2, 3, 3,
            3, 3, 3,
            4, 4, 4,
            5, 6, 6,
            6, 7, 7,
            7, 8, 8,
            9, 10, 10,
            10, 16, 16,
            18, 20, 20,
            20, 21, 22,
            23, 24, 24,
            24
        };
        return I[i];
    }

    static Scalar J_g(int i)
    {
        static const short int J[9] = {
            0, 1, -5,
            -4, -3, -2,
            -1, 2, 3
        };
        return J[i];
    }

    static Scalar J_r(int i)
    {
        static const short int J[43] = {
            0, 1, 2,
            3, 6, 1,
            2, 4, 7,
            36, 0, 1,
            3, 6, 35,
            1, 2, 3,
            7, 3, 16,
            35, 0, 11,
            25, 8, 36,
            13, 4, 10,
            14, 29, 50,
            57, 20, 35,
            48, 21, 53,
            39, 26, 40,
            58
        };
        return J[i];
    }

};

} // end namepace IAPWS
} // end namepace Dune

#endif
