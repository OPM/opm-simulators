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
 * \copydoc Opm::IAPWS::Region1
 */
#ifndef OPM_IAPWS_REGION1_HPP
#define OPM_IAPWS_REGION1_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
namespace IAPWS {
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
    template <class Evaluation>
    static bool isValid(const Evaluation& temperature, const Evaluation& pressure)
    {
        return temperature <= 623.15 && pressure <= 100e6;

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
    template <class Evaluation>
    static Evaluation tau(const Evaluation& temperature)
    { return 1386.0 / temperature; }

    /*!
     * \brief Returns the derivative of the reduced temperature to the
     *        temperature for IAPWS region 1 in \f$\mathrm{[1/K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    static Evaluation dtau_dT(const Evaluation& temperature)
    { return - 1386.0 / (temperature*temperature); }

    /*!
     * \brief Returns the reduced pressure for IAPWS region 1.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation pi(const Evaluation& pressure)
    { return pressure / 16.53e6; }

    /*!
     * \brief Returns the derivative of the reduced pressure to the
     *        pressure for IAPWS region 1 in \f$\mathrm{[1/Pa]}\f$.
     *
     * \param pressure temperature of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Scalar dpi_dp(const Evaluation& /*pressure*/)
    { return 1.0 / 16.53e6; }

    /*!
     * \brief Returns the derivative of the pressure to the
     *        reduced pressure for IAPWS region 1 in \f$\mathrm{[Pa]}\f$.
     *
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Scalar dp_dpi(const Evaluation& /*pressure*/)
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
    template <class Evaluation>
    static Evaluation gamma(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0;
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
    template <class Evaluation>
    static Evaluation dgamma_dtau(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                n(i) *
                pow(7.1 - pi_, static_cast<Scalar>(I(i))) *
                pow(tau_ - 1.222,  static_cast<Scalar>(J(i)-1)) *
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
    template <class Evaluation>
    static Evaluation dgamma_dpi(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                -n(i) *
                I(i) *
                pow(7.1 - pi_, static_cast<Scalar>(I(i) - 1)) *
                pow(tau_ - 1.222, static_cast<Scalar>(J(i)));
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
    template <class Evaluation>
    static Evaluation ddgamma_dtaudpi(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                -n(i) *
                I(i) *
                J(i) *
                pow(7.1 - pi_, static_cast<Scalar>(I(i) - 1)) *
                pow(tau_ - 1.222, static_cast<Scalar>(J(i) - 1));
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
    template <class Evaluation>
    static Evaluation ddgamma_ddpi(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                n(i) *
                I(i) *
                (I(i) - 1) *
                pow(7.1 - pi_, I(i) - 2) *
                pow(tau_ - 1.222, J(i));
        }

        return result;
    }

    /*!
     * \brief The second partial derivative of the Gibbs free energy to the
     *        normalized temperature for IAPWS region 1 (i.e. liquid) (dimensionless).
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    template <class Evaluation>
    static Evaluation ddgamma_ddtau(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation tau_ = tau(temperature);   /* reduced temperature */
        const Evaluation pi_ = pi(pressure);    /* reduced pressure */

        Evaluation result = 0.0;
        for (int i = 0; i < 34; i++) {
            result +=
                n(i) *
                pow(7.1 - pi_, I(i)) *
                J(i) *
                (J(i) - 1) *
                pow(tau_ - 1.222,  J(i) - 2);
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

    static Scalar I(int i)
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

    static Scalar J(int i)
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

} // namespace IAPWS
} // namespace Opm

#endif
