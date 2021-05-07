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
 * \copydoc Opm::Brine
 */
#ifndef OPM_BRINE_HPP
#define OPM_BRINE_HPP

#include <opm/material/components/Component.hpp>
#include <opm/material/common/MathToolbox.hpp>

namespace Opm {

/*!
 * \ingroup Components
 *
 * \brief A class for the brine fluid properties.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam H2O Static polymorphism: the Brine class can access all properties of the H2O class
 */
template <class Scalar, class H2O>
class Brine : public Component<Scalar, Brine<Scalar, H2O> >
{
public:
    //! The mass fraction of salt assumed to be in the brine.
    static Scalar salinity;

    /*!
     * \copydoc Component::name
     */
    static const char* name()
    { return "Brine"; }

    /*!
     * \copydoc H2O::gasIsIdeal
     */
    static bool gasIsIdeal()
    { return H2O::gasIsIdeal(); }

    /*!
     * \copydoc H2O::gasIsCompressible
     */
    static bool gasIsCompressible()
    { return H2O::gasIsCompressible(); }

    /*!
     * \copydoc H2O::liquidIsCompressible
     */
    static bool liquidIsCompressible()
    { return H2O::liquidIsCompressible(); }

    /*!
     * \copydoc Component::molarMass
     *
     * This assumes that the salt is pure NaCl.
     */
    static Scalar molarMass()
    {
        const Scalar M1 = H2O::molarMass();
        const Scalar M2 = 58e-3; // molar mass of NaCl [kg/mol]
        const Scalar X2 = salinity; // mass fraction of salt in brine
        return M1*M2/(M2 + X2*(M1 - M2));
    }

    /*!
     * \copydoc H2O::criticalTemperature
     */
    static Scalar criticalTemperature()
    { return H2O::criticalTemperature(); /* [K] */ }

    /*!
     * \copydoc H2O::criticalPressure
     */
    static Scalar criticalPressure()
    { return H2O::criticalPressure(); /* [N/m^2] */ }

    /*!
     * \copydoc H2O::tripleTemperature
     */
    static Scalar tripleTemperature()
    { return H2O::tripleTemperature(); /* [K] */ }

    /*!
     * \copydoc H2O::triplePressure
     */
    static Scalar triplePressure()
    { return H2O::triplePressure(); /* [N/m^2] */ }

    /*!
     * \copydoc H2O::vaporPressure
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& T)
    { return H2O::vaporPressure(T); /* [N/m^2] */ }

    /*!
     * \copydoc Component::gasEnthalpy
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    { return H2O::gasEnthalpy(temperature, pressure); /* [J/kg] */ }

    /*!
     * \copydoc Component::liquidEnthalpy
     *
     * Equations given in:
     * - Palliser & McKibbin 1997
     * - Michaelides 1981
     * - Daubert & Danner 1989
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure)
    {
        // Numerical coefficents from Palliser and McKibbin
        static const Scalar f[] = {
            2.63500e-1, 7.48368e-6, 1.44611e-6, -3.80860e-10
        };

        // Numerical coefficents from Michaelides for the enthalpy of brine
        static const Scalar a[4][3] = {
            { -9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667e-1 },
            { +0.17965e-2, +0.71924e-3, -0.4900e-4 }
        };

        const Evaluation& theta = temperature - 273.15;

        Evaluation S = salinity;
        const Evaluation& S_lSAT =
            f[0]
            + f[1]*theta
            + f[2]*pow(theta, 2)
            + f[3]*pow(theta, 3);

        // Regularization
        if (S > S_lSAT)
            S = S_lSAT;

        const Evaluation& hw = H2O::liquidEnthalpy(temperature, pressure)/1e3; // [kJ/kg]

        // From Daubert and Danner
        const Evaluation& h_NaCl =
            (3.6710e4*temperature
             + (6.2770e1/2)*temperature*temperature
             - (6.6670e-2/3)*temperature*temperature*temperature
             + (2.8000e-5/4)*pow(temperature, 4.0))/58.44e3
            - 2.045698e+02; // [kJ/kg]

        const Evaluation& m = S/(1-S)/58.44e-3;

        Evaluation d_h = 0;
        for (int i = 0; i<=3; ++i) {
            for (int j = 0; j <= 2; ++j) {
                d_h += a[i][j] * pow(theta, i) * pow(m, j);
            }
        }

        const Evaluation& delta_h = 4.184/(1e3 + (58.44 * m))*d_h;

        // Enthalpy of brine
        const Evaluation& h_ls = (1-S)*hw + S*h_NaCl + S*delta_h; // [kJ/kg]
        return h_ls*1e3; // convert to [J/kg]
    }


    /*!
     * \copydoc H2O::liquidHeatCapacity
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature,
                                         const Evaluation& pressure)
    {
        Scalar eps = scalarValue(temperature)*1e-8;
        return (liquidEnthalpy(temperature + eps, pressure) - liquidEnthalpy(temperature, pressure))/eps;
    }

    /*!
     * \copydoc H2O::gasHeatCapacity
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature,
                                      const Evaluation& pressure)
    { return H2O::gasHeatCapacity(temperature, pressure); }

    /*!
     * \copydoc H2O::gasInternalEnergy
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            pressure/gasDensity(temperature, pressure);
    }

    /*!
     * \copydoc H2O::liquidInternalEnergy
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& temperature,
                                           const Evaluation& pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \copydoc H2O::gasDensity
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    { return H2O::gasDensity(temperature, pressure); }

    /*!
     * \copydoc Component::liquidDensity
     *
     * Equations given in:
     * - Batzle & Wang (1992)
     * - cited by: Adams & Bachu in Geofluids (2002) 2, 257-271
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        Evaluation tempC = temperature - 273.15;
        Evaluation pMPa = pressure/1.0E6;

        const Evaluation rhow = H2O::liquidDensity(temperature, pressure);
        return
            rhow +
            1000*salinity*(
                0.668 +
                0.44*salinity +
                1.0E-6*(
                    300*pMPa -
                    2400*pMPa*salinity +
                    tempC*(
                        80.0 -
                        3*tempC -
                        3300*salinity -
                        13*pMPa +
                        47*pMPa*salinity)));
    }

    /*!
     * \copydoc H2O::gasPressure
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    { return H2O::gasPressure(temperature, density); }

    /*!
     * \copydoc H2O::liquidPressure
     */
    template <class Evaluation>
    static Evaluation liquidPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // We use the newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Evaluation pressure = 1.1*vaporPressure(temperature);
        Scalar eps = scalarValue(pressure)*1e-7;

        Evaluation deltaP = pressure*2;
        for (int i = 0;
             i < 5
                 && std::abs(scalarValue(pressure)*1e-9) < std::abs(scalarValue(deltaP));
             ++i)
        {
            const Evaluation& f = liquidDensity(temperature, pressure) - density;

            Evaluation df_dp = liquidDensity(temperature, pressure + eps);
            df_dp -= liquidDensity(temperature, pressure - eps);
            df_dp /= 2*eps;

            deltaP = - f/df_dp;

            pressure += deltaP;
        }

        return pressure;
    }

    /*!
     * \copydoc H2O::gasViscosity
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& pressure)
    { return H2O::gasViscosity(temperature, pressure); }

    /*!
     * \copydoc H2O::liquidViscosity
     *
     * Equation given in:
     * - Batzle & Wang (1992)
     * - cited by: Bachu & Adams (2002)
     *   "Equations of State for basin geofluids"
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        Evaluation T_C = temperature - 273.15;
        if(temperature <= 275.) // regularization
            T_C = 275.0;

        Evaluation A = (0.42*std::pow((std::pow(salinity, 0.8)-0.17), 2) + 0.045)*pow(T_C, 0.8);
        Evaluation mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity)*exp(-A);

        return mu_brine/1000.0; // convert to [Pa s] (todo: check if correct cP->Pa s is times 10...)
    }
};

/*!
 * \brief Default value for the salinity of the brine (dimensionless).
 */
template <class Scalar, class H2O>
Scalar Brine<Scalar, H2O>::salinity = 0.1; // also needs to be adapted in CO2 solubility table!

} // namespace Opm

#endif
