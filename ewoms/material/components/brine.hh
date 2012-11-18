// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::Brine
 */
#ifndef EWOMS_BRINE_HH
#define EWOMS_BRINE_HH

#include <ewoms/material/components/component.hh>

#include <cmath>

namespace Ewoms {

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
    static const char *name()
    { return "Brine"; }

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
    static Scalar vaporPressure(Scalar T)
    { return H2O::vaporPressure(T); /* [N/m^2] */ }

    /*!
     * \copydoc Component::gasEnthalpy
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    { return H2O::gasEnthalpy(temperature, pressure); /* [J/kg] */ }

    /*!
     * \copydoc Component::liquidEnthalpy
     *
     * Equations given in:
     * - Palliser & McKibbin 1997
     * - Michaelides 1981
     * - Daubert & Danner 1989
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        /*Numerical coefficents from PALLISER*/
        static const Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = {
            { -9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        Scalar theta, h_NaCl;
        Scalar m, h_ls, h_ls1, d_h;
        Scalar S_lSAT, delta_h;
        int i, j;
        Scalar hw;

        theta = temperature - 273.15;

        Scalar S = salinity;
        S_lSAT = f[0] + f[1]*theta + f[2]*std::pow(theta,2) + f[3]*std::pow(theta,3);
        /*Regularization*/
        if (S>S_lSAT) {
            S = S_lSAT;
        }

        hw = H2O::liquidEnthalpy(temperature, pressure)/1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*temperature + 0.5*(6.2770E1)*temperature*temperature - ((6.6670E-2)/3)*temperature*temperature*temperature
                        +((2.8000E-5)/4)*std::pow(temperature,4))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(S/(1-S));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * std::pow(theta, i) * std::pow(m, j);
            }
        }

        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine */

        h_ls1 =(1-S)*hw + S*h_NaCl + S*delta_h; /* kJ/kg */

        h_ls = h_ls1*1E3; /*J/kg*/

        return (h_ls);
    }


    /*!
     * \copydoc H2O::liquidHeatCapacity
     */
    static const Scalar liquidHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    {
        Scalar eps = temperature*1e-8;
        return (liquidEnthalpy(temperature + eps, pressure) - liquidEnthalpy(temperature, pressure))/eps;
    }

    /*!
     * \copydoc H2O::gasHeatCapacity
     */
    static const Scalar gasHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    { return H2O::gasHeatCapacity(temperature, pressure); }

    /*!
     * \copydoc H2O::gasInternalEnergy
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            pressure/gasDensity(temperature, pressure);
    }

    /*!
     * \copydoc H2O::liquidInternalEnergy
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }


    /*!
     * \copydoc H2O::gasDensity
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { return H2O::gasDensity(temperature, pressure); }

    /*!
     * \copydoc H2O::gasIsIdeal
     */
    static bool gasIsIdeal()
    { return H2O::gasIsIdeal(); }

    /*!
     * \copydoc Component::liquidDensity
     *
     * Equations given in:
     * - Batzle & Wang (1992)
     * - cited by: Adams & Bachu in Geofluids (2002) 2, 257-271
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        Scalar TempC = temperature - 273.15;
        Scalar pMPa = pressure/1.0E6;

        Scalar rhow = H2O::liquidDensity(temperature, pressure);
        return
            rhow +
            1000*salinity*(
                0.668 +
                0.44*salinity +
                1.0E-6*(
                    300*pMPa -
                    2400*pMPa*salinity +
                    TempC*(
                        80.0 -
                        3*TempC -
                        3300*salinity -
                        13*pMPa +
                        47*pMPa*salinity)));
    }

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
     * \copydoc H2O::gasPressure
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    { return H2O::gasPressure(temperature, density); }

    /*!
     * \copydoc H2O::liquidPressure
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        // We use the newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Scalar pressure = 1.1*vaporPressure(temperature);
        Scalar eps = pressure*1e-7;

        Scalar deltaP = pressure*2;
        for (int i = 0; i < 5 && std::abs(pressure*1e-9) < std::abs(deltaP); ++i) {
            Scalar f = liquidDensity(temperature, pressure) - density;

            Scalar df_dp;
            df_dp = liquidDensity(temperature, pressure + eps);
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
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { return H2O::gasViscosity(temperature, pressure); }

    /*!
     * \copydoc H2O::liquidViscosity
     *
     * Equation given in:
     * - Batzle & Wang (1992)
     * - cited by: Bachu & Adams (2002)
     *   "Equations of State for basin geofluids"
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        if(temperature <= 275.) // regularisation
        { temperature = 275; }
        Scalar T_C = temperature - 273.15;

        Scalar A = (0.42*std::pow((std::pow(salinity, 0.8)-0.17), 2) + 0.045)*std::pow(T_C, 0.8);
        Scalar mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity)*std::exp(-A);

        return mu_brine/1000.0; /* unit: Pa s */
    }
};

/*!
 * \brief Default value for the salinity of the brine (dimensionless).
 */
template <class Scalar, class H2O>
Scalar Brine<Scalar, H2O>::salinity = 0.1; // also needs to be adapted in CO2 solubility table!

} // end namepace

#endif
