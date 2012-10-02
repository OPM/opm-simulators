// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::CO2
 */
#ifndef DUMUX_CO2_HH
#define DUMUX_CO2_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/idealgas.hh>

#include <cmath>
#include <iostream>

namespace Dumux {

/*!
 * \brief A class for the CO2 fluid properties
 *
 * Under reservoir conditions, CO2 is typically in supercritical
 * state. These properties can be provided in tabulated form, which is
 * necessary for this component. The template is used by the
 * fluidsystem \c FluidSystems::BrineCO2. If thermodynamic precision
 * is not a top priority, the much simpler component \c Dumux::SimpleCO2 can be
 * used instead
 */
template <class Scalar, class CO2Tables>
class CO2 : public Component<Scalar, CO2<Scalar, CO2Tables> >
{
    static constexpr Scalar R = Constants<Scalar>::R;
    typedef typename Dumux::IdealGas<Scalar> IdealGas;
    
    static bool warningPrinted;

public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static const char *name()
    { return "CO2"; }

    /*!
     * \brief The mass in [kg] of one mole of CO2.
     */
    static Scalar molarMass()
    { return 44e-3; }

    /*!
     * \brief Returns the critical temperature [K] of CO2
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of CO2
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K]at CO2's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar minTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.minPress(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar maxTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.maxPress(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar minTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.minTemp(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar maxTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.maxTemp(); /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [N/m^2] of pure CO2
     *        at a given temperature.
     *
     * See:
     *
     * R. Span and W. Wagner: A New Equation of State for Carbon
     * Dioxide Covering the Fluid Region from the Triple‐Point
     * Temperature to 1100 K at Pressures up to 800 MPa. Journal of
     * Physical and Chemical Reference Data, 25 (6), pp. 1509-1596,
     * 1996
     */
    static Scalar vaporPressure(Scalar T)
    { 
        static const Scalar a[4] = 
            { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
        static const Scalar t[4] = 
            { 1.0, 1.5, 2.0, 4.0 };

        // this is on page 1524 of the reference
        Scalar exponent = 0;
        Scalar Tred = T/criticalTemperature();
        for (int i = 0; i < 5; ++i) {
            exponent += a[i]*std::pow(1 - Tred, t[i]);
        }
        exponent *= 1.0/Tred;
        
        return std::exp(exponent)*criticalPressure();
    }


    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return false; }

    /*!
     * \brief Specific enthalpy of gaseous CO2 [J/kg].
     */
    static Scalar gasEnthalpy(Scalar temperature,
                              Scalar pressure)
    {
#ifndef NDEBUG
        if ((temperature < criticalTemperature() or pressure < criticalPressure()) and !warningPrinted)
        {
            Dune::dwarn << "The CO2 exhibits subcritical values: Be aware "
                        << "to use tables with sufficient resolution!\n";
            warningPrinted=true;
        }
#endif

        return CO2Tables::tabulatedEnthalpy.eval(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of CO2 [J/kg].
     */
    static Scalar gasInternalEnergy(Scalar temperature,
                                    Scalar pressure)
    {
        Scalar h = gasEnthalpy(temperature, pressure);
        Scalar rho = gasDensity(temperature, pressure);

        return h - (pressure / rho);
    }

    /*!
     * \brief The density of CO2 at a given pressure and temperature [kg/m^3].
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
#ifndef NDEBUG
        if ((temperature < criticalTemperature() or pressure < criticalPressure()) and !warningPrinted)
        {
            Dune::dwarn << "Subcritical values: Be aware to use "
                        <<"Tables with sufficient resolution!"<< std::endl;
            warningPrinted=true;
        }
#endif

        return CO2Tables::tabulatedDensity.eval(temperature, pressure);
    }

    /*!
     * \brief The dynamic viscosity [Pa s] of CO2.
     *
     * Equations given in: - Vesovic et al., 1990
     *                        - Fenhour etl al., 1998
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const double a0 = 0.235156;
        static const double a1 = -0.491266;
        static const double a2 = 5.211155E-2;
        static const double a3 = 5.347906E-2;
        static const double a4 = -1.537102E-2;

        static const double d11 = 0.4071119E-2;
        static const double d21 = 0.7198037E-4;
        static const double d64 = 0.2411697E-16;
        static const double d81 = 0.2971072E-22;
        static const double d82 = -0.1627888E-22;

        static const double ESP = 251.196;

        double mu0, SigmaStar, TStar;
        double dmu, rho;
        double visco_CO2;

        if(temperature < 275.) // regularisation
        {
            temperature = 275;
            Dune::dgrave << "Temperature below 275K in viscosity function:"
                    << "Regularizing tempereature to 275K. " << std::endl;
        }


        TStar = temperature/ESP;

        /* mu0: viscosity in zero-density limit */
        SigmaStar = exp(a0 + a1*log(TStar)
                        + a2*log(TStar)*log(TStar)
                        + a3*log(TStar)*log(TStar)*log(TStar)
                        + a4*log(TStar)*log(TStar)*log(TStar)*log(TStar) );

        mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        /* dmu : excess viscosity at elevated density */
        rho = gasDensity(temperature, pressure); /* CO2 mass density [kg/m^3] */

        dmu = d11*rho + d21*rho*rho + d64*pow(rho,6)/(TStar*TStar*TStar)
            + d81*pow(rho,8) + d82*pow(rho,8)/TStar;

        /* dmucrit : viscosity increase near the critical point */

        // False (Lybke 2July2007)
        //e1 = 5.5930E-3;
        //e2 = 6.1757E-5;
        //e4 = 2.6430E-11;
        //dmucrit = e1*rho + e2*rho*rho + e4*rho*rho*rho;
        //visco_CO2 = (mu0 + dmu + dmucrit)/1.0E6;   /* conversion to [Pa s] */

        visco_CO2 = (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */

        return visco_CO2;
    };
};

template <class Scalar, class CO2Tables>
bool CO2<Scalar, CO2Tables>::warningPrinted = false;

} // end namepace

#endif
