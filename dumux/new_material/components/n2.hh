/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Properties of pure molecular nitrogen \f$N_2\f$.
 */
#ifndef DUMUX_N2_HH
#define DUMUX_N2_HH

#include <dumux/new_material/idealgas.hh>
#include <dune/common/exceptions.hh>

#include "component.hh"

#include <cmath>

namespace Dune
{

/*!
 * \brief Properties of pure molecular nitrogen \f$N_2\f$.
 *
 * See: 
 *
 * R. Span, E.W. Lemmon, et al.: "A Reference Equation of State for
 * the Thermodynamic Properties of Nitrogen for Temperatures from
 * 63.151 to 1000 K and Pressures to 2200 MPa", Journal of Physical
 * and Chemical Refefence Data, Vol. 29, No. 6, pp. 1361-1433
 */
template <class Scalar>
class N2 : public Component<Scalar, N2<Scalar> >
{
    typedef Component<Scalar, N2<Scalar> >  ParentType;
    typedef Dune::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for nitrogen.
     */
    static const char *name()
    { return "N2"; } 

    /*!
     * \brief The mass in [kg/mol] of one of molecular nitrogen.
     */
    static Scalar molarMass()
    { return 28.0134e-3;} 

    /*!
     * \brief Returns the critical temperature [K] of molecular nitrogen
     */
    static Scalar criticalTemperature()
    { return 126.192; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of molecular nitrogen
     */
    static Scalar criticalPressure()
    { return 3.39858e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K] at molecular nitrogen's triple point.
     */
    static Scalar tripleTemperature()
    { return 63.151; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at molecular nitrogen's triple point.
     */
    static Scalar triplePressure()
    { return 12.523e3; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [Pa] of pure molecular nitrogen
     *        at a given temperature.
     *
     * Taken from: 
     *
     * R. Span, E.W. Lemmon, et al.: "A Reference Equation of State
     * for the Thermodynamic Properties of Nitrogen for Temperatures
     * from 63.151 to 1000 K and Pressures to 2200 MPa", Journal of
     * Physical and Chemical Refefence Data, Vol. 29, No. 6,
     * pp. 1361-1433
     */
    static Scalar vaporPressure(Scalar T)
    { 
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // N2 is solid: We don't take sublimation into
                      // account

        // note: this is the ancillary equation given on page 1368
        Scalar sigma = 1.0 - T/criticalTemperature();
        Scalar sqrtSigma = sqrt(sigma);
        const Scalar N1 = -6.12445284;
        const Scalar N2 = 1.26327220;
        const Scalar N3 = -0.765910082;
        const Scalar N4 = -1.77570564;
        return
            criticalPressure() *
            exp(criticalTemperature()/T*
                (sigma*(N1 +
                        sqrtSigma*N2 + 
                        sigma*(sqrtSigma*N3 + 
                               sigma*sigma*sigma*N4))));
    }

    /*!
     * \brief The density [kg/m^3] of N2 at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), pressure, temperature);
    }
    
    /*!
     * \brief Specific enthalpy [J/kg] of pure nitrogen gas.
     */
    static const Scalar gasEnthalpy(Scalar temperature, 
                                    Scalar pressure)
    {
        const Scalar cvHat = 5./2; // constant for diatomic gases
        
        return temperature * (cvHat + 1) * IdealGas::R / molarMass();
    }

    /*!
     * \brief Specific enthalpy [J/kg] of pure liquid N2.
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(NotImplemented, "liquidEnthalpy for N2"); }

    /*!
     * \brief The dynamic viscosity [Pa s] of N2 at a given pressure and temperature.
     *
     * See: 
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 396-397, 664
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 89.8; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.0039; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]
        
        Scalar mu_r4 = 131.3 * dipole / std::sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;
        
        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Scalar Tstar = 1.2593 * temperature/Tc;
        Scalar Omega_v =
            1.16145*std::pow(Tstar, -0.14874) +
            0.52487*std::exp(- 0.77320*Tstar) +
            2.16178*std::exp(- 2.43787*Tstar);
        Scalar mu = 40.785*Fc*std::sqrt(M*temperature)/(std::pow(Vc, 2./3)*Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }

    /*!
     * \brief The dynamic liquid viscosity [N/m^3*s] of pure N2.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(NotImplemented, "liquidViscosity for N2"); }
};

} // end namepace

#endif
