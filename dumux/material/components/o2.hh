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
 * \brief Properties of pure molecular oxygen \f$O_2\f$.
 */
#ifndef DUMUX_O2_HH
#define DUMUX_O2_HH

#include <dumux/material/idealgas.hh>

#include "component.hh"

#include <cmath>

namespace Dumux
{

/*!
 * \brief Properties of pure molecular oxygen \f$O_2\f$.
 */
template <class Scalar>
class O2 : public Component<Scalar, O2<Scalar> >
{
    typedef Component<Scalar, O2<Scalar> >  ParentType;
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the O2.
     */
    static const char *name()
    { return "O2"; }

    /*!
     * \brief The mass in [kg] of one mole of molecular oxygen.
     */
    static Scalar molarMass()
    { return 32e-3; }

    /*!
     * \brief Returns the critical temperature of molecular oxygen
     */
    static Scalar criticalTemperature()
    { return 154.581; /* [K] */ }

    /*!
     * \brief Returns the critical pressure of molecular oxygen
     */
    static Scalar criticalPressure()
    { return 5.0804e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature at molecular oxygen's triple point.
     */
    static Scalar tripleTemperature()
    { return 54.359; /* [K] */ }

    /*!
     * \brief Returns the pressure at molecular oxygen's triple point.
     */
    static Scalar triplePressure()
    { return 148.0; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [N/m^2] of pure molecular oxygen
     *        at a given temperature.
     *
     * Taken from:
     *
     * R. Prydz: "An Improved Oxygen Vapor Pressure Representation",
     * Metrologia, Vol. 8, No. 1, pp. 1-4, 1972
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // O2 is solid: We don't take sublimation into account

        // vapor pressure between tripe and critical points.  See the
        // paper of Prydz for a discussion
        Scalar X =
            (1 - tripleTemperature()/T) /
            (1 - tripleTemperature()/criticalTemperature());
        const Scalar A = 7.568956;
        const Scalar B = 5.004836;
        const Scalar C = -2.137460;
        const Scalar D = 3.454481;
        const Scalar epsilon = 1.514;

        return
            triplePressure()*
            std::exp(X*(A +
                        X*(B + C*X) +
                        D*std::pow(1 - X,
                                   epsilon)));
    }

    /*!
     * \brief The density of pure O2 at a given pressure and temperature.
     *
     * \todo density liquid oxygen
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*
     * \brief The pressure of gaseous N2 at a given density and temperature [Pa].
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Specific enthalpy [J/kg] of pure oxygen gas.
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 154, 657, 665
     */
    static const Scalar gasEnthalpy(Scalar T,
                                    Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA = 28.11;
        const Scalar cpVapB = -3.680e-6;
        const Scalar cpVapC = 1.746e-5;
        const Scalar cpVapD = -1.065e-8;

        //Scalar cp =
        //    cpVapA + T*(cpVapB + T*(cpVapC + T*cpVapD));

        // calculate: \int_0^T c_p dT
        return
            1/molarMass()* // conversion from [J/mol] to [J/kg]

            T*(cpVapA + T*
               (cpVapB/2 + T*
                (cpVapC/3 + T*
                 (cpVapD/4))));
    }

    /*!
     * \brief The density [kg/m^3] of gaseous O2 at a given pressure and temperature.
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidDensity for O2"); }

    /*
     * \brief The pressure of liquid oxygen at a given density and
     *        temperature [Pa].
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    { DUNE_THROW(Dune::NotImplemented, "liquidPressure for O2"); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure liquid O2.
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for O2"); }

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
        const Scalar Vc = 73.4; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.025; // accentric factor
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
     * \brief The dynamic liquid viscosity [N/m^3*s] of pure H2.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidViscosity for O2"); }

};

} // end namepace

#endif
