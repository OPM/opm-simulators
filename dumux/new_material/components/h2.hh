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
 * \brief Properties of pure molecular nitrogen \f$H_2\f$.
 */
#ifndef DUMUX_H2_HH
#define DUMUX_H2_HH

#include <dumux/new_material/idealgas.hh>
#include <dune/common/exceptions.hh>

#include "component.hh"

#include <cmath>

namespace Dumux
{

/*!
 * \brief Properties of pure molecular hydrogen \f$H_2\f$.
 */
template <class Scalar>
class H2 : public Component<Scalar, H2<Scalar> >
{
    typedef Component<Scalar, H2<Scalar> >  ParentType;
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the H2.
     */
    static const char *name()
    { return "H2"; }

    /*!
     * \brief The mass in [kg/mol] of one of molecular hydrogen.
     */
    static Scalar molarMass()
    { return 1.0e-3; }

    /*!
     * \brief Returns the critical temperature [K] of molecular hydrogen
     */
    static Scalar criticalTemperature()
    { return 33.2; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of molecular hydrogen
     */
    static Scalar criticalPressure()
    { return 13.0e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K] at molecular hydrogen's triple point.
     */
    static Scalar tripleTemperature()
    { return 14.0; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at molecular hydrogen's triple point.
     */
    static Scalar triplePressure()
    { DUNE_THROW(Dune::NotImplemented, "triplePressure for H2"); }

    /*!
     * \brief The vapor pressure in [Pa] of pure molecular hydrogen
     *        at a given temperature.
     *
     * Taken from:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 208-209, 669
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        if (temperature > criticalTemperature())
            return criticalPressure();
        if (temperature < tripleTemperature())
            return 0; // H2 is solid: We don't take sublimation into
                      // account

        // TODO: the Gomez-Thodos approach would probably be better...

        // antoine equatuion
        const Scalar A = -7.76451;
        const Scalar B = 1.45838;
        const Scalar C = -2.77580;

        return 1e5 * std::exp(A - B/(temperature + C));
    }

    /*!
     * \brief The density [kg/m^3] of H2 at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
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
     * \brief Specific enthalpy [J/kg] of pure hydrogen gas.
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 154, 657, 665
     */
    static const Scalar gasEnthalpy(Scalar T,
                                    Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA =  27.14;
        const Scalar cpVapB =  9.273e-3;
        const Scalar cpVapC = -1.381e-5;
        const Scalar cpVapD =  7.645e-9;

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
     * \brief The density [kg/m^3] of liquid hydrogen at a given pressure and temperature.
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidDensity for H2"); }

    /*
     * \brief The pressure of liquid hydrogen at a given density and
     *        temperature [Pa].
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    { DUNE_THROW(Dune::NotImplemented, "liquidPressure for H2"); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure liquid H2 .
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for H2"); }

    /*!
     * \brief The dynamic viscosity [Pa s] of H2 at a given pressure and temperature.
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 396-397, 667
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 65.1; // critical specific volume [cm^3/mol]
        const Scalar omega = -0.218; // accentric factor
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
    { DUNE_THROW(Dune::NotImplemented, "liquidViscosity for H2"); }
};

} // end namepace

#endif
