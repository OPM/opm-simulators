// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
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
 * \ingroup Components
 *
 * \brief Material properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_H2O_HH
#define DUMUX_H2O_HH

#include <dumux/material/idealgas.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

#include "component.hh"

#include "iapws/common.hh"
#include "iapws/region1.hh"
#include "iapws/region2.hh"
#include "iapws/region4.hh"

#include <cmath>


namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief Material properties of pure water \f$H_2O\f$.
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
class H2O : public Component<Scalar, H2O<Scalar> >
{
    typedef Component<Scalar, H2O<Scalar> > ParentType;

    typedef IAPWS::Common<Scalar> Common;
    typedef IAPWS::Region1<Scalar> Region1;
    typedef IAPWS::Region2<Scalar> Region2;
    typedef IAPWS::Region4<Scalar> Region4;

    static const Scalar R = Common::R;  // specific gas constant of water
public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "H2O"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static Scalar molarMass()
    { return Common::molarMass; }

    /*!
     * \brief The acentric factor \f$\mathrm{[-]}\f$ of water.
     */
    static Scalar acentricFactor()
    { return Common::acentricFactor; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water
     */
    static Scalar criticalTemperature()
    { return Common::criticalTemperature; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static Scalar criticalPressure()
    { return Common::criticalPressure; }

    /*!
     * \brief Returns the molar volume \f$\mathrm{[m^3/mol]}\f$ of water at the critical point
     */
    static Scalar criticalMolarVolume()
    { return Common::criticalMolarVolume; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static Scalar tripleTemperature()
    { return Common::tripleTemperature; }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static Scalar triplePressure()
    { return Common::triplePressure; }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            T = criticalTemperature();
        if (T < tripleTemperature())
            T = tripleTemperature();

        return Region4::saturationPressure(T);
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the gas enthalpy at very low pressures. The
            // enthalpy of an ideal gas does not exhibit any
            // dependence on pressure, so we can just return the
            // specific enthalpy at the point of regularization, i.e.
            // the triple pressure - 100Pa
            return enthalpyRegion2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            Scalar dh_dp =
                R*temperature*
                Region2::tau(temperature)*
                Region2::dpi_dp(pv)*
                Region2::ddgamma_dtaudpi(temperature, pv);

            return
                enthalpyRegion2_(temperature, pv) +
                (pressure - pv)*dh_dp;
        };

        return enthalpyRegion2_(temperature, pressure);
    }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            Scalar dh_dp =
                R*temperature*
                Region1::tau(temperature)*
                Region1::dpi_dp(pv)*
                Region1::ddgamma_dtaudpi(temperature, pv);

            return
                enthalpyRegion1_(temperature, pv) +
                (pressure - pv)*dh_dp;
        };

        return enthalpyRegion1_(temperature, pressure);
    }

    /*!
     * \brief Specific isobaric heat capacity of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar gasHeatCap_p(Scalar temperature,
                                    Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Heat capacity of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            return heatCap_p_Region2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the heat cap at the vapor pressure to regularize
            return
                heatCap_p_Region2_(temperature, pv);
        };
        return heatCap_p_Region2_(temperature, pressure);
    }

    /*!
     * \brief Specific isobaric heat capacity of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidHeatCap_p(Scalar temperature,
                                       Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "heat Capacity of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the heat cap at the vapor pressure to regularize
            return
                heatCap_p_Region1_(temperature, pv);
        };

        return heatCap_p_Region1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }


        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the internal energy at the vapor pressure to
            // regularize

            /*
            // calculate the partial derivative of the internal energy
            // to the pressure at the vapor pressure.
            Scalar tau = Region1::tau(temperature);
            Scalar dgamma_dpi = Region1::dgamma_dpi(temperature, pv);
            Scalar ddgamma_dtaudpi = Region1::ddgamma_dtaudpi(temperature, pv);
            Scalar ddgamma_ddpi = Region1::ddgamma_ddpi(temperature, pv);
            Scalar pi = Region1::pi(pv);
            Scalar dpi_dp = Region1::dpi_dp(pv);
            Scalar du_dp =
                R*temperature*
                (tau*dpi_dp*ddgamma_dtaudpi + dpi_dp*dpi_dp*dgamma_dpi + pi*dpi_dp*ddgamma_ddpi);
            */

            // use a straight line for extrapolation. use forward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            static const Scalar eps = 1e-7;
            Scalar uv = internalEnergyRegion1_(temperature, pv);
            Scalar uvPEps = internalEnergyRegion1_(temperature, pv + eps);
            Scalar du_dp = (uvPEps - uv)/eps;
            return uv + du_dp*(pressure - pv);
        };

        return internalEnergyRegion1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of steam and water vapor \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
    */
    static Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the internal energy of gas at very low
            // pressures. The enthalpy of an ideal gas does not
            // exhibit any dependence on pressure, so we can just
            // return the specific enthalpy at the point of
            // regularization, i.e.  the triple pressure - 100Pa, and
            // subtract the work required to change the volume for an
            // ideal gas.
            return
                enthalpyRegion2_(temperature, triplePressure() - 100)
                -
                R*temperature; // = p*v   for an ideal gas!
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the internal energy at the vapor pressure to
            // regularize

            /*
            // calculate the partial derivative of the internal energy
            // to the pressure at the vapor pressure.
            Scalar tau = Region2::tau(temperature);
            Scalar dgamma_dpi = Region2::dgamma_dpi(temperature, pv);
            Scalar ddgamma_dtaudpi = Region2::ddgamma_dtaudpi(temperature, pv);
            Scalar ddgamma_ddpi = Region2::ddgamma_ddpi(temperature, pv);
            Scalar pi = Region2::pi(pv);
            Scalar dpi_dp = Region2::dpi_dp(pv);
            Scalar du_dp =
                R*temperature*
                (tau*dpi_dp*ddgamma_dtaudpi + dpi_dp*dpi_dp*dgamma_dpi + pi*dpi_dp*ddgamma_ddpi);

            // use a straight line for extrapolation
            Scalar uv = internalEnergyRegion2_(temperature, pv);
            return uv + du_dp*(pressure - pv);
            */

            // use a straight line for extrapolation. use backward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            static const Scalar eps = 1e-7;
            Scalar uv = internalEnergyRegion2_(temperature, pv);
            Scalar uvMEps = internalEnergyRegion2_(temperature, pv - eps);
            Scalar du_dp = (uv - uvMEps)/eps;
            return uv + du_dp*(pressure - pv);
        };

        return internalEnergyRegion2_(temperature, pressure);
    }

    /*!
     * \brief Specific isochoric heat capacity of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidHeatCap_v(Scalar temperature,
                                             Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Heat capacity of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }


        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the heat cap at the vapor pressure to regularize

            return heatCap_v_Region1_(temperature, pv);
        }

        return heatCap_v_Region1_(temperature, pressure);
    }

    /*!
     * \brief Specific isochoric heat capacity of steam and water vapor \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
    */
    static Scalar gasHeatCap_v(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Heat capacity of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            return
                heatCap_v_Region2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            return heatCap_v_Region2_(temperature, pv);
        };

        return heatCap_v_Region2_(temperature, pressure);
    }

    /*!
     * \brief The density of steam in \f$\mathrm{[kg/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Density of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the internal energy and enthalpy.
            Scalar rho0IAPWS = 1.0/volumeRegion2_(temperature,
                                                  triplePressure() - 100);
            Scalar rho0Id = IdealGas<Scalar>::density(molarMass(),
                                                      temperature,
                                                      triplePressure() - 100);
            return
                rho0IAPWS/rho0Id *
                IdealGas<Scalar>::density(molarMass(),
                                          temperature,
                                          pressure);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the density energy at the vapor pressure to
            // regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            const Scalar eps = pv*1e-8;
            Scalar v0 = volumeRegion2_(temperature, pv);
            Scalar v1 = volumeRegion2_(temperature, pv + eps);
            Scalar dv_dp = (v1 - v0)/eps;
            /*
            Scalar pi = Region2::pi(pv);
            Scalar dp_dpi = Region2::dp_dpi(pv);
            Scalar dgamma_dpi = Region2::dgamma_dpi(temperature, pv);
            Scalar ddgamma_ddpi = Region2::ddgamma_ddpi(temperature, pv);

            Scalar RT = R*temperature;
            Scalar dv_dp =
                RT/(dp_dpi*pv)
                *
                (dgamma_dpi + pi*ddgamma_ddpi - v0*dp_dpi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Scalar drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        };

        return 1.0/volumeRegion2_(temperature, pressure);
    }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(density);

        // We use the newton method for this. For the initial value we
        // assume steam to be an ideal gas
        Scalar pressure = IdealGas<Scalar>::pressure(temperature, density/molarMass());
        Scalar eps = pressure*1e-7;

        Scalar deltaP = pressure*2;
        Valgrind::CheckDefined(pressure);
        Valgrind::CheckDefined(deltaP);
        for (int i = 0; i < 5 && std::abs(pressure*1e-9) < std::abs(deltaP); ++i) {
            Scalar f = gasDensity(temperature, pressure) - density;

            Scalar df_dp;
            df_dp = gasDensity(temperature, pressure + eps);
            df_dp -= gasDensity(temperature, pressure - eps);
            df_dp /= 2*eps;

            deltaP = - f/df_dp;

            pressure += deltaP;
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(deltaP);
        }

        return pressure;
    }

    /*!
     * \brief The density of pure water in \f$\mathrm{[kg/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Density of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the density at the vapor pressure to regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            const Scalar eps = pv*1e-8;
            Scalar v0 = volumeRegion1_(temperature, pv);
            Scalar v1 = volumeRegion1_(temperature, pv + eps);
            Scalar dv_dp = (v1 - v0)/eps;

            /*
            Scalar v0 = volumeRegion1_(temperature, pv);
            Scalar pi = Region1::pi(pv);
            Scalar dp_dpi = Region1::dp_dpi(pv);
            Scalar dgamma_dpi = Region1::dgamma_dpi(temperature, pv);
            Scalar ddgamma_ddpi = Region1::ddgamma_ddpi(temperature, pv);

            Scalar RT = R*temperature;
            Scalar dv_dp =
                RT/(dp_dpi*pv)
                *
                (dgamma_dpi + pi*ddgamma_ddpi - v0*dp_dpi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Scalar drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        };

        return 1/volumeRegion1_(temperature, pressure);
    }

    /*!
     * \brief The pressure of liquid water in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * This method is only valid if pressure is below or at the vapor
     * pressure of water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar rho = gasDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    };

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        };

        Scalar rho = liquidDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    };

private:
    // the unregularized specific enthalpy for liquid water
    static Scalar enthalpyRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            Region1::tau(temperature) *
            Region1::dgamma_dtau(temperature, pressure) *
            R*temperature;
    };

    // the unregularized specific isobaric heat capacity
    static Scalar heatCap_p_Region1_(Scalar temperature, Scalar pressure)
    {
        return
            - pow(Region1::tau(temperature), 2 ) *
            Region1::ddgamma_ddtau(temperature, pressure) *
            R;
    };

    // the unregularized specific isochoric heat capacity
    static Scalar heatCap_v_Region1_(Scalar temperature, Scalar pressure)
    {
        double tau = Region1::tau(temperature);
        double num = Region1::dgamma_dpi(temperature, pressure) - tau * Region1::ddgamma_dtaudpi(temperature, pressure);
        double diff = pow(num, 2) / Region1::ddgamma_ddpi(temperature, pressure);

        return
            - pow(tau, 2 ) *
            Region1::ddgamma_ddtau(temperature, pressure) * R +
            diff;
    };

    // the unregularized specific internal energy for liquid water
    static Scalar internalEnergyRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            R * temperature *
            ( Region1::tau(temperature)*Region1::dgamma_dtau(temperature, pressure) -
              Region1::pi(pressure)*Region1::dgamma_dpi(temperature, pressure));
    };

    // the unregularized specific volume for liquid water
    static Scalar volumeRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            Region1::pi(pressure)*
            Region1::dgamma_dpi(temperature, pressure) *
            R * temperature / pressure;
    };

    // the unregularized specific enthalpy for steam
    static Scalar enthalpyRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            Region2::tau(temperature) *
            Region2::dgamma_dtau(temperature, pressure) *
            R*temperature;
    };

    // the unregularized specific internal energy for steam
    static Scalar internalEnergyRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            R * temperature *
            ( Region2::tau(temperature)*Region2::dgamma_dtau(temperature, pressure) -
              Region2::pi(pressure)*Region2::dgamma_dpi(temperature, pressure));
    };

    // the unregularized specific isobaric heat capacity
    static Scalar heatCap_p_Region2_(Scalar temperature, Scalar pressure)
    {
        return
            - pow(Region2::tau(temperature), 2 ) *
            Region2::ddgamma_ddtau(temperature, pressure) *
            R;
    };

    // the unregularized specific isochoric heat capacity
    static Scalar heatCap_v_Region2_(Scalar temperature, Scalar pressure)
    {
        double tau = Region2::tau(temperature);
        double pi = Region2::pi(pressure);
        double num = 1 + pi * Region2::dgamma_dpi(temperature, pressure) + tau * pi * Region2::ddgamma_dtaudpi(temperature, pressure);
        double diff = num * num / (1 - pi * pi * Region2::ddgamma_ddpi(temperature, pressure));
        return
            - pow(tau, 2 ) *
            Region2::ddgamma_ddtau(temperature, pressure) * R
            - diff;
    };

    // the unregularized specific volume for steam
    static Scalar volumeRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            Region2::pi(pressure)*
            Region2::dgamma_dpi(temperature, pressure) *
            R * temperature / pressure;
    };
}; // end class

} // end namepace

#endif
