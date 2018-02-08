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
 * \copydoc Opm::H2O
 */
#ifndef OPM_H2O_HPP
#define OPM_H2O_HPP

#include "iapws/Common.hpp"
#include "iapws/Region1.hpp"
#include "iapws/Region2.hpp"
#include "iapws/Region4.hpp"

#include "Component.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <cmath>
#include <cassert>
#include <sstream>

namespace Opm {

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
    typedef IAPWS::Common<Scalar> Common;
    typedef IAPWS::Region1<Scalar> Region1;
    typedef IAPWS::Region2<Scalar> Region2;
    typedef IAPWS::Region4<Scalar> Region4;

    static const Scalar Rs; // specific gas constant of water

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char* name()
    { return "H2O"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static const Scalar molarMass()
    { return Common::molarMass; }

    /*!
     * \brief The acentric factor \f$\mathrm{[-]}\f$ of water.
     */
    static const Scalar acentricFactor()
    { return Common::acentricFactor; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water
     */
    static const Scalar criticalTemperature()
    { return Common::criticalTemperature; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static const Scalar criticalPressure()
    { return Common::criticalPressure; }

    /*!
     * \brief Returns the molar volume \f$\mathrm{[m^3/mol]}\f$ of water at the critical point
     */
    static const Scalar criticalMolarVolume()
    { return Common::criticalMolarVolume; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static const Scalar tripleTemperature()
    { return Common::tripleTemperature; }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static const Scalar triplePressure()
    { return Common::triplePressure; }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param T Absolute temperature of the system in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    static Evaluation vaporPressure(Evaluation temperature)
    {
        if (temperature > criticalTemperature())
            temperature = criticalTemperature();
        if (temperature < tripleTemperature())
            temperature = tripleTemperature();

        return Region4::saturationPressure(temperature);
    }
    /*!
     * \brief The vapor temperature in \f$\mathrm{[Ka]}\f$ of pure water
     *        at a given pressure.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation vaporTemperature(const Evaluation& pressure)
    {
        if (pressure > criticalPressure())
            pressure = criticalPressure();
        if (pressure < triplePressure())
            pressure = triplePressure();

        return Region4::vaporTemperature(pressure);
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Enthalpy of steam is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;

            throw NumericalIssue(oss.str());
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the gas enthalpy at very low pressures. The
            // enthalpy of an ideal gas does not exhibit any
            // dependence on pressure, so we can just return the
            // specific enthalpy at the point of regularization, i.e.
            // the triple pressure - 100Pa
            return enthalpyRegion2_<Evaluation>(temperature, triplePressure() - 100);
        }
        Evaluation pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            Evaluation dh_dp =
                Rs*temperature*
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
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Enthalpy of water is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;

            throw NumericalIssue(oss.str());
        }

        // regularization
        const Evaluation& pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            const Evaluation& dh_dp =
                Rs * temperature*
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
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature,
                                      const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Heat capacity of steam is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;

            throw NumericalIssue(oss.str());
        }

        // regularization
        if (pressure < triplePressure() - 100)
            return heatCap_p_Region2_(temperature, Evaluation(triplePressure() - 100));
        const Evaluation& pv = vaporPressure(temperature);
        if (pressure > pv)
            // the pressure is too high, in this case we use the heat
            // cap at the vapor pressure to regularize
            return heatCap_p_Region2_(temperature, pv);

        return heatCap_p_Region2_(temperature, pressure);
    }

    /*!
     * \brief Specific isobaric heat capacity of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature,
                                         const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Heat capacity of water is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        }

        // regularization
        const Evaluation& pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the heat capacity at the
            // vapor pressure to regularize
            return heatCap_p_Region1_(temperature, pv);
        };

        return heatCap_p_Region1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& temperature,
                                           const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Internal Energy of water is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;

            throw NumericalIssue(oss.str());
        }


        // regularization
        Scalar pv = vaporPressure<Scalar>(Opm::scalarValue(temperature));
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
            Rs*temperature*
            (tau*dpi_dp*ddgamma_dtaudpi + dpi_dp*dpi_dp*dgamma_dpi + pi*dpi_dp*ddgamma_ddpi);
            */

            // use a straight line for extrapolation. use forward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            Scalar eps = 1e-7;
            const Evaluation& uv = internalEnergyRegion1_(temperature, Evaluation(pv));
            const Evaluation& uvPEps = internalEnergyRegion1_(temperature, Evaluation(pv + eps));
            const Evaluation& du_dp = (uvPEps - uv)/eps;
            return uv + du_dp*(pressure - pv);
        };

        return internalEnergyRegion1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of steam and water vapor \f$\mathrm{[J/kg]}\f$.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss <<"Internal energy of steam is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
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
                enthalpyRegion2_(temperature, Evaluation(triplePressure() - 100.0))
                -
                Rs*temperature; // = p*v   for an ideal gas!
        }
        Scalar pv = vaporPressure(Opm::scalarValue(temperature));
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
            Rs*temperature*
            (tau*dpi_dp*ddgamma_dtaudpi + dpi_dp*dpi_dp*dgamma_dpi + pi*dpi_dp*ddgamma_ddpi);

            // use a straight line for extrapolation
            Scalar uv = internalEnergyRegion2_(temperature, pv);
            return uv + du_dp*(pressure - pv);
            */

            // use a straight line for extrapolation. use backward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            Scalar eps = 1e-7;
            const Evaluation& uv = internalEnergyRegion2_(temperature, Evaluation(pv));
            const Evaluation& uvMEps = internalEnergyRegion2_(temperature, Evaluation(pv - eps));
            const Evaluation& du_dp = (uv - uvMEps)/eps;
            return uv + du_dp*(pressure - pv);
        };

        return internalEnergyRegion2_(temperature, pressure);
    }

    /*!
     * \brief Specific isochoric heat capacity of liquid water \f$\mathrm{[J/m^3]}\f$.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacityConstVolume(const Evaluation& temperature,
                                                    const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Heat capacity of water is only implemented for temperatures below 623.15K and "
                      "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;

            throw NumericalIssue(oss.str());
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
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacityConstVolume(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Heat capacity of steam is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
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
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return true; }

    /*!
     * \brief The density of steam in \f$\mathrm{[kg/m^3]}\f$ at a given pressure and temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Density of steam is only implemented for temperatures below 623.15K and "
                      "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        }

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the internal energy and enthalpy.
            const Evaluation& rho0IAPWS =
                1.0/volumeRegion2_(temperature,
                                   Evaluation(triplePressure() - 100));
            const Evaluation& rho0Id =
                IdealGas<Scalar>::density(Evaluation(molarMass()),
                                          temperature,
                                          Evaluation(triplePressure() - 100));
            return
                rho0IAPWS/rho0Id
                *IdealGas<Scalar>::density(Evaluation(molarMass()),
                                           temperature,
                                           pressure);
        }
        Evaluation pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the density energy at the vapor pressure to
            // regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            Scalar eps = Opm::scalarValue(pv)*1e-8;
            Evaluation v0 = volumeRegion2_(temperature, pv);
            Evaluation v1 = volumeRegion2_(temperature, pv + eps);
            Evaluation dv_dp = (v1 - v0)/eps;
            /*
              Scalar pi = Region2::pi(pv);
              Scalar dp_dpi = Region2::dp_dpi(pv);
              Scalar dgamma_dpi = Region2::dgamma_dpi(temperature, pv);
              Scalar ddgamma_ddpi = Region2::ddgamma_ddpi(temperature, pv);

              Scalar RT = Rs*temperature;
              Scalar dv_dp =
              RT/(dp_dpi*pv)
              *
              (dgamma_dpi + pi*ddgamma_ddpi - v0*dp_dpi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Evaluation drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        };

        return 1.0/volumeRegion2_(temperature, pressure);
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return false; }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param density Density in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, Scalar density)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(density);

        // We use the newton method for this. For the initial value we
        // assume steam to be an ideal gas
        Evaluation pressure = IdealGas<Scalar>::pressure(temperature, density/molarMass());
        Scalar eps = pressure*1e-7;

        Evaluation deltaP = pressure*2;
        Valgrind::CheckDefined(pressure);
        Valgrind::CheckDefined(deltaP);
        for (int i = 0; i < 5 && std::abs(Opm::scalarValue(pressure)*1e-9) < std::abs(Opm::scalarValue(deltaP)); ++i) {
            Evaluation f = gasDensity(temperature, pressure) - density;

            Evaluation df_dp;
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
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Density of water is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        }

        // regularization
        Evaluation pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the density at the vapor pressure to regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            Scalar eps = Opm::scalarValue(pv)*1e-8;
            Evaluation v0 = volumeRegion1_(temperature, pv);
            Evaluation v1 = volumeRegion1_(temperature, pv + eps);
            Evaluation dv_dp = (v1 - v0)/eps;

            /*
              Scalar v0 = volumeRegion1_(temperature, pv);
              Scalar pi = Region1::pi(pv);
              Scalar dp_dpi = Region1::dp_dpi(pv);
              Scalar dgamma_dpi = Region1::dgamma_dpi(temperature, pv);
              Scalar ddgamma_ddpi = Region1::ddgamma_ddpi(temperature, pv);

              Scalar RT = Rs*temperature;
              Scalar dv_dp =
              RT/(dp_dpi*pv)
              *
              (dgamma_dpi + pi*ddgamma_ddpi - v0*dp_dpi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Evaluation drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        };

        return 1/volumeRegion1_(temperature, pressure);
    }

    /*!
     * \brief The pressure of liquid water in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param density Density of the fluid in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidPressure(const Evaluation& temperature, Scalar density)
    {
        // We use the Newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Evaluation pressure = 1.1*vaporPressure(temperature);
        Scalar eps = Opm::scalarValue(pressure)*1e-7;

        Evaluation deltaP = pressure*2;
        for (int i = 0; i < 5 && std::abs(Opm::scalarValue(pressure)*1e-9) < std::abs(Opm::scalarValue(deltaP)); ++i) {
            Evaluation f = liquidDensity(temperature, pressure) - density;

            Evaluation df_dp;
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
     * This method is only valid if pressure is below or at the vapor
     * pressure of water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Viscosity of steam is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        }

        Evaluation rho = gasDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     *
     * \param temperature Absolute temperature of the fluid in \f$\mathrm{[K]}\f$
     * \param pressure Phase pressure in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            std::ostringstream oss;
            oss << "Viscosity of water is only implemented for temperatures below 623.15K and "
                << "pressures below 100MPa. (T = " << temperature << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        };

        const Evaluation& rho = liquidDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m K)]}\f$ of water (IAPWS) .
     *
     * Implementation taken from:
     * freesteam - IAPWS-IF97 steam tables library
     * Copyright (C) 2004-2009  John Pye
     *
     * Appendix B: Recommended Interpolating equation for Industrial Use
     * see http://www.iapws.org/relguide/thcond.pdf
     *
     * \param temperature Absolute temperature in K
     * \param pressure Phase pressure of the phase in Pa
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& temperature,  const Evaluation& pressure)
    {
        const Evaluation& rho = liquidDensity(temperature, pressure);
        return Common::thermalConductivityIAPWS(temperature, rho);
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m K)]}\f$ of water (IAPWS) .
     *
     * Implementation taken from:
     * freesteam - IAPWS-IF97 steam tables library
     * Copyright (C) 2004-2009  John Pye
     *
     * Appendix B: Recommended Interpolating equation for Industrial Use
     * see http://www.iapws.org/relguide/thcond.pdf
     *
     * \param temperature Absolute temperature in K
     * \param pressure Phase pressure of the phase in Pa
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& rho = gasDensity(temperature, pressure);
        return Common::thermalConductivityIAPWS(temperature, rho);
    }

private:
    // the unregularized specific enthalpy for liquid water
    template <class Evaluation>
    static Evaluation enthalpyRegion1_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Region1::tau(temperature) *
            Region1::dgamma_dtau(temperature, pressure) *
            Rs*temperature;
    }

    // the unregularized specific isobaric heat capacity
    template <class Evaluation>
    static Evaluation heatCap_p_Region1_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            - Opm::pow(Region1::tau(temperature), 2.0) *
            Region1::ddgamma_ddtau(temperature, pressure) *
            Rs;
    }

    // the unregularized specific isochoric heat capacity
    template <class Evaluation>
    static Evaluation heatCap_v_Region1_(const Evaluation& temperature, const Evaluation& pressure)
    {
        double tau = Region1::tau(temperature);
        double num = Region1::dgamma_dpi(temperature, pressure) - tau * Region1::ddgamma_dtaudpi(temperature, pressure);
        double diff = std::pow(num, 2) / Region1::ddgamma_ddpi(temperature, pressure);

        return
            - std::pow(tau, 2 ) *
            Region1::ddgamma_ddtau(temperature, pressure) * Rs +
            diff;
    }

    // the unregularized specific internal energy for liquid water
    template <class Evaluation>
    static Evaluation internalEnergyRegion1_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Rs * temperature *
            ( Region1::tau(temperature)*Region1::dgamma_dtau(temperature, pressure) -
              Region1::pi(pressure)*Region1::dgamma_dpi(temperature, pressure));
    }

    // the unregularized specific volume for liquid water
    template <class Evaluation>
    static Evaluation volumeRegion1_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Region1::pi(pressure)*
            Region1::dgamma_dpi(temperature, pressure) *
            Rs * temperature / pressure;
    }

    // the unregularized specific enthalpy for steam
    template <class Evaluation>
    static Evaluation enthalpyRegion2_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Region2::tau(temperature) *
            Region2::dgamma_dtau(temperature, pressure) *
            Rs*temperature;
    }

    // the unregularized specific internal energy for steam
    template <class Evaluation>
    static Evaluation internalEnergyRegion2_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Rs * temperature *
            ( Region2::tau(temperature)*Region2::dgamma_dtau(temperature, pressure) -
              Region2::pi(pressure)*Region2::dgamma_dpi(temperature, pressure));
    }

    // the unregularized specific isobaric heat capacity
    template <class Evaluation>
    static Evaluation heatCap_p_Region2_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            - Opm::pow(Region2::tau(temperature), 2 ) *
            Region2::ddgamma_ddtau(temperature, pressure) *
            Rs;
    }

    // the unregularized specific isochoric heat capacity
    template <class Evaluation>
    static Evaluation heatCap_v_Region2_(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& tau = Region2::tau(temperature);
        const Evaluation& pi = Region2::pi(pressure);
        const Evaluation& num = 1 + pi * Region2::dgamma_dpi(temperature, pressure) + tau * pi * Region2::ddgamma_dtaudpi(temperature, pressure);
        const Evaluation& diff = num * num / (1 - pi * pi * Region2::ddgamma_ddpi(temperature, pressure));
        return
            - std::pow(tau, 2 ) *
            Region2::ddgamma_ddtau(temperature, pressure) * Rs
            - diff;
    }

    // the unregularized specific volume for steam
    template <class Evaluation>
    static Evaluation volumeRegion2_(const Evaluation& temperature, const Evaluation& pressure)
    {
        return
            Region2::pi(pressure)*
            Region2::dgamma_dpi(temperature, pressure) *
            Rs * temperature / pressure;
    }
}; // end class

template <class Scalar>
const Scalar H2O<Scalar>::Rs = Common::Rs;
} // namespace Opm

#endif
