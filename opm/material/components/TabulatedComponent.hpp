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
 * \copydoc Opm::TabulatedComponent
 */
#ifndef OPM_TABULATED_COMPONENT_HPP
#define OPM_TABULATED_COMPONENT_HPP

#include <cmath>
#include <limits>
#include <cassert>
#include <iostream>

#include <opm/material/common/MathToolbox.hpp>

namespace Opm {
/*!
 * \ingroup Components
 *
 * \brief A generic class which tabulates all thermodynamic properties
 *        of a given component.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam RawComponent The component which ought to be tabulated
 * \tparam useVaporPressure If true, tabulate all quantities along the
 *                          vapor pressure curve, if false use the
 *                          pressure range [p_min, p_max]
 */
template <class ScalarT, class RawComponent, bool useVaporPressure=true>
class TabulatedComponent
{
public:
    typedef ScalarT Scalar;

    static const bool isTabulated = true;

    /*!
     * \brief Initialize the tables.
     *
     * \param tempMin The minimum of the temperature range in \f$\mathrm{[K]}\f$
     * \param tempMax The maximum of the temperature range in \f$\mathrm{[K]}\f$
     * \param nTemp The number of entries/steps within the temperature range
     * \param pressMin The minimum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param nPress The number of entries/steps within the pressure range
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        tempMin_ = tempMin;
        tempMax_ = tempMax;
        nTemp_ = nTemp;
        pressMin_ = pressMin;
        pressMax_ = pressMax;
        nPress_ = nPress;
        nDensity_ = nPress_;

        // allocate the arrays
        vaporPressure_ = new Scalar[nTemp_];
        minGasDensity__ = new Scalar[nTemp_];
        maxGasDensity__ = new Scalar[nTemp_];
        minLiquidDensity__ = new Scalar[nTemp_];
        maxLiquidDensity__ = new Scalar[nTemp_];

        gasEnthalpy_ = new Scalar[nTemp_*nPress_];
        liquidEnthalpy_ = new Scalar[nTemp_*nPress_];
        gasHeatCapacity_ = new Scalar[nTemp_*nPress_];
        liquidHeatCapacity_ = new Scalar[nTemp_*nPress_];
        gasDensity_ = new Scalar[nTemp_*nPress_];
        liquidDensity_ = new Scalar[nTemp_*nPress_];
        gasViscosity_ = new Scalar[nTemp_*nPress_];
        liquidViscosity_ = new Scalar[nTemp_*nPress_];
        gasThermalConductivity_ = new Scalar[nTemp_*nPress_];
        liquidThermalConductivity_ = new Scalar[nTemp_*nPress_];
        gasPressure_ = new Scalar[nTemp_*nDensity_];
        liquidPressure_ = new Scalar[nTemp_*nDensity_];

        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        Scalar NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // fill the temperature-pressure arrays
        for (unsigned iT = 0; iT < nTemp_; ++ iT) {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            try { vaporPressure_[iT] = RawComponent::vaporPressure(temperature); }
            catch (const std::exception&) { vaporPressure_[iT] = NaN; }

            Scalar pgMax = maxGasPressure_(iT);
            Scalar pgMin = minGasPressure_(iT);

            // fill the temperature, pressure gas arrays
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (pgMax - pgMin)/(nPress_ - 1) + pgMin;

                unsigned i = iT + iP*nTemp_;

                try { gasEnthalpy_[i] = RawComponent::gasEnthalpy(temperature, pressure); }
                catch (const std::exception&) { gasEnthalpy_[i] = NaN; }

                try { gasHeatCapacity_[i] = RawComponent::gasHeatCapacity(temperature, pressure); }
                catch (const std::exception&) { gasHeatCapacity_[i] = NaN; }

                try { gasDensity_[i] = RawComponent::gasDensity(temperature, pressure); }
                catch (const std::exception&) { gasDensity_[i] = NaN; }

                try { gasViscosity_[i] = RawComponent::gasViscosity(temperature, pressure); }
                catch (const std::exception&) { gasViscosity_[i] = NaN; }

                try { gasThermalConductivity_[i] = RawComponent::gasThermalConductivity(temperature, pressure); }
                catch (const std::exception&) { gasThermalConductivity_[i] = NaN; }
            };

            Scalar plMin = minLiquidPressure_(iT);
            Scalar plMax = maxLiquidPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (plMax - plMin)/(nPress_ - 1) + plMin;

                unsigned i = iT + iP*nTemp_;

                try { liquidEnthalpy_[i] = RawComponent::liquidEnthalpy(temperature, pressure); }
                catch (const std::exception&) { liquidEnthalpy_[i] = NaN; }

                try { liquidHeatCapacity_[i] = RawComponent::liquidHeatCapacity(temperature, pressure); }
                catch (const std::exception&) { liquidHeatCapacity_[i] = NaN; }

                try { liquidDensity_[i] = RawComponent::liquidDensity(temperature, pressure); }
                catch (const std::exception&) { liquidDensity_[i] = NaN; }

                try { liquidViscosity_[i] = RawComponent::liquidViscosity(temperature, pressure); }
                catch (const std::exception&) { liquidViscosity_[i] = NaN; }

                try { liquidThermalConductivity_[i] = RawComponent::liquidThermalConductivity(temperature, pressure); }
                catch (const std::exception&) { liquidThermalConductivity_[i] = NaN; }
            }
        }

        // fill the temperature-density arrays
        for (unsigned iT = 0; iT < nTemp_; ++ iT) {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            // calculate the minimum and maximum values for the gas
            // densities
            minGasDensity__[iT] = RawComponent::gasDensity(temperature, minGasPressure_(iT));
            if (iT < nTemp_ - 1)
                maxGasDensity__[iT] = RawComponent::gasDensity(temperature, maxGasPressure_(iT + 1));
            else
                maxGasDensity__[iT] = RawComponent::gasDensity(temperature, maxGasPressure_(iT));

            // fill the temperature, density gas arrays
            for (unsigned iRho = 0; iRho < nDensity_; ++ iRho) {
                Scalar density =
                    Scalar(iRho)/(nDensity_ - 1) *
                    (maxGasDensity__[iT] - minGasDensity__[iT])
                    +
                    minGasDensity__[iT];

                unsigned i = iT + iRho*nTemp_;

                try { gasPressure_[i] = RawComponent::gasPressure(temperature, density); }
                catch (const std::exception&) { gasPressure_[i] = NaN; };
            };

            // calculate the minimum and maximum values for the liquid
            // densities
            minLiquidDensity__[iT] = RawComponent::liquidDensity(temperature, minLiquidPressure_(iT));
            if (iT < nTemp_ - 1)
                maxLiquidDensity__[iT] = RawComponent::liquidDensity(temperature, maxLiquidPressure_(iT + 1));
            else
                maxLiquidDensity__[iT] = RawComponent::liquidDensity(temperature, maxLiquidPressure_(iT));

            // fill the temperature, density liquid arrays
            for (unsigned iRho = 0; iRho < nDensity_; ++ iRho) {
                Scalar density =
                    Scalar(iRho)/(nDensity_ - 1) *
                    (maxLiquidDensity__[iT] - minLiquidDensity__[iT])
                    +
                    minLiquidDensity__[iT];

                unsigned i = iT + iRho*nTemp_;

                try { liquidPressure_[i] = RawComponent::liquidPressure(temperature, density); }
                catch (const std::exception&) { liquidPressure_[i] = NaN; };
            };
        }
    }

    /*!
     * \brief A human readable name for the component.
     */
    static const char* name()
    { return RawComponent::name(); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static Scalar molarMass()
    { return RawComponent::molarMass(); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    static Scalar criticalTemperature()
    { return RawComponent::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    static Scalar criticalPressure()
    { return RawComponent::criticalPressure(); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature()
    { return RawComponent::tripleTemperature(); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    { return RawComponent::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     *
     * \param T temperature of component
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& temperature)
    {
        const Evaluation& result = interpolateT_(vaporPressure_, temperature);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::vaporPressure(temperature);
        return result;
    }

    /*!
     * \brief Specific enthalpy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateGasTP_(gasEnthalpy_,
                                                     temperature,
                                                     pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasEnthalpy(temperature, pressure);
        return result;
    }

    /*!
     * \brief Specific enthalpy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateLiquidTP_(liquidEnthalpy_,
                                                        temperature,
                                                        pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidEnthalpy(temperature, pressure);
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the gas \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateGasTP_(gasHeatCapacity_,
                                                     temperature,
                                                     pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasHeatCapacity(temperature, pressure);
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the liquid \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateLiquidTP_(liquidHeatCapacity_,
                                                        temperature,
                                                        pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidHeatCapacity(temperature, pressure);
        return result;
    }

    /*!
     * \brief Specific internal energy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature, const Evaluation& pressure)
    { return gasEnthalpy(temperature, pressure) - pressure/gasDensity(temperature, pressure); }

    /*!
     * \brief Specific internal energy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& temperature, const Evaluation& pressure)
    { return liquidEnthalpy(temperature, pressure) - pressure/liquidDensity(temperature, pressure); }

    /*!
     * \brief The pressure of gas in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, Scalar density)
    {
        const Evaluation& result = interpolateGasTRho_(gasPressure_,
                                                       temperature,
                                                       density);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasPressure(temperature,
                                             density);
        return result;
    }

    /*!
     * \brief The pressure of liquid in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidPressure(const Evaluation& temperature, Scalar density)
    {
        const Evaluation& result = interpolateLiquidTRho_(liquidPressure_,
                                                          temperature,
                                                          density);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidPressure(temperature,
                                                density);
        return result;
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return RawComponent::gasIsCompressible(); }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return RawComponent::liquidIsCompressible(); }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return RawComponent::gasIsIdeal(); }


    /*!
     * \brief The density of gas at a given pressure and temperature
     *        \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateGasTP_(gasDensity_,
                                                     temperature,
                                                     pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasDensity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The density of liquid at a given pressure and
     *        temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateLiquidTP_(liquidDensity_,
                                                        temperature,
                                                        pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidDensity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateGasTP_(gasViscosity_,
                                                     temperature,
                                                     pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasViscosity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateLiquidTP_(liquidViscosity_,
                                                        temperature,
                                                        pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidViscosity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The thermal conductivity of gaseous water \f$\mathrm{[W / (m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateGasTP_(gasThermalConductivity_,
                                                     temperature,
                                                     pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::gasThermalConductivity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The thermal conductivity of liquid water \f$\mathrm{[W / (m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& temperature, const Evaluation& pressure)
    {
        const Evaluation& result = interpolateLiquidTP_(liquidThermalConductivity_,
                                                        temperature,
                                                        pressure);
        if (std::isnan(Opm::scalarValue(result)))
            return RawComponent::liquidThermalConductivity(temperature, pressure);
        return result;
    }

private:
    // returns an interpolated value depending on temperature
    template <class Evaluation>
    static Evaluation interpolateT_(const Scalar* values, const Evaluation& T)
    {
        Evaluation alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1)
            return std::numeric_limits<Scalar>::quiet_NaN();

        size_t iT = static_cast<size_t>(Opm::scalarValue(alphaT));
        alphaT -= iT;

        return
            values[iT    ]*(1 - alphaT) +
            values[iT + 1]*(    alphaT);
    }

    // returns an interpolated value for liquid depending on
    // temperature and pressure
    template <class Evaluation>
    static Evaluation interpolateLiquidTP_(const Scalar* values, const Evaluation& T, const Evaluation& p)
    {
        Evaluation alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1)
            return std::numeric_limits<Scalar>::quiet_NaN();

        size_t iT = static_cast<size_t>(Opm::scalarValue(alphaT));
        alphaT -= iT;

        Evaluation alphaP1 = pressLiquidIdx_(p, iT);
        Evaluation alphaP2 = pressLiquidIdx_(p, iT + 1);

        size_t iP1 =
            static_cast<size_t>(
                std::max<int>(0, std::min(static_cast<int>(nPress_) - 2,
                                          static_cast<int>(Opm::scalarValue(alphaP1)))));
        size_t iP2 =
            static_cast<size_t>(
                std::max(0, std::min(static_cast<int>(nPress_) - 2,
                                     static_cast<int>(Opm::scalarValue(alphaP2)))));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    // returns an interpolated value for gas depending on
    // temperature and pressure
    template <class Evaluation>
    static Evaluation interpolateGasTP_(const Scalar* values, const Evaluation& T, const Evaluation& p)
    {
        Evaluation alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1)
            return std::numeric_limits<Scalar>::quiet_NaN();

        size_t iT =
            static_cast<size_t>(
                std::max(0, std::min(static_cast<int>(nTemp_) - 2,
                                     static_cast<int>(Opm::scalarValue(alphaT)))));
        alphaT -= iT;

        Evaluation alphaP1 = pressGasIdx_(p, iT);
        Evaluation alphaP2 = pressGasIdx_(p, iT + 1);
        size_t iP1 =
            static_cast<size_t>(
                std::max(0, std::min(static_cast<int>(nPress_) - 2,
                                     static_cast<int>(Opm::scalarValue(alphaP1)))));
        size_t iP2 =
            static_cast<size_t>(
                std::max(0, std::min(static_cast<int>(nPress_) - 2,
                                     static_cast<int>(Opm::scalarValue(alphaP2)))));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    // returns an interpolated value for gas depending on
    // temperature and density
    template <class Evaluation>
    static Evaluation interpolateGasTRho_(const Scalar* values, const Evaluation& T, const Evaluation& rho)
    {
        Evaluation alphaT = tempIdx_(T);
        unsigned iT = std::max(0,
                               std::min(static_cast<int>(nTemp_ - 2),
                                        static_cast<int>(alphaT)));
        alphaT -= iT;

        Evaluation alphaP1 = densityGasIdx_(rho, iT);
        Evaluation alphaP2 = densityGasIdx_(rho, iT + 1);
        unsigned iP1 =
            std::max(0,
                     std::min(static_cast<int>(nDensity_ - 2),
                              static_cast<int>(alphaP1)));
        unsigned iP2 =
            std::max(0,
                     std::min(static_cast<int>(nDensity_ - 2),
                              static_cast<int>(alphaP2)));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    // returns an interpolated value for liquid depending on
    // temperature and density
    template <class Evaluation>
    static Evaluation interpolateLiquidTRho_(const Scalar* values, const Evaluation& T, const Evaluation& rho)
    {
        Evaluation alphaT = tempIdx_(T);
        unsigned iT = std::max<int>(0, std::min<int>(nTemp_ - 2, static_cast<int>(alphaT)));
        alphaT -= iT;

        Evaluation alphaP1 = densityLiquidIdx_(rho, iT);
        Evaluation alphaP2 = densityLiquidIdx_(rho, iT + 1);
        unsigned iP1 = std::max<int>(0, std::min<int>(nDensity_ - 2, static_cast<int>(alphaP1)));
        unsigned iP2 = std::max<int>(0, std::min<int>(nDensity_ - 2, static_cast<int>(alphaP2)));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }


    // returns the index of an entry in a temperature field
    template <class Evaluation>
    static Evaluation tempIdx_(const Evaluation& temperature)
    {
        return (nTemp_ - 1)*(temperature - tempMin_)/(tempMax_ - tempMin_);
    }

    // returns the index of an entry in a pressure field
    template <class Evaluation>
    static Evaluation pressLiquidIdx_(const Evaluation& pressure, size_t tempIdx)
    {
        Scalar plMin = minLiquidPressure_(tempIdx);
        Scalar plMax = maxLiquidPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - plMin)/(plMax - plMin);
    }

    // returns the index of an entry in a temperature field
    template <class Evaluation>
    static Evaluation pressGasIdx_(const Evaluation& pressure, size_t tempIdx)
    {
        Scalar pgMin = minGasPressure_(tempIdx);
        Scalar pgMax = maxGasPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - pgMin)/(pgMax - pgMin);
    }

    // returns the index of an entry in a density field
    template <class Evaluation>
    static Evaluation densityLiquidIdx_(const Evaluation& density, size_t tempIdx)
    {
        Scalar densityMin = minLiquidDensity_(tempIdx);
        Scalar densityMax = maxLiquidDensity_(tempIdx);
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    // returns the index of an entry in a density field
    template <class Evaluation>
    static Evaluation densityGasIdx_(const Evaluation& density, size_t tempIdx)
    {
        Scalar densityMin = minGasDensity_(tempIdx);
        Scalar densityMax = maxGasDensity_(tempIdx);
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    // returns the minimum tabulized liquid pressure at a given
    // temperature index
    static Scalar minLiquidPressure_(size_t tempIdx)
    {
        if (!useVaporPressure)
            return pressMin_;
        else
            return std::max<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1);
    }

    // returns the maximum tabulized liquid pressure at a given
    // temperature index
    static Scalar maxLiquidPressure_(size_t tempIdx)
    {
        if (!useVaporPressure)
            return pressMax_;
        else
            return std::max<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    // returns the minumum tabulized gas pressure at a given
    // temperature index
    static Scalar minGasPressure_(size_t tempIdx)
    {
        if (!useVaporPressure)
            return pressMin_;
        else
            return std::min<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1 );
    }

    // returns the maximum tabulized gas pressure at a given
    // temperature index
    static Scalar maxGasPressure_(size_t tempIdx)
    {
        if (!useVaporPressure)
            return pressMax_;
        else
            return std::min<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }


    // returns the minimum tabulized liquid density at a given
    // temperature index
    static Scalar minLiquidDensity_(size_t tempIdx)
    { return minLiquidDensity__[tempIdx]; }

    // returns the maximum tabulized liquid density at a given
    // temperature index
    static Scalar maxLiquidDensity_(size_t tempIdx)
    { return maxLiquidDensity__[tempIdx]; }

    // returns the minumum tabulized gas density at a given
    // temperature index
    static Scalar minGasDensity_(size_t tempIdx)
    { return minGasDensity__[tempIdx]; }

    // returns the maximum tabulized gas density at a given
    // temperature index
    static Scalar maxGasDensity_(size_t tempIdx)
    { return maxGasDensity__[tempIdx]; }

    // 1D fields with the temperature as degree of freedom
    static Scalar* vaporPressure_;

    static Scalar* minLiquidDensity__;
    static Scalar* maxLiquidDensity__;

    static Scalar* minGasDensity__;
    static Scalar* maxGasDensity__;

    // 2D fields with the temperature and pressure as degrees of
    // freedom
    static Scalar* gasEnthalpy_;
    static Scalar* liquidEnthalpy_;

    static Scalar* gasHeatCapacity_;
    static Scalar* liquidHeatCapacity_;

    static Scalar* gasDensity_;
    static Scalar* liquidDensity_;

    static Scalar* gasViscosity_;
    static Scalar* liquidViscosity_;

    static Scalar* gasThermalConductivity_;
    static Scalar* liquidThermalConductivity_;

    // 2D fields with the temperature and density as degrees of
    // freedom
    static Scalar* gasPressure_;
    static Scalar* liquidPressure_;

    // temperature, pressure and density ranges
    static Scalar tempMin_;
    static Scalar tempMax_;
    static unsigned nTemp_;

    static Scalar pressMin_;
    static Scalar pressMax_;
    static unsigned nPress_;

    static Scalar densityMin_;
    static Scalar densityMax_;
    static unsigned nDensity_;
};

template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::vaporPressure_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::minLiquidDensity__;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::maxLiquidDensity__;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::minGasDensity__;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::maxGasDensity__;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasEnthalpy_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidEnthalpy_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasHeatCapacity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidHeatCapacity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasDensity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidDensity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasViscosity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidViscosity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasThermalConductivity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidThermalConductivity_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::gasPressure_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar* TabulatedComponent<Scalar, RawComponent, useVaporPressure>::liquidPressure_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::tempMin_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::tempMax_;
template <class Scalar, class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<Scalar, RawComponent, useVaporPressure>::nTemp_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::pressMin_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::pressMax_;
template <class Scalar, class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<Scalar, RawComponent, useVaporPressure>::nPress_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::densityMin_;
template <class Scalar, class RawComponent, bool useVaporPressure>
Scalar TabulatedComponent<Scalar, RawComponent, useVaporPressure>::densityMax_;
template <class Scalar, class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<Scalar, RawComponent, useVaporPressure>::nDensity_;


} // namespace Opm

#endif
