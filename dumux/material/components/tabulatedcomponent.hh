// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
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
 * \copydoc Dumux::TabulatedComponent
 */
#ifndef DUMUX_TABULATED_COMPONENT_HH
#define DUMUX_TABULATED_COMPONENT_HH

#include <cmath>
#include <limits>
#include <cassert>
#include <iostream>

#include <dumux/common/exceptions.hh>

namespace Dumux
{

/*!
 * \ingroup Components
 *
 * \brief A generic class which tabulates all thermodynamic properties
 *        of a given component.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 *
 * \tparam Scalar  The type used for scalar values
 * \tparam Scalar  The component which ought to be tabulated
 * \tparam useVaporPressure If true, tabulate all quantities along the
 *                          vapor pressure curve, if false use the pressure range [p_min, p_max]
 */
template <class Scalar, class RawComponent, bool useVaporPressure=true>
class TabulatedComponent
{
public:
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
#ifndef NDEBUG
        initialized_  = true;
        warningPrinted_ = false;
#endif
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
            catch (Dune::NotImplemented) { vaporPressure_[iT] = NaN; }
            catch (NumericalProblem e) { vaporPressure_[iT] = NaN; };

            Scalar pgMax = maxGasPressure_(iT);
            Scalar pgMin = minGasPressure_(iT);

            // fill the temperature, pressure gas arrays
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (pgMax - pgMin)/(nPress_ - 1) + pgMin;

                unsigned i = iT + iP*nTemp_;

                try { gasEnthalpy_[i] = RawComponent::gasEnthalpy(temperature, pressure); }
                catch (Dune::NotImplemented) { gasEnthalpy_[i] = NaN; }
                catch (NumericalProblem) { gasEnthalpy_[i] = NaN; };

                try { gasHeatCapacity_[i] = RawComponent::gasHeatCapacity(temperature, pressure); }
                catch (Dune::NotImplemented) { gasHeatCapacity_[i] = NaN; }
                catch (NumericalProblem) { gasHeatCapacity_[i] = NaN; };

                try { gasDensity_[i] = RawComponent::gasDensity(temperature, pressure); }
                catch (Dune::NotImplemented) { gasDensity_[i] = NaN; }
                catch (NumericalProblem) { gasDensity_[i] = NaN; };

                try { gasViscosity_[i] = RawComponent::gasViscosity(temperature, pressure); }
                catch (Dune::NotImplemented) { gasViscosity_[i] = NaN; }
                catch (NumericalProblem) { gasViscosity_[i] = NaN; };

                try { gasThermalConductivity_[i] = RawComponent::gasThermalConductivity(temperature, pressure); }
                catch (Dune::NotImplemented) { gasThermalConductivity_[i] = NaN; }
                catch (NumericalProblem) { gasThermalConductivity_[i] = NaN; };
            };

            Scalar plMin = minLiquidPressure_(iT);
            Scalar plMax = maxLiquidPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (plMax - plMin)/(nPress_ - 1) + plMin;

                unsigned i = iT + iP*nTemp_;

                try { liquidEnthalpy_[i] = RawComponent::liquidEnthalpy(temperature, pressure); }
                catch (Dune::NotImplemented) { liquidEnthalpy_[i] = NaN; }
                catch (NumericalProblem) { liquidEnthalpy_[i] = NaN; };

                try { liquidHeatCapacity_[i] = RawComponent::liquidHeatCapacity(temperature, pressure); }
                catch (Dune::NotImplemented) { liquidHeatCapacity_[i] = NaN; }
                catch (NumericalProblem) { liquidHeatCapacity_[i] = NaN; };

                try { liquidDensity_[i] = RawComponent::liquidDensity(temperature, pressure); }
                catch (Dune::NotImplemented) { liquidDensity_[i] = NaN; }
                catch (NumericalProblem) { liquidDensity_[i] = NaN; };

                try { liquidViscosity_[i] = RawComponent::liquidViscosity(temperature, pressure); }
                catch (Dune::NotImplemented) { liquidViscosity_[i] = NaN; }
                catch (NumericalProblem) { liquidViscosity_[i] = NaN; };

                try { liquidThermalConductivity_[i] = RawComponent::liquidThermalConductivity(temperature, pressure); }
                catch (Dune::NotImplemented) { liquidThermalConductivity_[i] = NaN; }
                catch (NumericalProblem) { liquidThermalConductivity_[i] = NaN; };
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
                catch (NumericalProblem) { gasPressure_[i] = NaN; };
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
                catch (NumericalProblem) { liquidPressure_[i] = NaN; };
            };
        }
    }

    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { return RawComponent::name(); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static constexpr Scalar molarMass()
    { return RawComponent::molarMass(); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    static constexpr Scalar criticalTemperature()
    { return RawComponent::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    static constexpr Scalar criticalPressure()
    { return RawComponent::criticalPressure(); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static constexpr Scalar tripleTemperature()
    { return RawComponent::tripleTemperature(); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static constexpr Scalar triplePressure()
    { return RawComponent::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     *
     * \param T temperature of component
     */
    static Scalar vaporPressure(Scalar T)
    {
        Scalar result = interpolateT_(vaporPressure_, T);
        if (std::isnan(result)) {
            return RawComponent::vaporPressure(T);
        }
        return result;
    }

    /*!
     * \brief Specific enthalpy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasEnthalpy_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            printWarning_("gasEnthalpy", temperature, pressure);
            return RawComponent::gasEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific enthalpy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidEnthalpy_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            printWarning_("liquidEnthalpy", temperature, pressure);
            return RawComponent::liquidEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the gas \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasHeatCapacity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            printWarning_("gasHeatCapacity", temperature, pressure);
            return RawComponent::gasHeatCapacity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the liquid \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidHeatCapacity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            printWarning_("liquidHeatCapacity", temperature, pressure);
            return RawComponent::liquidHeatCapacity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific internal energy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        Scalar result =
            gasEnthalpy(temperature, pressure) - pressure/gasDensity(temperature, pressure);
        return result;
    }

    /*!
     * \brief Specific internal energy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        Scalar result =
            liquidEnthalpy(temperature, pressure) - pressure/liquidDensity(temperature, pressure);
        return result;
    }

    /*!
     * \brief The pressure of gas in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateGasTRho_(gasPressure_,
                                            temperature,
                                            density);
        if (std::isnan(result)) {
            printWarning_("gasPressure", temperature, density);
            return RawComponent::gasPressure(temperature,
                                             density);
        }
        return result;
    }

    /*!
     * \brief The pressure of liquid in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateLiquidTRho_(liquidPressure_,
                                               temperature,
                                               density);
        if (std::isnan(result)) {
            printWarning_("liquidPressure", temperature, density);
            return RawComponent::liquidPressure(temperature,
                                                density);
        }
        return result;
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return RawComponent::gasIsCompressible(); }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return RawComponent::liquidIsCompressible(); }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return RawComponent::gasIsIdeal(); }


    /*!
     * \brief The density of gas at a given pressure and temperature
     *        \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasDensity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            printWarning_("gasDensity", temperature, pressure);
            return RawComponent::gasDensity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The density of liquid at a given pressure and
     *        temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidDensity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            printWarning_("liquidDensity", temperature, pressure);
            return RawComponent::liquidDensity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasViscosity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            printWarning_("gasViscosity", temperature, pressure);
            return RawComponent::gasViscosity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidViscosity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            printWarning_("liquidViscosity",temperature, pressure);
            return RawComponent::liquidViscosity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The thermal conductivity of gaseous water \f$\mathrm{[W / (m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasThermalConductivity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            printWarning_("gasThermalConductivity", temperature, pressure);
            return RawComponent::gasThermalConductivity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The thermal conductivity of liquid water \f$\mathrm{[W / (m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidThermalConductivity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            printWarning_("liquidThermalConductivity", temperature, pressure);
            return RawComponent::liquidThermalConductivity(temperature, pressure);
        }
        return result;
    }


private:
    // prints a warning if the result is not in range or the table has
    // not been initialized
    static void printWarning_(const char *quantity, Scalar arg1, Scalar arg2)
    {
#ifndef NDEBUG
        if (warningPrinted_)
            return;

        if (!initialized_)
            std::cerr << "TABULATED COMPONENT '" << name()
                      << "' WAS NOT INITIALIZED! "
                      << "PLEASE CALL FluidSystem::init()\n";
        else
            std::cerr << "FORWARD METHOD CALL "<<quantity<<"("<<arg1<<", "<<arg2<<") OF COMPONENT '"<<name()<<"'. TABULATION TOO SMALL?\n";
        warningPrinted_ = true;
#endif
    }


    // returns an interpolated value depending on temperature
    static Scalar interpolateT_(const Scalar *values, Scalar T)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1)
            return std::numeric_limits<Scalar>::quiet_NaN();

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        return
            values[iT    ]*(1 - alphaT) +
            values[iT + 1]*(    alphaT);
    }

    // returns an interpolated value for liquid depending on
    // temperature and pressure
    static Scalar interpolateLiquidTP_(const Scalar *values, Scalar T, Scalar p)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iT = std::max<int>(0, std::min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = pressLiquidIdx_(p, iT);
        Scalar alphaP2 = pressLiquidIdx_(p, iT + 1);

        unsigned iP1 = std::max<int>(0, std::min<int>(nPress_ - 2, (int) alphaP1));
        unsigned iP2 = std::max<int>(0, std::min<int>(nPress_ - 2, (int) alphaP2));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

#if 0 && !defined NDEBUG
        if(!(0 <= alphaT && alphaT <= 1.0))
            DUNE_THROW(NumericalProblem, "Temperature out of range: "
                       << "T=" << T << " range: [" << tempMin_ << ", " << tempMax_ << "]");
        if(!(0 <= alphaP1 && alphaP1 <= 1.0))
            DUNE_THROW(NumericalProblem, "First liquid pressure out of range: "
                       << "p=" << p << " range: [" << minLiquidPressure_(tempIdx_(T)) << ", " << maxLiquidPressure_(tempIdx_(T)) << "]");
        if(!(0 <= alphaP2 && alphaP2 <= 1.0))
            DUNE_THROW(NumericalProblem, "Second liquid pressure out of range: "
                       << "p=" << p << " range: [" << minLiquidPressure_(tempIdx_(T) + 1) << ", " << maxLiquidPressure_(tempIdx_(T) + 1) << "]");
#endif

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    // returns an interpolated value for gas depending on
    // temperature and pressure
    static Scalar interpolateGasTP_(const Scalar *values, Scalar T, Scalar p)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iT = std::max<int>(0, std::min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = pressGasIdx_(p, iT);
        Scalar alphaP2 = pressGasIdx_(p, iT + 1);
        unsigned iP1 = std::max<int>(0, std::min<int>(nPress_ - 2, (int) alphaP1));
        unsigned iP2 = std::max<int>(0, std::min<int>(nPress_ - 2, (int) alphaP2));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

#if 0 && !defined NDEBUG
        if(!(0 <= alphaT && alphaT <= 1.0))
            DUNE_THROW(NumericalProblem, "Temperature out of range: "
                       << "T=" << T << " range: [" << tempMin_ << ", " << tempMax_ << "]");
        if(!(0 <= alphaP1 && alphaP1 <= 1.0))
            DUNE_THROW(NumericalProblem, "First gas pressure out of range: "
                       << "p=" << p << " range: [" << minGasPressure_(tempIdx_(T)) << ", " << maxGasPressure_(tempIdx_(T)) << "]");
        if(!(0 <= alphaP2 && alphaP2 <= 1.0))
            DUNE_THROW(NumericalProblem, "Second gas pressure out of range: "
                       << "p=" << p << " range: [" << minGasPressure_(tempIdx_(T) + 1) << ", " << maxGasPressure_(tempIdx_(T) + 1) << "]");
#endif

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    // returns an interpolated value for gas depending on
    // temperature and density
    static Scalar interpolateGasTRho_(const Scalar *values, Scalar T, Scalar rho)
    {
        Scalar alphaT = tempIdx_(T);
        unsigned iT = std::max<int>(0, std::min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = densityGasIdx_(rho, iT);
        Scalar alphaP2 = densityGasIdx_(rho, iT + 1);
        unsigned iP1 = std::max<int>(0, std::min<int>(nDensity_ - 2, (int) alphaP1));
        unsigned iP2 = std::max<int>(0, std::min<int>(nDensity_ - 2, (int) alphaP2));
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
    static Scalar interpolateLiquidTRho_(const Scalar *values, Scalar T, Scalar rho)
    {
        Scalar alphaT = tempIdx_(T);
        unsigned iT = std::max<int>(0, std::min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = densityLiquidIdx_(rho, iT);
        Scalar alphaP2 = densityLiquidIdx_(rho, iT + 1);
        unsigned iP1 = std::max<int>(0, std::min<int>(nDensity_ - 2, (int) alphaP1));
        unsigned iP2 = std::max<int>(0, std::min<int>(nDensity_ - 2, (int) alphaP2));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return
            values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
            values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
            values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
            values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }


    // returns the index of an entry in a temperature field
    static Scalar tempIdx_(Scalar temperature)
    {
        return (nTemp_ - 1)*(temperature - tempMin_)/(tempMax_ - tempMin_);
    }

    // returns the index of an entry in a pressure field
    static Scalar pressLiquidIdx_(Scalar pressure, unsigned tempIdx)
    {
        Scalar plMin = minLiquidPressure_(tempIdx);
        Scalar plMax = maxLiquidPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - plMin)/(plMax - plMin);
    }

    // returns the index of an entry in a temperature field
    static Scalar pressGasIdx_(Scalar pressure, unsigned tempIdx)
    {
        Scalar pgMin = minGasPressure_(tempIdx);
        Scalar pgMax = maxGasPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - pgMin)/(pgMax - pgMin);
    }

    // returns the index of an entry in a density field
    static Scalar densityLiquidIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minLiquidDensity_(tempIdx);
        Scalar densityMax = maxLiquidDensity_(tempIdx);
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    // returns the index of an entry in a density field
    static Scalar densityGasIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minGasDensity_(tempIdx);
        Scalar densityMax = maxGasDensity_(tempIdx);
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    // returns the minimum tabulized liquid pressure at a given
    // temperature index
    static Scalar minLiquidPressure_(int tempIdx)
    {
        if (!useVaporPressure)
            return pressMin_;
        else
            return std::max<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1);
    }

    // returns the maximum tabulized liquid pressure at a given
    // temperature index
    static Scalar maxLiquidPressure_(int tempIdx)
    {
        if (!useVaporPressure)
            return pressMax_;
        else
            return std::max<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    // returns the minumum tabulized gas pressure at a given
    // temperature index
    static Scalar minGasPressure_(int tempIdx)
    {
        if (!useVaporPressure)
            return pressMin_;
        else
            return std::min<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1 );
    }

    // returns the maximum tabulized gas pressure at a given
    // temperature index
    static Scalar maxGasPressure_(int tempIdx)
    {
        if (!useVaporPressure)
            return pressMax_;
        else
            return std::min<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }


    // returns the minimum tabulized liquid density at a given
    // temperature index
    static Scalar minLiquidDensity_(int tempIdx)
    { return minLiquidDensity__[tempIdx]; }

    // returns the maximum tabulized liquid density at a given
    // temperature index
    static Scalar maxLiquidDensity_(int tempIdx)
    { return maxLiquidDensity__[tempIdx]; }

    // returns the minumum tabulized gas density at a given
    // temperature index
    static Scalar minGasDensity_(int tempIdx)
    { return minGasDensity__[tempIdx]; }

    // returns the maximum tabulized gas density at a given
    // temperature index
    static Scalar maxGasDensity_(int tempIdx)
    { return maxGasDensity__[tempIdx]; }


#ifndef NDEBUG
    // specifies whether the table was initialized
    static bool initialized_;
    // specifies whether some warning was printed
    static bool warningPrinted_;
#endif

    // 1D fields with the temperature as degree of freedom
    static Scalar *vaporPressure_;

    static Scalar *minLiquidDensity__;
    static Scalar *maxLiquidDensity__;

    static Scalar *minGasDensity__;
    static Scalar *maxGasDensity__;

    // 2D fields with the temperature and pressure as degrees of
    // freedom
    static Scalar *gasEnthalpy_;
    static Scalar *liquidEnthalpy_;

    static Scalar *gasHeatCapacity_;
    static Scalar *liquidHeatCapacity_;

    static Scalar *gasDensity_;
    static Scalar *liquidDensity_;

    static Scalar *gasViscosity_;
    static Scalar *liquidViscosity_;

    static Scalar *gasThermalConductivity_;
    static Scalar *liquidThermalConductivity_;

    // 2D fields with the temperature and density as degrees of
    // freedom
    static Scalar *gasPressure_;
    static Scalar *liquidPressure_;

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

#ifndef NDEBUG
template <class Scalar, class RawComponent, bool useVaporPressure>
bool TabulatedComponent<Scalar, RawComponent, useVaporPressure>::initialized_ = false;

template <class Scalar, class RawComponent, bool useVaporPressure>
bool TabulatedComponent<Scalar, RawComponent, useVaporPressure>::warningPrinted_ = false;
#endif

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


} // end namepace

#endif
