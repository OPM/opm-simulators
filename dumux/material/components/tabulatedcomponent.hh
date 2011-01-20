// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \ingroup Components
 * \brief Tabulates all thermodynamic properties of a given
 *        untabulated chemical species.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 */
#ifndef DUMUX_TABULATED_COMPONENT_HH
#define DUMUX_TABULATED_COMPONENT_HH

#include <dumux/common/exceptions.hh>

#include <boost/math/special_functions/fpclassify.hpp>

#include <iostream>

namespace Dumux
{

/*!
 * \ingroup Components
 *
 * \brief  Tabulates all thermodynamic properties of a given
 *        untabulated chemical species.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 *
 * \tparam Scalar  The type used for scalar values
 * \tparam Scalar  The component which ought to be tabulated
 * \tparam verbose If set to true, a warning will be printed each time
 *                 a request can not be fulfilled from the tabulated
 *                 arrays. This is quite useful for debugging
 *                 purposes.
 */
template <class Scalar, class RawComponent, bool verbose=true>
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
        gasDensity_ = new Scalar[nTemp_*nPress_];
        liquidDensity_ = new Scalar[nTemp_*nPress_];
        gasViscosity_ = new Scalar[nTemp_*nPress_];
        liquidViscosity_ = new Scalar[nTemp_*nPress_];
        gasPressure_ = new Scalar[nTemp_*nDensity_];
        liquidPressure_ = new Scalar[nTemp_*nDensity_];

        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        Scalar NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // fill the temperature-pressure arrays
        for (unsigned iT = 0; iT < nTemp_; ++ iT) {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            try { vaporPressure_[iT] = RawComponent::vaporPressure(temperature); }
            catch (NumericalProblem e) { vaporPressure_[iT] = NaN; };

            Scalar pgMax = maxGasPressure_(iT);
            Scalar pgMin = minGasPressure_(iT);

            // fill the temperature, pressure gas arrays
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (pgMax - pgMin)/(nPress_ - 1) + pgMin;

                unsigned i = iT + iP*nTemp_;

                try { gasEnthalpy_[i] = RawComponent::gasEnthalpy(temperature, pressure); }
                catch (NumericalProblem) { gasEnthalpy_[i] = NaN; };

                try { gasDensity_[i] = RawComponent::gasDensity(temperature, pressure); }
                catch (NumericalProblem) { gasDensity_[i] = NaN; };

                try { gasViscosity_[i] = RawComponent::gasViscosity(temperature, pressure); }
                catch (NumericalProblem) { gasViscosity_[i] = NaN; };
            };

            Scalar plMin = minLiquidPressure_(iT);
            Scalar plMax = maxLiquidPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (plMax - plMin)/(nPress_ - 1) + plMin;

                unsigned i = iT + iP*nTemp_;

                try { liquidEnthalpy_[i] = RawComponent::liquidEnthalpy(temperature, pressure); }
                catch (NumericalProblem) { liquidEnthalpy_[i] = NaN; };

                try { liquidDensity_[i] = RawComponent::liquidDensity(temperature, pressure); }
                catch (NumericalProblem) { liquidDensity_[i] = NaN; };

                try { liquidViscosity_[i] = RawComponent::liquidViscosity(temperature, pressure); }
                catch (NumericalProblem) { liquidViscosity_[i] = NaN; };
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
    static Scalar vaporPressure(Scalar T)
    {
        Scalar result = interpolateT_(vaporPressure_, T);
        if (std::isnan(result)) {
            return RawComponent::vaporPressure(T);
        }
        return result;
    };

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
    };

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
    };

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
    };

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
    };

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
        else if (verbose)
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

        unsigned iT = std::max<long long>(0, std::min<long long>(nTemp_ - 2, (long long) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = pressLiquidIdx_(p, iT);
        Scalar alphaP2 = pressLiquidIdx_(p, iT + 1);

        unsigned iP1 = std::max<long long>(0, std::min<long long>(nPress_ - 2, (long long) alphaP1));
        unsigned iP2 = std::max<long long>(0, std::min<long long>(nPress_ - 2, (long long) alphaP2));
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
    static Scalar interpolateGasTP_(const Scalar *values, Scalar T, Scalar p)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            // std::cerr << __LINE__ << " T: " << T << "\n";
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iT = std::max<long long>(0, std::min<long long>(nTemp_ - 2, (long long) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = pressGasIdx_(p, iT);
        Scalar alphaP2 = pressGasIdx_(p, iT + 1);
        unsigned iP1 = std::max<long long>(0, std::min<long long>(nPress_ - 2, (long long) alphaP1));
        unsigned iP2 = std::max<long long>(0, std::min<long long>(nPress_ - 2, (long long) alphaP2));
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
    static Scalar interpolateGasTRho_(const Scalar *values, Scalar T, Scalar rho)
    {
        Scalar alphaT = tempIdx_(T);
        unsigned iT = std::max<long long>(0, std::min<long long>(nTemp_ - 2, (long long) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = densityGasIdx_(rho, iT);
        Scalar alphaP2 = densityGasIdx_(rho, iT + 1);
        unsigned iP1 = std::max<long long>(0, std::min<long long>(nDensity_ - 2, (long long) alphaP1));
        unsigned iP2 = std::max<long long>(0, std::min<long long>(nDensity_ - 2, (long long) alphaP2));
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
        unsigned iT = std::max<long long>(0, std::min<long long>(nTemp_ - 2, (long long) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = densityLiquidIdx_(rho, iT);
        Scalar alphaP2 = densityLiquidIdx_(rho, iT + 1);
        unsigned iP1 = std::max<long long>(0, std::min<long long>(nDensity_ - 2, (long long) alphaP1));
        unsigned iP2 = std::max<long long>(0, std::min<long long>(nDensity_ - 2, (long long) alphaP2));
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
    { return std::max<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1); }

    // returns the maximum tabulized liquid pressure at a given
    // temperature index
    static Scalar maxLiquidPressure_(int tempIdx)
    { return std::max<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1); }

    // returns the minumum tabulized gas pressure at a given
    // temperature index
    static Scalar minGasPressure_(int tempIdx)
    { return std::min<Scalar>(pressMin_, vaporPressure_[tempIdx] / 1.1 ); }

    // returns the maximum tabulized gas pressure at a given
    // temperature index
    static Scalar maxGasPressure_(int tempIdx)
    { return std::min<Scalar>(pressMax_, vaporPressure_[tempIdx] * 1.1); }


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

    static Scalar *gasDensity_;
    static Scalar *liquidDensity_;

    static Scalar *gasViscosity_;
    static Scalar *liquidViscosity_;

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
template <class Scalar, class RawComponent, bool verbose>
bool TabulatedComponent<Scalar, RawComponent, verbose>::initialized_ = false;

template <class Scalar, class RawComponent, bool verbose>
bool TabulatedComponent<Scalar, RawComponent, verbose>::warningPrinted_ = false;
#endif

template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::vaporPressure_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::minLiquidDensity__;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::maxLiquidDensity__;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::minGasDensity__;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::maxGasDensity__;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::gasEnthalpy_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::liquidEnthalpy_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::gasDensity_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::liquidDensity_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::gasViscosity_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::liquidViscosity_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::gasPressure_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::liquidPressure_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::tempMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::tempMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nTemp_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::pressMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::pressMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nPress_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::densityMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar TabulatedComponent<Scalar, RawComponent, verbose>::densityMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nDensity_;


} // end namepace

#endif
