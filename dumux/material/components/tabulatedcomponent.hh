// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Tabulates all thermodynamic properties of a given
 *        untabulated chemical species.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 */
#ifndef DUMUX_TABULATED_COMPONENT_HH
#define DUMUX_TABULATED_COMPONENT_HH

#include <boost/math/special_functions/fpclassify.hpp>

namespace Dumux
{

/*!
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
    /*!
     * \brief Initialize the tables.
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
        gasInternalEnergy_ = new Scalar[nTemp_*nPress_];
        liquidInternalEnergy_ = new Scalar[nTemp_*nPress_];
        gasDensity_ = new Scalar[nTemp_*nPress_];
        liquidDensity_ = new Scalar[nTemp_*nPress_];
        gasViscosity_ = new Scalar[nTemp_*nPress_];
        liquidViscosity_ = new Scalar[nTemp_*nPress_];
        gasPressure_ = new Scalar[nTemp_*nDensity_];
        liquidPressure_ = new Scalar[nTemp_*nDensity_];

        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        Scalar NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // fill the arrays
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

                try { gasInternalEnergy_[i] = RawComponent::gasInternalEnergy(temperature, pressure); }
                catch (NumericalProblem) { gasInternalEnergy_[i] = NaN; };

                try { gasDensity_[i] = RawComponent::gasDensity(temperature, pressure); }
                catch (NumericalProblem) { gasDensity_[i] = NaN; };

                try { gasViscosity_[i] = RawComponent::gasViscosity(temperature, pressure); }
                catch (NumericalProblem) { gasViscosity_[i] = NaN; };
            };

            // calculate the minimum and maximum values for the gas
            // densities
            minGasDensity__[iT] = RawComponent::gasDensity(temperature, pgMin);
            maxGasDensity__[iT] = RawComponent::gasDensity(temperature, pgMax);

            // fill the temperature, density gas arrays
            for (unsigned iRho = 0; iRho < nDensity_; ++ iRho) {
                Scalar density =
                    iRho * (maxGasDensity__[iT] - minGasDensity__[iT])
                    /
                    (nDensity_ - 1)
                    +
                    minGasDensity__[iT];

                unsigned i = iT + iRho*nTemp_;

                try { gasPressure_[i] = RawComponent::gasPressure(temperature, density); }
                catch (NumericalProblem) { gasPressure_[i] = NaN; };
            };

            Scalar plMin = minLiquidPressure_(iT);
            Scalar plMax = maxLiquidPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (plMax - plMin)/(nPress_ - 1) + plMin;

                unsigned i = iT + iP*nTemp_;

                try { liquidEnthalpy_[i] = RawComponent::liquidEnthalpy(temperature, pressure); }
                catch (NumericalProblem) { liquidEnthalpy_[i] = NaN; };

                try { liquidInternalEnergy_[i] = RawComponent::liquidInternalEnergy(temperature, pressure); }
                catch (NumericalProblem) { liquidInternalEnergy_[i] = NaN; };

                try { liquidDensity_[i] = RawComponent::liquidDensity(temperature, pressure); }
                catch (NumericalProblem) { liquidDensity_[i] = NaN; };

                try { liquidViscosity_[i] = RawComponent::liquidViscosity(temperature, pressure); }
                catch (NumericalProblem) { liquidViscosity_[i] = NaN; };
            }

            // calculate the minimum and maximum values for the liquid
            // densities
            minLiquidDensity__[iT] = RawComponent::liquidDensity(temperature, plMin);
            maxLiquidDensity__[iT] = RawComponent::liquidDensity(temperature, plMax);

            // fill the temperature, density liquid arrays
            for (unsigned iRho = 0; iRho < nDensity_; ++ iRho) {
                Scalar density =
                    iRho * (maxLiquidDensity__[iT] - minLiquidDensity__[iT])
                    /
                    (nDensity_ - 1)
                    +
                    minLiquidDensity__[iT];

                unsigned i = iT + iRho*nTemp_;

                try { liquidPressure_[i] = RawComponent::liquidPressure(temperature, density); }
                catch (NumericalProblem) { liquidPressure_[i] = NaN; };
            };
        }
    }

    /*!
     * \brief A human readable name for the compoent.
     */
    static const char *name()
    { return RawComponent::name(); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static Scalar molarMass()
    { return RawComponent::molarMass(); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static Scalar criticalTemperature()
    { return RawComponent::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static Scalar criticalPressure()
    { return RawComponent::criticalPressure(); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static Scalar tripleTemperature()
    { return RawComponent::tripleTemperature(); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static Scalar triplePressure()
    { return RawComponent::triplePressure(); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
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
     * \brief Specific enthalpy of the gas [J/kg].
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasEnthalpy_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward gasEnthalpy("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::gasEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific internal energy of the liquid [J/kg].
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidEnthalpy_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward liquidEnthalpy("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::liquidEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific internal energy of the gas [J/kg].
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasInternalEnergy_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward gasInternalEnergy("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::gasInternalEnergy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific enthalpy of the liquid [J/kg].
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidInternalEnergy_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward liquidInternalEnergy("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::liquidInternalEnergy(temperature, pressure);
        }
        return result;
    }


    /*!
     * \brief The pressure of gas at a given density and temperature [Pa].
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateGasTRho_(liquidPressure_,
                                            temperature,
                                            density);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward gasPressure("<<temperature<<", "<<density<<")\n";
            return RawComponent::gasPressure(temperature,
                                             density);
        }
        return result;
    };

    /*!
     * \brief The pressure of liquid at a given density and temperature [Pa].
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateLiquidTRho_(liquidPressure_,
                                               temperature,
                                               density);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward liquidPressure("<<temperature<<", "<<density<<")\n";
            return RawComponent::liquidPressure(temperature,
                                                density);
        }
        return result;
    };

    /*!
     * \brief The density of gas at a given pressure and temperature
     *        [kg/m^3].
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasDensity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward gasDensity("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::gasDensity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The density of liquid at a given pressure and
     *        temperature [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidDensity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward liquidDensity("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::liquidDensity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of gas.
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateGasTP_(gasViscosity_,
                                          temperature,
                                          pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward gasViscosity("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::gasViscosity(temperature, pressure);
        }
        return result;
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of liquid.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidViscosity_,
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            if (verbose)
                std::cerr << "forward liquidViscosity("<<temperature<<", "<<pressure<<")\n";
            return RawComponent::liquidViscosity(temperature, pressure);
        }
        return result;
    };

private:
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

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = pressLiquidIdx_(p, iT);
        Scalar alphaP2 = pressLiquidIdx_(p, iT + 1);

        unsigned iP1 = std::min(nPress_ - 2, (unsigned) alphaP1);
        alphaP1 -= iP1;
        unsigned iP2 = std::min(nPress_ - 2, (unsigned) alphaP2);
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

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = pressGasIdx_(p, iT);
        Scalar alphaP2 = pressGasIdx_(p, iT + 1);

        unsigned iP1 = std::min(nPress_ - 2, (unsigned) alphaP1);
        alphaP1 -= iP1;
        unsigned iP2 = std::min(nPress_ - 2, (unsigned) alphaP2);
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
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            // std::cerr << __LINE__ << " T: " << T << "\n";
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = densityGasIdx_(rho, iT);
        Scalar alphaP2 = densityGasIdx_(rho, iT + 1);

        unsigned iP1 = std::min(nDensity_ - 2, (unsigned) alphaP1);
        alphaP1 -= iP1;
        unsigned iP2 = std::min(nDensity_ - 2, (unsigned) alphaP2);
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
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            // std::cerr << __LINE__ << " T: " << T << "\n";
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = densityLiquidIdx_(rho, iT);
        Scalar alphaP2 = densityLiquidIdx_(rho, iT + 1);

        unsigned iP1 = std::min(nDensity_ - 2, (unsigned) alphaP1);
        alphaP1 -= iP1;
        unsigned iP2 = std::min(nDensity_ - 2, (unsigned) alphaP2);
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
        Scalar pressMax = maxGasPressure_(tempIdx);
        return (nPress_ - 1)*(pressure - pressMin_)/(pressMax - pressMin_);
    }

    // returns the index of an entry in a density field
    static Scalar densityLiquidIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minLiquidDensity_(tempIdx);
        Scalar densityMax = maxLiquidDensity_(tempIdx);
        return (nDensity_ - 1)*(density - densityMin)/(densityMax - densityMin);
    }

    // returns the index of an entry in a density field
    static Scalar densityGasIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minGasDensity_(tempIdx);
        Scalar densityMax = maxGasDensity_(tempIdx);
        return (nDensity_ - 1)*(density - densityMin)/(densityMax - densityMin);
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

    static Scalar *gasInternalEnergy_;
    static Scalar *liquidInternalEnergy_;

    static Scalar *gasDensity_;
    static Scalar *liquidDensity_;

    static Scalar *gasViscosity_;
    static Scalar *liquidViscosity_;

    // 2D fields with the temperature and density as degrees of
    // freedom
    static Scalar *gasPressure_;
    static Scalar *liquidPressure_;

    // temperature, pressure and density ranges
    static Scalar   tempMin_;
    static Scalar   tempMax_;
    static unsigned nTemp_;

    static Scalar   pressMin_;
    static Scalar   pressMax_;
    static unsigned nPress_;

    static Scalar   densityMin_;
    static Scalar   densityMax_;
    static unsigned nDensity_;
};

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
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::gasInternalEnergy_;
template <class Scalar, class RawComponent, bool verbose>
Scalar* TabulatedComponent<Scalar, RawComponent, verbose>::liquidInternalEnergy_;
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
Scalar  TabulatedComponent<Scalar, RawComponent, verbose>::tempMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar  TabulatedComponent<Scalar, RawComponent, verbose>::tempMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nTemp_;
template <class Scalar, class RawComponent, bool verbose>
Scalar   TabulatedComponent<Scalar, RawComponent, verbose>::pressMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar   TabulatedComponent<Scalar, RawComponent, verbose>::pressMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nPress_;
template <class Scalar, class RawComponent, bool verbose>
Scalar   TabulatedComponent<Scalar, RawComponent, verbose>::densityMin_;
template <class Scalar, class RawComponent, bool verbose>
Scalar   TabulatedComponent<Scalar, RawComponent, verbose>::densityMax_;
template <class Scalar, class RawComponent, bool verbose>
unsigned TabulatedComponent<Scalar, RawComponent, verbose>::nDensity_;


} // end namepace

#endif
