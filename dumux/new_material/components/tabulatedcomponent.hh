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
 */
#ifndef DUMUX_TABULATED_COMPONENT_HH
#define DUMUX_TABULATED_COMPONENT_HH

#include <boost/math/special_functions/fpclassify.hpp>

namespace Dune
{

/*!
 * \brief Abstract base class of a pure chemical species.
 */
template <class Scalar, class RawComponent>
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
                  
        // allocate the arrays
        vaporPressure_ = new Scalar[nTemp_];

        gasEnthalpy_ = new Scalar[nTemp_*nPress_];
        liquidEnthalpy_ = new Scalar[nTemp_*nPress_];
        gasDensity_ = new Scalar[nTemp_*nPress_];
        liquidDensity_ = new Scalar[nTemp_*nPress_];
        gasViscosity_ = new Scalar[nTemp_*nPress_];
        liquidViscosity_ = new Scalar[nTemp_*nPress_];

        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        Scalar NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // fill the arrays
        for (unsigned iT = 0; iT < nTemp_; ++ iT) {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            try { vaporPressure_[iT] = RawComponent::vaporPressure(temperature); }
            catch (NumericalProblem e) { vaporPressure_[iT] = NaN; };
            
            Scalar pgMax = maxGasPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (pgMax - pressMin_)/(nPress_ - 1) + pressMin_;

                unsigned i = iT + iP*nTemp_;

                try { gasEnthalpy_[i] = RawComponent::gasEnthalpy(temperature, pressure); }
                catch (NumericalProblem) { gasEnthalpy_[i] = NaN; };

                try { gasDensity_[i] = RawComponent::gasDensity(temperature, pressure); }
                catch (NumericalProblem) { gasDensity_[i] = NaN; };

                try { gasViscosity_[i] = RawComponent::gasViscosity(temperature, pressure); }
                catch (NumericalProblem) { gasViscosity_[i] = NaN; };
            };
            
            Scalar plMin = minLiquidPressure_(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP) {
                Scalar pressure = iP * (pressMax_ - plMin)/(nPress_ - 1) + plMin;

                unsigned i = iT + iP*nTemp_;
                
                try { liquidEnthalpy_[i] = RawComponent::liquidEnthalpy(temperature, pressure); }
                catch (NumericalProblem) { liquidEnthalpy_[i] = NaN; };
                
                try { liquidDensity_[i] = RawComponent::liquidDensity(temperature, pressure); }
                catch (NumericalProblem) { liquidDensity_[i] = NaN; };

                try { liquidViscosity_[i] = RawComponent::liquidViscosity(temperature, pressure); }
                catch (NumericalProblem) { liquidViscosity_[i] = NaN; };
            }
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
            return RawComponent::gasEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific enthalpy of the liquid [J/kg].
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateLiquidTP_(liquidEnthalpy_, 
                                             temperature,
                                             pressure);
        if (std::isnan(result)) {
            return RawComponent::liquidEnthalpy(temperature, pressure);
        }
        return result;
    }

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
            // std::cout << __LINE__ << " T: " << T << "\n";
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        
        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = pressLiquidIdx_(p, iT);
        if (alphaP1 < 0 || alphaP1 >= nPress_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        Scalar alphaP2 = pressLiquidIdx_(p, iT + 1);
        if (alphaP2 < 0 || alphaP2 >= nPress_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        
        unsigned iP1 = (unsigned) alphaP1;
        alphaP1 -= iP1;
        unsigned iP2 = (unsigned) alphaP2;
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
            // std::cout << __LINE__ << " T: " << T << "\n";
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        
        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        Scalar alphaP1 = pressGasIdx_(p, iT);
        if (alphaP1 < 0 || alphaP1 >= nPress_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        Scalar alphaP2 = pressGasIdx_(p, iT + 1);
        if (alphaP2 < 0 || alphaP2 >= nPress_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        unsigned iP1 = (unsigned) alphaP1;
        alphaP1 -= iP1;
        unsigned iP2 = (unsigned) alphaP2;
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
        Scalar pressMin = minLiquidPressure_(tempIdx);
        return (nPress_ - 1)*(pressure - pressMin)/(pressMax_ - pressMin);
    }

    // returns the index of an entry in a temperature field
    static Scalar pressGasIdx_(Scalar pressure, unsigned tempIdx)
    {
        Scalar pressMax = maxGasPressure_(tempIdx);
        return (nPress_ - 1)*(pressure - pressMin_)/(pressMax - pressMin_);
    }
    
    // returns the minimum tabulized liquid pressure at a given
    // temperature index
    static Scalar minLiquidPressure_(int tempIdx)
    { return std::max(pressMin_, vaporPressure_[tempIdx]/1.2); }

    // returns the maximum tabulized gas pressure at a given
    // temperature index
    static Scalar maxGasPressure_(int tempIdx)
    { return std::min(pressMax_, vaporPressure_[tempIdx]*1.3); }

    static Scalar *vaporPressure_;

    static Scalar *gasEnthalpy_;
    static Scalar *liquidEnthalpy_;

    static Scalar *gasDensity_;
    static Scalar *liquidDensity_;

    static Scalar *gasViscosity_;
    static Scalar *liquidViscosity_;

    static Scalar   tempMin_;
    static Scalar   tempMax_;
    static unsigned nTemp_;

    static Scalar   pressMin_;
    static Scalar   pressMax_;
    static unsigned nPress_;
};

template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::vaporPressure_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::gasEnthalpy_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::liquidEnthalpy_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::gasDensity_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::liquidDensity_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::gasViscosity_;
template <class Scalar, class RawComponent>
Scalar* TabulatedComponent<Scalar, RawComponent>::liquidViscosity_;
template <class Scalar, class RawComponent>
Scalar  TabulatedComponent<Scalar, RawComponent>::tempMin_;
template <class Scalar, class RawComponent>
Scalar  TabulatedComponent<Scalar, RawComponent>::tempMax_;
template <class Scalar, class RawComponent>
unsigned TabulatedComponent<Scalar, RawComponent>::nTemp_;
template <class Scalar, class RawComponent>
Scalar   TabulatedComponent<Scalar, RawComponent>::pressMin_;
template <class Scalar, class RawComponent>
Scalar   TabulatedComponent<Scalar, RawComponent>::pressMax_;
template <class Scalar, class RawComponent>
unsigned TabulatedComponent<Scalar, RawComponent>::nPress_;


} // end namepace

#endif
