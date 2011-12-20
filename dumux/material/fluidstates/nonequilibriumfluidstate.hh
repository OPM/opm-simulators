// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
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
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
#ifndef DUMUX_GENERIC_FLUID_STATE_HH
#define DUMUX_GENERIC_FLUID_STATE_HH

#include <dumux/common/valgrind.hh>

#include <cmath>
#include <algorithm>

namespace Dumux
{
/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
template <class Scalar, class FluidSystem>
class NonEquilibriumFluidState
{
public:
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    NonEquilibriumFluidState()
    { 
        // set the composition to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)  {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                moleFraction_[phaseIdx][compIdx] = 0;
            
            averageMolarMass_[phaseIdx] = 0;
            sumMoleFractions_[phaseIdx] = 0;
        }
        
        // make everything undefined so that valgrind will complain
        Valgrind::SetUndefined(*this);
    }

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the saturation of a phase []
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief The mole fraction of a component in a phase []
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return moleFraction_[phaseIdx][compIdx]; }


    /*!
     * \brief The mass fraction of a component in a phase []
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return
            sumMoleFractions_[phaseIdx]
            * moleFraction_[phaseIdx][compIdx]
            * FluidSystem::molarMass(compIdx)
            / std::max(1e-10, averageMolarMass_[phaseIdx]);
    }

    /*!
     * \brief The mean molar mass of a fluid phase [kg/mol]
     *
     * The average mass is the mean molar mass of a molecule of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \[ \bar M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa \]
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return averageMolarMass_[phaseIdx]; }

    /*!
     * \brief The concentration of a component in a phase [mol/m^3]
     *
     * This quantity is often called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return pressure_[phaseIdx]*fugacityCoefficient_[phaseIdx][compIdx]*moleFraction_[phaseIdx][compIdx]; }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * \brief The mass density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(int phaseIdx) const
    { return density_[phaseIdx]/averageMolarMass(phaseIdx); }


    /*!
     * \brief The temperature of a fluid phase [K]
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; };
    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return enthalpy_[phaseIdx] - pressure(phaseIdx)/density(phaseIdx); }

    /*!
     * \brief The dynamic viscosity of a fluid phase [Pa s]
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; };

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature [K] of a fluid phase
     */
    void setTemperature(int phaseIdx, Scalar value)
    { temperature_[phaseIdx] = value; }

    /*!
     * \brief Set the fluid pressure of a phase [Pa]
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

    /*!
     * \brief Set the saturation of a phase []
     */
    void setSaturation(int phaseIdx, Scalar value)
    { saturation_[phaseIdx] = value; }

    /*!
     * \brief Set the mole fraction of a component in a phase []
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    { 
        Valgrind::CheckDefined(value);

        if (std::isfinite(averageMolarMass_[phaseIdx])) {
            Scalar delta = value - moleFraction_[phaseIdx][compIdx];
            
            moleFraction_[phaseIdx][compIdx] = value;
            
            sumMoleFractions_[phaseIdx] += delta;
            averageMolarMass_[phaseIdx] += delta*FluidSystem::molarMass(compIdx);
        }
        else { 
            moleFraction_[phaseIdx][compIdx] = value;
            
            // re-calculate the mean molar mass
            sumMoleFractions_[phaseIdx] = 0.0;
            averageMolarMass_[phaseIdx] = 0.0;
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compJIdx];
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compJIdx]*FluidSystem::molarMass(compJIdx);
            }
        }
        
        Valgrind::SetDefined(sumMoleFractions_[phaseIdx]);
        Valgrind::SetDefined(averageMolarMass_[phaseIdx]);
    }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(int phaseIdx, int compIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Set the density of a phase [kg / m^3]
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the specific enthalpy of a phase [J/m^3]
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { enthalpy_[phaseIdx] = value; }

    /*!
     * \brief Set the dynamic viscosity of a phase [Pa s]
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            averageMolarMass_[phaseIdx] = 0;
            sumMoleFractions_[phaseIdx] = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFraction_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
                fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
                averageMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
                sumMoleFractions_[phaseIdx] += moleFraction_[phaseIdx][compIdx];
            }

            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            density_[phaseIdx] = fs.density(phaseIdx);
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
            temperature_[phaseIdx] = fs.temperature(phaseIdx);
        }
    };

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
#if HAVE_VALGRIND && ! defined NDEBUG
        for (int i = 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                Valgrind::CheckDefined(fugacityCoefficient_[i][j]);
                Valgrind::CheckDefined(moleFraction_[i][j]);
            }
            Valgrind::CheckDefined(averageMolarMass_[i]);
            Valgrind::CheckDefined(pressure_[i]);
            Valgrind::CheckDefined(saturation_[i]);
            Valgrind::CheckDefined(density_[i]);
            Valgrind::CheckDefined(temperature_[i]);
            //Valgrind::CheckDefined(internalEnergy_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }
#endif // HAVE_VALGRIND
    }

protected:
    Scalar moleFraction_[numPhases][numComponents];
    Scalar fugacityCoefficient_[numPhases][numComponents];

    Scalar averageMolarMass_[numPhases];
    Scalar sumMoleFractions_[numPhases];
    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar density_[numPhases];
    Scalar enthalpy_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_[numPhases];
};

} // end namepace Dumux

#endif
