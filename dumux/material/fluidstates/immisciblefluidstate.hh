/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_IMMISCIBLE_FLUID_STATE_HH
#define DUMUX_IMMISCIBLE_FLUID_STATE_HH
#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
template <class Scalar, class StaticParameters>
class ImmiscibleFluidState
{
public:
    enum { numComponents = StaticParameters::numComponents };
    enum { numPhases = StaticParameters::numPhases };
    static_assert(numPhases == numComponents,
                  "The number of phases must be equal to the number of "
                  "components if immiscibility is assumed!");

    ImmiscibleFluidState()
    { Valgrind::SetUndefined(*this); }

    template <class FluidState>
    ImmiscibleFluidState(FluidState &fs)
    { assign(fs); }

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
    Scalar moleFrac(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * \brief The mass fraction of a component in a phase []
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?1.0:0.0; }

    /*!
     * \brief The sum of all component mole fractions in a phase []
     *
     * We define this to be the same as the sum of all mass fractions.
     */
    Scalar sumMoleFrac(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief The sum of all component mass fractions in a phase []
     *
     * We define this to be the same as the sum of all mole fractions.
     */
    Scalar sumMassFrac(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief The mean molar mass of a fluid phase [kg/mol]
     */
    Scalar meanMolarMass(int phaseIdx) const
    { return StaticParameters::molarMass(/*compIdx=*/phaseIdx); }

    /*!
     * \brief The concentration of a component in a phase [mol/m^3]
     *
     * This is often called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFrac(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     *
     * To avoid numerical issues with code that assumes miscibility,
     * we return a fugacity of 0 for components which do not mix with
     * the specified phase. (Actually it undefined, but for finite
     * fugacity coefficients, the only way to get components
     * completely out of a phase is 0 to feed it zero fugacity.)
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { 
        if (phaseIdx == compIdx)
            return fugacityCoeff(phaseIdx, compIdx)*moleFrac(phaseIdx, compIdx)*pressure(phaseIdx);
        else
            return 0.0;
    };
    
    /*!
     * \brief The fugacity coefficient of a component in a phase [Pa]
     *
     * Since we assume immiscibility, the fugacity coefficients for
     * the components which are not miscible with the phase is
     * infinite. Beware that this will very likely break your code if
     * you don't keep that in mind.
     */
    Scalar fugacityCoeff(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return fugacityCoeff_[phaseIdx];
        else
            return std::numeric_limits<Scalar>::infinity();
    }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(int phaseIdx) const
    { return molarVolume_[phaseIdx]; }

    /*!
     * \brief The mass density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { return molarDensity(phaseIdx)*meanMolarMass(phaseIdx); }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(int phaseIdx) const
    { return 1/molarVolume(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase [K]
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; }
    
    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return enthalpy(phaseIdx) - pressure(phaseIdx)/density(phaseIdx); }

    /*!
     * \brief The dynamic viscosity of a fluid phase [Pa s]
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }


    /*****************************************************
     * Access to fluid properties which only make sense
     * if assuming thermodynamic equilibrium
     *****************************************************/

    /*!
     * \brief The temperature within the domain [K]
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief The fugacity of a component
     *
     * This assumes chemical equilibrium.
     */
    Scalar fugacity(int compIdx) const
    { return fugacity(0, compIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the 
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    
    /*!
     * \brief Retrieve all parameters from an arbitrary fluid 
     *        state.
     *
     * \note If the other fluid state object is inconsistent with the
     *       thermodynamic equilibrium, the result of this method is
     *       undefined.
     */
    template <class FluidState>
    void assign(const FluidState &fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fugacityCoeff_[phaseIdx] = fs.fugacityCoeff(phaseIdx, phaseIdx);
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
            molarVolume_[phaseIdx] = fs.molarVolume(phaseIdx);
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
            viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
        }
        temperature_ = fs.temperature(0);
    };

    /*!
     * \brief Set the temperature [K] of a fluid phase
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }   

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
     * \brief Set the fugacity coefficient of a component in a phase []
     *
     * The phase is always index the same as the component, since we
     * assume immiscibility.
     */
    void setFugacityCoeff(int compIdx, Scalar value)
    { fugacityCoeff_[compIdx] = value; }   

    /*!
     * \brief Set the molar volume of a phase [m^3/mol]
     */
    void setMolarVolume(int phaseIdx, Scalar value)
    { molarVolume_[phaseIdx] = value; }   

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
            //for (int j = 0; j < numComponents; ++j) {
            //    Valgrind::CheckDefined(fugacityCoeff_[i][j]);
            //}
            Valgrind::CheckDefined(pressure_[i]);
            Valgrind::CheckDefined(saturation_[i]);
            Valgrind::CheckDefined(molarVolume_[i]);
            //Valgrind::CheckDefined(enthalpy_[i]);
            Valgrind::CheckDefined(viscosity_[i]);
        }

        Valgrind::CheckDefined(temperature_);
#endif // HAVE_VALGRIND
    }

protected:
    Scalar fugacityCoeff_[numComponents];

    Scalar pressure_[numPhases];
    Scalar saturation_[numPhases];
    Scalar molarVolume_[numPhases];
    Scalar enthalpy_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar temperature_;
};

} // end namepace Dumux

#endif
