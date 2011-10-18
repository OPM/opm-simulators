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
 * \brief A GenericFluidState which keeps track on which variables changed.
 */
#ifndef DUMUX_TRACKING_GENERIC_FLUID_STATE_HH
#define DUMUX_TRACKING_GENERIC_FLUID_STATE_HH

#include "genericfluidstate.hh"

#include <dumux/common/valgrind.hh>

#include <cmath>
#include <algorithm>

namespace Dumux
{
/*!
 * \brief A GenericFluidState which keeps track on which variables changed.
 */
template <class Scalar, class StaticParameters>
class TrackingGenericFluidState
{
    typedef Dumux::GenericFluidState<Scalar, StaticParameters> GenericFluidState;
public:
    enum { numComponents = StaticParameters::numComponents };
    enum { numPhases = StaticParameters::numPhases };

    TrackingGenericFluidState()
    {
        for (int i = 0; i < numPhases; ++i)
            changed_[i].raw = ~((unsigned) 0x00);
    }

    /*****************************************************
     * Generic access to fluid properties (No assumptions 
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the saturation of a phase []
     */
    Scalar saturation(int phaseIdx) const
    { return fs_.saturation(phaseIdx); }

    /*!
     * \brief The mole fraction of a component in a phase []
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    { return fs_.moleFrac(phaseIdx, compIdx); }

    /*!
     * \brief The mass fraction of a component in a phase []
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    { return fs_.massFrac(phaseIdx, compIdx); }

    /*!
     * \brief The sum of all component mole fractions in a phase []
     *
     * We define this to be the same as the sum of all mass fractions.
     */
    Scalar sumMoleFrac(int phaseIdx) const
    { return fs_.sumMoleFrac(phaseIdx); }

    /*!
     * \brief The sum of all component mass fractions in a phase []
     *
     * We define this to be the same as the sum of all mole fractions.
     */
    Scalar sumMassFrac(int phaseIdx) const
    { return fs_.sumMassFrac(phaseIdx); }

    /*!
     * \brief The mean molar mass of a fluid phase [kg/mol]
     */
    Scalar meanMolarMass(int phaseIdx) const
    { return fs_.meanMolarMass(phaseIdx); }

    /*!
     * \brief The molar concentration of a component in a phase [mol/m^3]
     *
     * This is often just called "concentration", but there are many
     * other (though less common) measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return fs_.molarity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fs_.fugacity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoeff(int phaseIdx, int compIdx) const
    { return fs_.fugacityCoeff(phaseIdx, compIdx); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(int phaseIdx) const
    { return fs_.molarVolume(phaseIdx); }

    /*!
     * \brief The mass density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { return fs_.density(phaseIdx); }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(int phaseIdx) const
    { return fs_.molarDensity(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase [K]
     */
    Scalar temperature(int phaseIdx) const
    { return fs_.temperature(phaseIdx); }
    
    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    Scalar pressure(int phaseIdx) const
    { return fs_.pressure(phaseIdx); }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { return fs_.enthalpy(phaseIdx); }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fs_.internalEnergy(phaseIdx); }

    /*!
     * \brief The dynamic viscosity of a fluid phase [Pa s]
     */
    Scalar viscosity(int phaseIdx) const
    { return fs_.viscosity(phaseIdx); };

    /*****************************************************
     * Setter methods. Note that these are not part of the 
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature [K] of a fluid phase
     */
    void setTemperature(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.temperature = 1;
        fs_.setTemperature(phaseIdx, value);
    };

    /*!
     * \brief Set the fluid pressure of a phase [Pa]
     */
    void setPressure(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.pressure = 1;
        fs_.setPressure(phaseIdx, value);
    };

    /*!
     * \brief Set the saturation of a phase []
     */
    void setSaturation(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.saturation = 1;
        fs_.setSaturation(phaseIdx, value);
    };

    /*!
     * \brief Set the mole fraction of a component in a phase []
     */
    void setMoleFrac(int phaseIdx, int compIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.moleFrac |= 1 << compIdx;
        fs_.setMoleFrac(phaseIdx, compIdx, value);
    };


    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoeff(int phaseIdx, int compIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.fugacityCoeff |= 1 << compIdx;
        fs_.setFugacityCoeff(phaseIdx, compIdx, value);
    }

    /*!
     * \brief Set the molar volume of a phase [m^3/mol]
     */
    void setMolarVolume(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.molarVolume = 1;
        fs_.setMolarVolume(phaseIdx, value);
    }

    /*!
     * \brief Set the specific enthalpy of a phase [J/m^3]
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.enthalpy = 1;
        fs_.setEnthalpy(phaseIdx, value);
    }   

    /*!
     * \brief Set the dynamic viscosity of a phase [Pa s]
     */
    void setViscosity(int phaseIdx, Scalar value)
    { 
        changed_[phaseIdx].fields.viscosity = 1;
        fs_.setViscosity(phaseIdx, value);
    }   

    /*!
     * \brief Calculatate the mean molar mass of a phase given that
     *        all mole fractions have been set
     */
    void updateMeanMolarMass(int phaseIdx)
    { fs_.updateMeanMolarMass(phaseIdx); }
    

    /*!
     * \brief Reset the dirty tracking for a phase.
     */
    void setPhaseClean(int phaseIdx)
    {
        // set all tracking bits for the phase to zero
        changed_[phaseIdx].raw = 0;
    }

    /*!
     * \brief Reset the dirty tracking of all phases.
     */
    void setClean()
    {
        for (int i = 0; i < numPhases; ++i)
            changed_[i].raw = 0;
    }

    /*!
     * \brief Return if all phases are clean.
     */
    bool isClean() const
    { 
        for (int i = 0; i < numPhases; ++i)
            if (changed_[i].raw != 0)
                return false;
        return true;
    };

    /*!
     * \brief Find out whether the values for a phase been since the
     *        last call to cleanPhase().
     */
    bool phaseClean(int phaseIdx) const
    { return changed_[phaseIdx].raw == 0; }


    /*!
     * \brief Find out whether a specific phase's temperature has been
     *        changed.
     */
    bool temperatureClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.temperature == 0; }

    /*!
     * \brief Find out whether a specific phase's saturation has been
     *        changed.
     */
    bool saturationClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.saturation == 0; }

    /*!
     * \brief Find out whether a specific phase's pressure has been
     *        changed.
     */
    bool pressureClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.pressure == 0; }

    /*!
     * \brief Find out whether a specific phase's enthalpy has been
     *        changed.
     */
    bool enthalpyClean(int phaseIdx) const
    { return changed_[phaseIdx].enthalpy == 0; }

    /*!
     * \brief Find out whether a specific phase's temperature has been
     *        changed.
     */
    bool molarVolumeClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.molarVolume == 0; }

    /*!
     * \brief Find out whether some of the mole fractions of a phase
     *        have been changed since the last call to
     *        setCleanPhase().
     */
    bool moleFracsClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.moleFrac == 0; }

    /*!
     * \brief Find out whether a specific mole fraction of a component
     *        in a phase has been changed since the last call to
     *        setCleanPhase().
     */
    bool moleFracsClean(int phaseIdx, int compIdx) const
    { return (changed_[phaseIdx].fields.moleFrac  & (1 << compIdx)) == 0; }

    /*!
     * \brief Find out whether at least one fugacity of a component in
     *        a phase has been changed since the last call to
     *        setCleanPhase().
     */
    bool fugacityCoeffsClean(int phaseIdx) const
    { return changed_[phaseIdx].fields.fugacitCoeffs == 0; }

    /*!
     * \brief Find out whether a specific fugacity of a component in a
     *        phase has been changed since the last call to
     *        setCleanPhase().
     */
    bool fugacityCoeffClean(int phaseIdx, int compIdx) const
    { return (changed_[phaseIdx].fields.fugacityCoeffs  & (1 << compIdx)) == 0; }
    

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
        Valgrind::CheckDefined(changed_);
        fs_.checkDefined();
#endif // HAVE_VALGRIND
    }

protected:
    static_assert(sizeof(unsigned)*8 >= 5 + 2*numComponents,
                  "Too many components to use the bitfield trick");
    union {
        unsigned raw;
        struct {
            unsigned char moleFrac : numComponents;
            unsigned char fugacityCoeff : numComponents;

            unsigned char molarVolume : 1;
            unsigned char temperature : 1;
            unsigned char saturation : 1;
            unsigned char pressure : 1;
            unsigned char enthalpy : 1;
            unsigned char viscosity : 1;
        } __attribute__((__packed__)) fields;
    } changed_[numPhases];
    GenericFluidState fs_;
};

} // end namepace Dumux

#endif
