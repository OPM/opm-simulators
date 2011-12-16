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
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
#ifndef DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH
#define DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH

#include <dumux/common/valgrind.hh>

#include <tr1/array>

namespace Dumux
{
/*!
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
template <class Scalar, class FluidState>
class PressureOverlayFluidState
{
public:
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    /*!
     * \brief Constructor
     *
     * The overlay fluid state copies the pressures from the argument,
     * so it initially behaves exactly like the underlying fluid
     * state.
     */
    PressureOverlayFluidState(const FluidState &fs)
        : fs_(&fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
    }

    // copy constructor
    PressureOverlayFluidState(const PressureOverlayFluidState &fs)
        : fs_(fs.fs_)
        , pressure_(fs.pressure_)
    {
    }

    // assignment operator
    PressureOverlayFluidState &operator=(const PressureOverlayFluidState &fs)
    {
        fs_ = fs.fs_;
        pressure_ = fs.pressure_;
        return *this;
    }

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the saturation of a phase []
     */
    Scalar saturation(int phaseIdx) const
    { return fs_->saturation(phaseIdx); }

    /*!
     * \brief The mole fraction of a component in a phase []
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief The mass fraction of a component in a phase []
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return fs_->massFraction(phaseIdx, compIdx); }

    /*!
     * \brief The average molar mass of a fluid phase [kg/mol]
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fs_->averageMolarMass(phaseIdx); }

    /*!
     * \brief The concentration of a component in a phase [mol/m^3]
     *
     * This is usually just called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return fs_->molarity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fs_->fugacity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity coefficient of a component in a phase [Pa]
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fs_->fugacityCoefficient(phaseIdx, compIdx); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(int phaseIdx) const
    { return fs_->molarVolume(phaseIdx); }

    /*!
     * \brief The mass density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { return fs_->density(phaseIdx); }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(int phaseIdx) const
    { return fs_->molarDensity(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase [K]
     */
    Scalar temperature(int phaseIdx) const
    { return fs_->temperature(phaseIdx); }

    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { return fs_->enthalpy(phaseIdx); }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fs_->internalEnergy(phaseIdx); }

    /*!
     * \brief The dynamic viscosity of a fluid phase [Pa s]
     */
    Scalar viscosity(int phaseIdx) const
    { return fs_->viscosity(phaseIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the pressure [Pa] of a fluid phase
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

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
        Valgrind::CheckDefined(pressure_);
    }

protected:
    const FluidState *fs_;
    std::tr1::array<Scalar, numPhases> pressure_;
};

} // end namepace Dumux

#endif
