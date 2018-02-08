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
 * \copydoc Opm::TemperatureOverlayFluidState
 */
#ifndef OPM_TEMPERATURE_OVERLAY_FLUID_STATE_HPP
#define OPM_TEMPERATURE_OVERLAY_FLUID_STATE_HPP

#include <opm/material/common/Valgrind.hpp>

#include <utility>

namespace Opm {

/*!
 * \brief This is a fluid state which allows to set the fluid
 *        temperatures and takes all other quantities from an other
 *        fluid state.
 */
template <class FluidState>
class TemperatureOverlayFluidState
{
public:
    typedef typename FluidState::Scalar Scalar;

    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    /*!
     * \brief Constructor
     *
     * The overlay fluid state copies the temperature from the first
     * fluid phase of the argument, so it initially behaves exactly
     * like the underlying fluid state, provided that the underlying
     * fluid state is in thermal equilibrium.
     */
    TemperatureOverlayFluidState(const FluidState& fs)
        : fs_(&fs)
    {
        temperature_ = fs.temperature(/*phaseIdx=*/0);
    }

    TemperatureOverlayFluidState(Scalar T, const FluidState& fs)
        : temperature_(T), fs_(&fs)
    { }

    // copy constructor
    TemperatureOverlayFluidState(const TemperatureOverlayFluidState& fs)
        : fs_(fs.fs_)
        , temperature_(fs.temperature_)
    {}

    // assignment operator
    TemperatureOverlayFluidState& operator=(const TemperatureOverlayFluidState& fs)
    {
        fs_ = fs.fs_;
        temperature_ = fs.temperature_;
        return *this;
    }

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the saturation of a phase []
     */
    auto saturation(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().saturation(phaseIdx))
    { return fs_->saturation(phaseIdx); }

    /*!
     * \brief The mole fraction of a component in a phase []
     */
    auto moleFraction(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(std::declval<FluidState>().moleFraction(phaseIdx, compIdx))
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief The mass fraction of a component in a phase []
     */
    auto massFraction(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(std::declval<FluidState>().massFraction(phaseIdx, compIdx))
    { return fs_->massFraction(phaseIdx, compIdx); }

    /*!
     * \brief The average molar mass of a fluid phase [kg/mol]
     *
     * The average mass is the mean molar mass of a molecule of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[ \bar M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa \f]
     */
    auto averageMolarMass(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().averageMolarMass(phaseIdx))
    { return fs_->averageMolarMass(phaseIdx); }

    /*!
     * \brief The molar concentration of a component in a phase [mol/m^3]
     *
     * This quantity is usually called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    auto molarity(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(std::declval<FluidState>().molarity(phaseIdx, compIdx))
    { return fs_->molarity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    auto fugacity(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(std::declval<FluidState>().fugacity(phaseIdx, compIdx))
    { return fs_->fugacity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    auto fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(std::declval<FluidState>().fugacityCoefficient(phaseIdx, compIdx))
    { return fs_->fugacityCoefficient(phaseIdx, compIdx); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    auto molarVolume(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().molarVolume(phaseIdx))
    { return fs_->molarVolume(phaseIdx); }

    /*!
     * \brief The mass density of a fluid phase [kg/m^3]
     */
    auto density(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().density(phaseIdx))
    { return fs_->density(phaseIdx); }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    auto molarDensity(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().molarDensity(phaseIdx))
    { return fs_->molarDensity(phaseIdx); }

    /*!
     * \brief The temperature of a fluid phase [K]
     */
    const Scalar& temperature(unsigned /*phaseIdx*/) const
    { return temperature_; }

    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    auto pressure(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().pressure(phaseIdx))
    { return fs_->pressure(phaseIdx); }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    auto enthalpy(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().enthalpy(phaseIdx))
    { return fs_->enthalpy(phaseIdx); }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    auto internalEnergy(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().internalEnergy(phaseIdx))
    { return fs_->internalEnergy(phaseIdx); }

    /*!
     * \brief The dynamic viscosity of a fluid phase [Pa s]
     */
    auto viscosity(unsigned phaseIdx) const
        -> decltype(std::declval<FluidState>().viscosity(phaseIdx))
    { return fs_->viscosity(phaseIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature [K] of a fluid phase
     */
    void setTemperature(const Scalar& value)
    { temperature_ = value; }

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
        Valgrind::CheckDefined(temperature_);
    }

protected:
    const FluidState* fs_;
    Scalar temperature_;
};

} // namespace Opm

#endif
