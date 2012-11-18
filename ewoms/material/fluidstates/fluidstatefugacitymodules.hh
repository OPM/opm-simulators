// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief Modules for the ModularFluidState which represent fugacity/chemical potential.
 */
#ifndef EWOMS_FLUID_STATE_FUGACITY_MODULES_HH
#define EWOMS_FLUID_STATE_FUGACITY_MODULES_HH

#include <ewoms/common/valgrind.hh>

#include <dune/common/exceptions.hh>

#include <algorithm>
#include <limits>

namespace Ewoms
{
/*!
 * \brief Module for the modular fluid state which stores the
 *        phase fugacity coefficients explicitly.
 */
template <class Scalar,
          class FluidSystem,
          class Implementation>
class FluidStateExplicitFugacityModule
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

public:
    FluidStateExplicitFugacityModule()
    {
        Valgrind::SetUndefined(fugacityCoefficient_);
    }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return asImp_().pressure(phaseIdx)*fugacityCoefficient_[phaseIdx][compIdx]*asImp_().moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(int phaseIdx, int compIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
            }
        }
    }

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
        Valgrind::CheckDefined(fugacityCoefficient_);
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar fugacityCoefficient_[numPhases][numComponents];
};

/*!
 * \brief Module for the modular fluid state which stores the phase
 *        fugacity coefficients explicitly assuming immiscibility.
 */
template <class Scalar,
          class FluidSystem,
          class Implementation>
class FluidStateImmiscibleFugacityModule
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    static_assert((int) numPhases == (int) numComponents,
                  "The number of phases must be the same as the number of (pseudo-) components if you assume immiscibility");

public:
    FluidStateImmiscibleFugacityModule()
    { Valgrind::SetUndefined(fugacityCoefficient_); }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return (phaseIdx == compIdx)?fugacityCoefficient_[phaseIdx]:std::numeric_limits<Scalar>::infinity(); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return asImp_().pressure(phaseIdx)*fugacityCoefficient(phaseIdx, compIdx)*asImp_().moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(int phaseIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fugacityCoefficient_[phaseIdx] = fs.fugacityCoefficient(phaseIdx, /*compIdx=*/phaseIdx);
        }
    }

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
        Valgrind::CheckDefined(fugacityCoefficient_);
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar fugacityCoefficient_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not store the
 *        fugacitys but throws Dune::InvalidState instead.
 */
template <class Scalar,
          class FluidSystem,
          class Implementation>
class FluidStateNullFugacityModule
{
public:
    FluidStateNullFugacityModule()
    { }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Fugacity coefficients are not provided by this fluid state"); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Fugacities coefficients are not provided by this fluid state"); }

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    { }
};


} // end namespace Ewoms

#endif
