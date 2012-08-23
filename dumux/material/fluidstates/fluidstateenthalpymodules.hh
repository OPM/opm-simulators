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
 * \brief Modules for the ModularFluidState which represent enthalpy.
 */
#ifndef DUMUX_FLUID_STATE_ENTHALPY_MODULES_HH
#define DUMUX_FLUID_STATE_ENTHALPY_MODULES_HH

#include <dumux/common/valgrind.hh>

#include <dune/common/exceptions.hh>

#include <algorithm>

namespace Dumux
{
/*!
 * \brief Module for the modular fluid state which stores the
 *       enthalpys explicitly.
 */
template <class Scalar, 
          class FluidSystem,
          class Implementation>
class FluidStateExplicitEnthalpyModule
{
    enum { numPhases = FluidSystem::numPhases };

public:
    FluidStateExplicitEnthalpyModule()
    { Valgrind::SetUndefined(enthalpy_); }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return enthalpy_[phaseIdx] - asImp_().pressure(phaseIdx)/asImp_().density(phaseIdx); }


    /*!
     * \brief Set the specific enthalpy of a phase [J/kg]
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    { enthalpy_[phaseIdx] = value; }


    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
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
        Valgrind::CheckDefined(enthalpy_);
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar enthalpy_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not store the
 *        enthalpies but throws Dune::InvalidState instead.
 */
template <class Scalar, 
          class FluidSystem,
          class Implementation>
class FluidStateNullEnthalpyModule
{
public:
    FluidStateNullEnthalpyModule()
    { }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(int phaseIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Internal energy is not provided by this fluid state"); }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    Scalar enthalpy(int phaseIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Enthalpy is not provided by this fluid state"); }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    { }

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

} // end namepace Dumux

#endif
