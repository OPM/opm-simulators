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
 * \brief Modules for the ModularFluidState which represent saturation.
 */
#ifndef EWOMS_FLUID_STATE_SATURATION_MODULES_HH
#define EWOMS_FLUID_STATE_SATURATION_MODULES_HH

#include <ewoms/common/valgrind.hh>

#include <dune/common/exceptions.hh>

#include <algorithm>

namespace Ewoms {

/*!
 * \brief Module for the modular fluid state which stores the
 *       saturations explicitly.
 */
template <class Scalar,
          class FluidSystem,
          class Implementation>
class FluidStateExplicitSaturationModule
{
    enum { numPhases = FluidSystem::numPhases };

public:
    FluidStateExplicitSaturationModule()
    { Valgrind::SetUndefined(saturation_); }

    /*!
     * \brief The saturation of a fluid phase [-]
     */
    Scalar saturation(int phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Set the saturation of a phase [-]
     */
    void setSaturation(int phaseIdx, Scalar value)
    { saturation_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            saturation_[phaseIdx] = fs.saturation(phaseIdx);
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
        Valgrind::CheckDefined(saturation_);
    }

protected:
    Scalar saturation_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        saturations but throws Dune::InvalidState instead.
 */
template <class Scalar,
          class FluidSystem,
          class Implementation>
class FluidStateNullSaturationModule
{
public:
    FluidStateNullSaturationModule()
    { }

    /*!
     * \brief The saturation of a fluid phase [-]
     */
    Scalar saturation(int phaseIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Saturation is not provided by this fluid state"); }

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

} // end namespace Ewoms

#endif
