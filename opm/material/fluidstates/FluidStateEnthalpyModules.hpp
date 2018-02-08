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
 *
 * \brief Modules for the ModularFluidState which represent enthalpy.
 */
#ifndef OPM_FLUID_STATE_ENTHALPY_MODULES_HPP
#define OPM_FLUID_STATE_ENTHALPY_MODULES_HPP

#include <opm/material/common/Exceptions.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <algorithm>

namespace Opm {
/*!
 * \brief Module for the modular fluid state which stores the
 *       enthalpies explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateExplicitEnthalpyModule
{
public:
    FluidStateExplicitEnthalpyModule()
    { Valgrind::SetUndefined(enthalpy_); }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    const Scalar& enthalpy(unsigned phaseIdx) const
    { return enthalpy_[phaseIdx]; }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    Scalar internalEnergy(unsigned phaseIdx) const
    { return enthalpy_[phaseIdx] - asImp_().pressure(phaseIdx)/asImp_().density(phaseIdx); }

    /*!
     * \brief Set the specific enthalpy of a phase [J/kg]
     */
    void setEnthalpy(unsigned phaseIdx, const Scalar& value)
    { enthalpy_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            enthalpy_[phaseIdx] = Opm::decay<Scalar>(fs.enthalpy(phaseIdx));
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
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar enthalpy_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not store the
 *        enthalpies but returns 0 instead.
 *
 * Also, the returned values are marked as undefined in Valgrind... */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateNullEnthalpyModule
{
public:
    FluidStateNullEnthalpyModule()
    { }

    /*!
     * \brief The specific internal energy of a fluid phase [J/kg]
     */
    const Scalar& internalEnergy(unsigned /* phaseIdx */) const
    {
        static Scalar tmp = 0;
        Valgrind::SetUndefined(tmp);
        return tmp;
    }

    /*!
     * \brief The specific enthalpy of a fluid phase [J/kg]
     */
    const Scalar& enthalpy(unsigned /* phaseIdx */) const
    {
        static Scalar tmp = 0;
        Valgrind::SetUndefined(tmp);
        return tmp;
    }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& /* fs */)
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

} // namespace Opm

#endif
