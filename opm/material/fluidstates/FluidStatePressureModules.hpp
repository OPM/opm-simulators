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
 * \brief Modules for the ModularFluidState which represent pressure.
 */
#ifndef OPM_FLUID_STATE_PRESSURE_MODULES_HPP
#define OPM_FLUID_STATE_PRESSURE_MODULES_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <algorithm>

namespace Opm {

/*!
 * \brief Module for the modular fluid state which stores the
 *       pressures explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateExplicitPressureModule
{
public:
    FluidStateExplicitPressureModule()
    { Valgrind::SetUndefined(pressure_); }

    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    const Scalar& pressure(unsigned phaseIdx) const
    { return pressure_[phaseIdx]; }


    /*!
     * \brief Set the pressure of a phase [Pa]
     */
    void setPressure(unsigned phaseIdx, const Scalar& value)
    { pressure_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            pressure_[phaseIdx] = Opm::decay<Scalar>(fs.pressure(phaseIdx));
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
        Valgrind::CheckDefined(pressure_);
    }

protected:
    Scalar pressure_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        pressures but throws std::logic_error instead.
 */
template <class Scalar>
class FluidStateNullPressureModule
{
public:
    FluidStateNullPressureModule()
    { }

    /*!
     * \brief The pressure of a fluid phase [Pa]
     */
    const Scalar& pressure(unsigned /* phaseIdx */) const
    { throw std::logic_error("Pressure is not provided by this fluid state"); }


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
