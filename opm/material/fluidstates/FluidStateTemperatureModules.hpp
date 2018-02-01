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
 * \brief Modules for the ModularFluidState which represent temperature.
 */
#ifndef OPM_FLUID_STATE_TEMPERATURE_MODULES_HPP
#define OPM_FLUID_STATE_TEMPERATURE_MODULES_HPP

#include <opm/material/common/Valgrind.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <algorithm>
#include <cassert>

namespace Opm {

/*!
 * \brief Module for the modular fluid state which stores the
 *       temperatures explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateExplicitTemperatureModule
{
public:
    FluidStateExplicitTemperatureModule()
    { Valgrind::SetUndefined(temperature_); }

    /*!
     * \brief The temperature of a fluid phase [-]
     */
    const Scalar& temperature(unsigned phaseIdx) const
    { return temperature_[phaseIdx]; }

    /*!
     * \brief Set the temperature of a phase [-]
     */
    void setTemperature(unsigned phaseIdx, const Scalar& value)
    { temperature_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            temperature_[phaseIdx] = fs.temperature(phaseIdx);
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
        Valgrind::CheckDefined(temperature_);
    }

protected:
    Scalar temperature_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which stores the
 *        temperatures explicitly and assumes thermal equilibrium.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateEquilibriumTemperatureModule
{
public:
    FluidStateEquilibriumTemperatureModule()
    { Valgrind::SetUndefined(temperature_); }

    /*!
     * \brief The temperature of a fluid phase [-]
     */
    const Scalar& temperature(unsigned /*phaseIdx*/) const
    { return temperature_; }

    /*!
     * \brief Set the temperature of a phase [-]
     */
    void setTemperature(const Scalar& value)
    { temperature_ = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        temperature_ = Opm::decay<Scalar>(fs.temperature(/*phaseIdx=*/0));

#ifndef NDEBUG
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(std::abs(Opm::scalarValue(fs.temperature(phaseIdx))
                            - Opm::scalarValue(temperature_)) < 1e-30);
        }
#endif
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
        Valgrind::CheckDefined(temperature_);
    }

protected:
    Scalar temperature_;
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        temperatures but throws std::logic_error instead.
 */
template <class Scalar>
class FluidStateNullTemperatureModule
{
public:
    FluidStateNullTemperatureModule()
    { }

    /*!
     * \brief The temperature of a fluid phase [-]
     */
    const Scalar& temperature(unsigned /* phaseIdx */) const
    { throw std::runtime_error("Temperature is not provided by this fluid state"); }

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
