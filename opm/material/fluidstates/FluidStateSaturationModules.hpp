// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \brief Modules for the ModularFluidState which represent saturation.
 */
#ifndef OPM_FLUID_STATE_SATURATION_MODULES_HPP
#define OPM_FLUID_STATE_SATURATION_MODULES_HPP

#include <opm/material/common/ErrorMacros.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <algorithm>

namespace Opm {

/*!
 * \brief Module for the modular fluid state which stores the
 *       saturations explicitly.
 */
template <class Scalar,
          int numPhases,
          class Implementation>
class FluidStateExplicitSaturationModule
{
public:
    FluidStateExplicitSaturationModule()
    { Valgrind::SetUndefined(saturation_); }

    /*!
     * \brief The saturation of a fluid phase [-]
     */
    const Scalar& saturation(unsigned phaseIdx) const
    { return saturation_[phaseIdx]; }

    /*!
     * \brief Set the saturation of a phase [-]
     */
    void setSaturation(unsigned phaseIdx, const Scalar& value)
    { saturation_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        typedef typename FluidState::Scalar FsScalar;
        typedef Opm::MathToolbox<FsScalar> FsToolbox;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            saturation_[phaseIdx] = FsToolbox::template toLhs<Scalar>(fs.saturation(phaseIdx));
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
 *        saturations but throws std::logic_error instead.
 */
template <class Scalar>
class FluidStateNullSaturationModule
{
public:
    FluidStateNullSaturationModule()
    { }

    /*!
     * \brief The saturation of a fluid phase [-]
     */
    const Scalar& saturation(int /* phaseIdx */) const
    { OPM_THROW(std::runtime_error, "Saturation is not provided by this fluid state"); }

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
