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
 * \brief Modules for the ModularFluidState which represent fugacity/chemical potential.
 */
#ifndef OPM_FLUID_STATE_FUGACITY_MODULES_HPP
#define OPM_FLUID_STATE_FUGACITY_MODULES_HPP

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <algorithm>
#include <limits>

namespace Opm {

/*!
 * \brief Module for the modular fluid state which stores the
 *        phase fugacity coefficients explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          unsigned numComponents,
          class Implementation>
class FluidStateExplicitFugacityModule
{
public:
    FluidStateExplicitFugacityModule()
    {
        Valgrind::SetUndefined(fugacityCoefficient_);
    }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    const Scalar& fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(unsigned phaseIdx, unsigned compIdx) const
    { return asImp_().pressure(phaseIdx)*fugacityCoefficient_[phaseIdx][compIdx]*asImp_().moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(unsigned phaseIdx, unsigned compIdx, const Scalar& value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
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
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar fugacityCoefficient_[numPhases][numComponents];
};

/*!
 * \brief Module for the modular fluid state which stores the phase
 *        fugacity coefficients explicitly assuming immiscibility.
 */
template <class Scalar,
          unsigned numPhases,
          unsigned numComponents,
          class Implementation>
class FluidStateImmiscibleFugacityModule
{
    static_assert(numPhases == numComponents,
                  "The number of phases must be the same as the number of (pseudo-) components if you assume immiscibility");

public:
    FluidStateImmiscibleFugacityModule()
    { Valgrind::SetUndefined(fugacityCoefficient_); }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    Scalar fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return (phaseIdx == compIdx)?fugacityCoefficient_[phaseIdx]:std::numeric_limits<Scalar>::infinity(); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(unsigned phaseIdx, unsigned compIdx) const
    { return asImp_().pressure(phaseIdx)*fugacityCoefficient(phaseIdx, compIdx)*asImp_().moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(unsigned phaseIdx, const Scalar& value)
    { fugacityCoefficient_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar fugacityCoefficient_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not store the
 *        fugacities but throws std::logic_error instead.
 */
template <class Scalar>
class FluidStateNullFugacityModule
{
public:
    FluidStateNullFugacityModule()
    { }

    /*!
     * \brief The fugacity coefficient of a component in a phase []
     */
    const Scalar& fugacityCoefficient(unsigned /* phaseIdx */, unsigned /* compIdx */) const
    { throw std::logic_error("Fugacity coefficients are not provided by this fluid state"); }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    const Scalar& fugacity(unsigned /* phaseIdx */, unsigned /* compIdx */) const
    { throw std::logic_error("Fugacities coefficients are not provided by this fluid state"); }

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
