// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Modules for the ModularFluidState which represent density.
 */
#ifndef DUMUX_FLUID_STATE_DENSITY_MODULES_HH
#define DUMUX_FLUID_STATE_DENSITY_MODULES_HH

#include <dumux/common/valgrind.hh>

#include <dune/common/exceptions.hh>

#include <algorithm>

namespace Dumux
{
/*!
 * \brief Module for the modular fluid state which stores the
 *       densitys explicitly.
 */
template <class Scalar, 
          class FluidSystem,
          class Implementation>
class FluidStateExplicitDensityModule
{
    enum { numPhases = FluidSystem::numPhases };

public:
    FluidStateExplicitDensityModule()
    { Valgrind::SetUndefined(density_); }

    /*!
     * \brief The density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(int phaseIdx) const
    { return density_[phaseIdx]/asImp_().averageMolarMass(phaseIdx); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(int phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * \brief Set the density of a phase [kg/m^3]
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = fs.density(phaseIdx);
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
        Valgrind::CheckDefined(density_);
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar density_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        densitys but throws Dune::InvalidState instead.
 */
template <class Scalar, 
          class FluidSystem,
          class Implementation>
class FluidStateNullDensityModule
{
public:
    FluidStateNullDensityModule()
    { }

    /*!
     * \brief The density of a fluid phase [kg/m^3]
     */
    Scalar density(int phaseIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Density is not provided by this fluid state"); }


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
