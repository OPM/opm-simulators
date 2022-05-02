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
 * \copydoc Opm::FlashPrimaryVariables
 */
#ifndef EWOMS_FLASH_PRIMARY_VARIABLES_HH
#define EWOMS_FLASH_PRIMARY_VARIABLES_HH

#include "flashindices.hh"
#include "flashproperties.hh"

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/common/energymodule.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <iostream>

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief Represents the primary variables used by the compositional
 *        flow model based on flash calculations.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class FlashPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    // primary variable indices
    enum { z0Idx = Indices::z0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    using Toolbox = typename Opm::MathToolbox<Evaluation>;
    using EnergyModule = Opm::EnergyModule<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>;

public:
    FlashPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    FlashPrimaryVariables(Scalar value) : ParentType(value)
    {
        Opm::Valgrind::CheckDefined(value);
        Opm::Valgrind::SetDefined(*this);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    FlashPrimaryVariables(const FlashPrimaryVariables& value) = default;
    FlashPrimaryVariables& operator=(const FlashPrimaryVariables& value) = default;

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams&,
                                bool = false)
    {
        // there is no difference between naive and mass conservative
        // assignment in the flash model. (we only need the total
        // concentrations.)
        assignNaive(fluidState);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // reset everything
        (*this) = 0.0;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // determine the phase presence.
        Dune::FieldVector<Scalar, numComponents> z(0.0);
        Scalar sumMoles = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar tmp = Opm::getValue(fluidState.molarity(phaseIdx, compIdx) * fluidState.saturation(phaseIdx));
                z[compIdx] += Opm::max(tmp, 1e-8);
                sumMoles += tmp;
            }
        }
        z /= sumMoles;

        for (int i = 0; i < numComponents - 1; ++i)
            (*this)[z0Idx + i] = z[i];
        (*this)[pressure0Idx] = Opm::getValue(fluidState.pressure(0));
    }

    /*!
     * \brief Prints the names of the primary variables and their values.
     *
     * \param os The \c std::ostream which should be used for the output.
     */
    void print(std::ostream& os = std::cout) const
    {
        os << "(p_" << FluidSystem::phaseName(0) << " = "
           << this->operator[](pressure0Idx);
        for (unsigned compIdx = 0; compIdx < numComponents - 2; ++compIdx) {
            os << ", z_" << FluidSystem::componentName(compIdx) << " = "
               << this->operator[](z0Idx + compIdx);
        }
        os << ")" << std::flush;
    }
};

} // namespace Opm

#endif
