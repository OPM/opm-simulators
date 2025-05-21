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
#ifndef OPM_PTFLASH_PRIMARY_VARIABLES_HH
#define OPM_PTFLASH_PRIMARY_VARIABLES_HH

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/energymodule.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/ptflash/flashindices.hh>

#include <ostream>

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
    enum { water0Idx = Indices::water0Idx };

    static constexpr bool waterEnabled = Indices::waterEnabled;

    // phase indices
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    using Toolbox = MathToolbox<Evaluation>;
    using EnergyModule = ::Opm::EnergyModule<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>;

public:
    FlashPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    FlashPrimaryVariables(const FlashPrimaryVariables& value) = default;
    FlashPrimaryVariables& operator=(const FlashPrimaryVariables& value) = default;

    using ParentType::operator=; //!< Import base class assignment operators.

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

        // assign components total fraction
        for (int i = 0; i < numComponents - 1; ++i)
            (*this)[z0Idx + i] = getValue(fluidState.moleFraction(i));

        // assign pressure
        (*this)[pressure0Idx] = getValue(fluidState.pressure(oilPhaseIdx));

        // assign water saturation
        if constexpr (waterEnabled) {
            (*this)[water0Idx] = getValue(fluidState.saturation(waterPhaseIdx));
        }
    }

    /*!
     * \brief Prints the names of the primary variables and their values.
     *
     * \param os The \c std::ostream which should be used for the output.
     */
    void print(std::ostream& os) const
    {
        os << "(p_" << FluidSystem::phaseName(FluidSystem::oilPhaseIdx) << " = "
           << (*this)[pressure0Idx];
        for (unsigned compIdx = 0; compIdx < numComponents - 2; ++compIdx) {
            os << ", z_" << FluidSystem::componentName(compIdx) << " = "
               << (*this)[z0Idx + compIdx];
        }
        if constexpr (waterEnabled) {
            os << ", S_w = " << (*this)[water0Idx];
        }
        os << ")" << std::flush;
    }
};

} // namespace Opm

#endif
