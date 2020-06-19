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
 * \copydoc Opm::NcpPrimaryVariables
 */
#ifndef EWOMS_NCP_PRIMARY_VARIABLES_HH
#define EWOMS_NCP_PRIMARY_VARIABLES_HH

#include "ncpproperties.hh"

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/common/energymodule.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Math.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup NcpModel
 *
 * \brief Represents the primary variables used by the compositional
 *        multi-phase NCP model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class NcpPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { fugacity0Idx = Indices::fugacity0Idx };

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;

    using NcpFlash = Opm::NcpFlash<Scalar, FluidSystem>;
    using Toolbox = Opm::MathToolbox<Evaluation>;

public:
    NcpPrimaryVariables() : ParentType()
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    NcpPrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    NcpPrimaryVariables(const NcpPrimaryVariables& value) = default;
    NcpPrimaryVariables& operator=(const NcpPrimaryVariables& value) = default;

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium = false)
    {
        using FsToolbox = Opm::MathToolbox<typename FluidState::Scalar>;

#ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            assert(fluidState.temperature(0) == fluidState.temperature(phaseIdx));
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        Opm::CompositionalFluidState<Scalar, FluidSystem> fsFlash;

        // use the externally given fluid state as initial value for
        // the flash calculation
        fsFlash.assign(fluidState);

        // calculate the phase densities
        paramCache.updateAll(fsFlash);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fsFlash, paramCache, phaseIdx);
            fsFlash.setDensity(phaseIdx, rho);
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    FsToolbox::value(fsFlash.saturation(phaseIdx))
                    * FsToolbox::value(fsFlash.molarity(phaseIdx, compIdx));
            }
        }

        // run the flash calculation
        NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState, unsigned refPhaseIdx = 0)
    {
        using FsToolbox = Opm::MathToolbox<typename FluidState::Scalar>;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // assign fugacities.
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fluidState, refPhaseIdx);
        Scalar pRef = FsToolbox::value(fluidState.pressure(refPhaseIdx));
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            // we always compute the fugacities because they are quite exotic quantities
            // and this easily forgotten to be specified
            Scalar fugCoeff =
                FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fluidState,
                                                                              paramCache,
                                                                              refPhaseIdx,
                                                                              compIdx);
            (*this)[fugacity0Idx + compIdx] =
                fugCoeff*fluidState.moleFraction(refPhaseIdx, compIdx)*pRef;
        }

        // assign pressure of first phase
        (*this)[pressure0Idx] = FsToolbox::value(fluidState.pressure(/*phaseIdx=*/0));

        // assign first M - 1 saturations
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            (*this)[saturation0Idx + phaseIdx] = FsToolbox::value(fluidState.saturation(phaseIdx));
    }
};

} // namespace Opm

#endif
