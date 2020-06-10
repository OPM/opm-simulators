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
 * \copydoc Opm::NcpIntensiveQuantities
 */
#ifndef EWOMS_NCP_INTENSIVE_QUANTITIES_HH
#define EWOMS_NCP_INTENSIVE_QUANTITIES_HH

#include "ncpproperties.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/CompositionFromFugacities.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
/*!
 * \ingroup NcpModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase NCP model.
 */
template <class TypeTag>
class NcpIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { dimWorld = GridView::dimensionworld };

    using CompositionFromFugacitiesSolver = Opm::CompositionFromFugacities<Scalar, FluidSystem, Evaluation>;
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, /*storeEnthalpy=*/enableEnergy>;
    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DiffusionIntensiveQuantities = Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;

public:
    NcpIntensiveQuantities()
    {}

    NcpIntensiveQuantities(const NcpIntensiveQuantities& other) = default;

    NcpIntensiveQuantities& operator=(const NcpIntensiveQuantities& other)  = default;

    /*!
     * \brief IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx,
                unsigned dofIdx,
                unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        ParentType::checkDefined();

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // set the phase saturations
        Evaluation sumSat = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            const Evaluation& val = priVars.makeEvaluation(saturation0Idx + phaseIdx, timeIdx);
            fluidState_.setSaturation(phaseIdx, val);
            sumSat += val;
        }
        fluidState_.setSaturation(numPhases - 1, 1.0 - sumSat);
        Opm::Valgrind::CheckDefined(sumSat);

        // set the fluid phase temperature
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        // retrieve capillary pressure parameters
        const auto& problem = elemCtx.problem();
        const MaterialLawParams& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        // calculate capillary pressures
        Evaluation capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState_);
        // add to the pressure of the first fluid phase
        const Evaluation& pressure0 = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, pressure0 + (capPress[phaseIdx] - capPress[0]));

        ComponentVector fug;
        // retrieve component fugacities
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            fug[compIdx] = priVars.makeEvaluation(fugacity0Idx + compIdx, timeIdx);

        // calculate phase compositions
        const auto *hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // initial guess
            if (hint) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    // use the hint for the initial mole fraction!
                    const Evaluation& moleFracIJ = hint->fluidState().moleFraction(phaseIdx, compIdx);
                    fluidState_.setMoleFraction(phaseIdx, compIdx, moleFracIJ);
                }
            }
            else // !hint
                CompositionFromFugacitiesSolver::guessInitial(fluidState_, phaseIdx, fug);

            // calculate the phase composition from the component
            // fugacities
            CompositionFromFugacitiesSolver::solve(fluidState_, paramCache, phaseIdx, fug);
        }

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(porosity_);

        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);

        // dynamic viscosities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx]/mu;
        }

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        checkDefined();
    }

    /*!
     * \brief ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * \brief IntensiveQuantities::checkDefined
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        ParentType::checkDefined();

        Opm::Valgrind::CheckDefined(porosity_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        fluidState_.checkDefined();
#endif
    }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Opm

#endif
