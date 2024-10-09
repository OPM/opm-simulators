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
 * \copydoc Opm::FlashIntensiveQuantities
 */
#ifndef EWOMS_FLASH_INTENSIVE_QUANTITIES_HH
#define EWOMS_FLASH_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/models/common/diffusionmodule.hh>
#include <opm/models/common/energymodule.hh>

#include <opm/models/flash/flashindices.hh>
#include <opm/models/flash/flashparameters.hh>
#include <opm/models/flash/flashproperties.hh>

namespace Opm {

/*!
 * \ingroup FlashModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the intensive quantities of the flash-based compositional multi-phase model
 */
template <class TypeTag>
class FlashIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;

    // primary variable indices
    enum { cTot0Idx = Indices::cTot0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { dimWorld = GridView::dimensionworld };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;

    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;

    FlashIntensiveQuantities()
    { }

    FlashIntensiveQuantities(const FlashIntensiveQuantities& other) = default;

    FlashIntensiveQuantities& operator=(const FlashIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();
        Scalar flashTolerance = Parameters::Get<Parameters::FlashTolerance<Scalar>>();

        // extract the total molar densities of the components
        ComponentVector cTotal;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            cTotal[compIdx] = priVars.makeEvaluation(cTot0Idx + compIdx, timeIdx);

        const auto *hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        if (hint) {
            // use the same fluid state as the one of the hint, but
            // make sure that we don't overwrite the temperature
            // specified by the primary variables
            Evaluation T = fluidState_.temperature(/*phaseIdx=*/0);
            fluidState_.assign(hint->fluidState());
            fluidState_.setTemperature(T);
        }
        else
            FlashSolver::guessInitial(fluidState_, cTotal);

        // compute the phase compositions, densities and pressures
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        const MaterialLawParams& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        FlashSolver::template solve<MaterialLaw>(fluidState_,
                                                 materialParams,
                                                 paramCache,
                                                 cTotal,
                                                 flashTolerance);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        // set the phase viscosities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState_, phaseIdx);

            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Opm::Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        /////////////
        // calculate the remaining quantities
        /////////////

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    std::array<Evaluation,numPhases> relativePermeability_;
    std::array<Evaluation,numPhases> mobility_;
};

} // namespace Opm

#endif
