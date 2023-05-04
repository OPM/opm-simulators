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
 * \copydoc Opm::PvsIntensiveQuantities
 */
#ifndef EWOMS_PVS_INTENSIVE_QUANTITIES_HH
#define EWOMS_PVS_INTENSIVE_QUANTITIES_HH

#include "pvsproperties.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>

namespace Opm {
/*!
 * \ingroup PvsModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase primary
 *        variable switching model.
 */
template <class TypeTag>
class PvsIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { switch0Idx = Indices::switch0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { dimWorld = GridView::dimensionworld };

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using MiscibleMultiPhaseComposition = Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem, Evaluation>;
    using ComputeFromReferencePhase = Opm::ComputeFromReferencePhase<Scalar, FluidSystem, Evaluation>;

    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    using EvalPhaseVector = Dune::FieldVector<Evaluation, numPhases>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem>;

    PvsIntensiveQuantities()
    { }

    PvsIntensiveQuantities(const PvsIntensiveQuantities& other) = default;

    PvsIntensiveQuantities& operator=(const PvsIntensiveQuantities& other) = default;

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();

        /////////////
        // set the saturations
        /////////////
        Evaluation sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setSaturation(phaseIdx, priVars.explicitSaturationValue(phaseIdx, timeIdx));
            Opm::Valgrind::CheckDefined(fluidState_.saturation(phaseIdx));
            sumSat += fluidState_.saturation(phaseIdx);
        }
        Opm::Valgrind::CheckDefined(priVars.implicitSaturationIdx());
        Opm::Valgrind::CheckDefined(sumSat);
        fluidState_.setSaturation(priVars.implicitSaturationIdx(), 1.0 - sumSat);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        EvalPhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // set the absolute phase pressures in the fluid state
        const Evaluation& p0 = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, p0 + (pC[phaseIdx] - pC[0]));

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        unsigned lowestPresentPhaseIdx = priVars.lowestPresentPhaseIdx();
        unsigned numNonPresentPhases = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!priVars.phaseIsPresent(phaseIdx))
                ++numNonPresentPhases;
        }

        // now comes the tricky part: calculate phase compositions
        if (numNonPresentPhases == numPhases - 1) {
            // only one phase is present, i.e. the primary variables
            // contain the complete composition of the phase
            Evaluation sumx = 0.0;
            for (unsigned compIdx = 1; compIdx < numComponents; ++compIdx) {
                const Evaluation& x = priVars.makeEvaluation(switch0Idx + compIdx - 1, timeIdx);
                fluidState_.setMoleFraction(lowestPresentPhaseIdx, compIdx, x);
                sumx += x;
            }

            // set the mole fraction of the first component
            fluidState_.setMoleFraction(lowestPresentPhaseIdx, 0, 1 - sumx);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_, paramCache,
                                             lowestPresentPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        else {
            // create the auxiliary constraints
            unsigned numAuxConstraints = numComponents + numNonPresentPhases - numPhases;
            Opm::MMPCAuxConstraint<Evaluation> auxConstraints[numComponents];

            unsigned auxIdx = 0;
            unsigned switchIdx = 0;
            for (; switchIdx < numPhases - 1; ++switchIdx) {
                unsigned compIdx = switchIdx + 1;
                unsigned switchPhaseIdx = switchIdx;
                if (switchIdx >= lowestPresentPhaseIdx)
                    switchPhaseIdx += 1;

                if (!priVars.phaseIsPresent(switchPhaseIdx)) {
                    auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx,
                                               priVars.makeEvaluation(switch0Idx + switchIdx, timeIdx));
                    ++auxIdx;
                }
            }

            for (; auxIdx < numAuxConstraints; ++auxIdx, ++switchIdx) {
                unsigned compIdx = numPhases - numNonPresentPhases + auxIdx;
                auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx,
                                           priVars.makeEvaluation(switch0Idx + switchIdx, timeIdx));
            }

            // both phases are present, i.e. phase compositions are a result of the the
            // gas <-> liquid equilibrium. This is the job of the
            // "MiscibleMultiPhaseComposition" constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_, paramCache,
                                                 priVars.phasePresence(),
                                                 auxConstraints,
                                                 numAuxConstraints,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);
        }

#ifndef NDEBUG
        // make valgrind happy and set the enthalpies to NaN
        if (!enableEnergy) {
            Scalar myNan = std::numeric_limits<Scalar>::quiet_NaN();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState_.setEnthalpy(phaseIdx, myNan);
        }
#endif

        /////////////
        // calculate the remaining quantities
        /////////////

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        // mobilities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobility_[phaseIdx] =
                relativePermeability_[phaseIdx] / fluidState().viscosity(phaseIdx);

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

        fluidState_.checkDefined();
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
    { return mobility_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Evaluation porosity_;
    DimMatrix intrinsicPerm_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Opm

#endif
