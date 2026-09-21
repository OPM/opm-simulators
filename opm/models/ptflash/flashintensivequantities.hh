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
#ifndef OPM_FLASH_INTENSIVE_QUANTITIES_HH
#define OPM_FLASH_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/Constants.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/PTFlashMethod.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/models/flash/flashproperties.hh>

#include <opm/models/ptflash/flashindices.hh>
#include <opm/models/ptflash/flashparameters.hh>

#include <fmt/format.h>

#include <array>
#include <iterator>
#include <string>

namespace Opm {

/*!
 * \ingroup FlashModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the intensive quantities of the ptflash-based compositional multi-phase model
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
    enum { z0Idx = Indices::z0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    enum { dimWorld = GridView::dimensionworld };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { water0Idx = Indices::water0Idx};

    static constexpr bool waterEnabled = Indices::waterEnabled;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    /// Minimum hydrocarbon share of the pore space, including cells that would
    /// otherwise hold water alone. Component-storage derivatives scale with this
    /// share, so it cannot vanish. The value matches the composition floor, and
    /// the displaced water remains below the precision of single-precision
    /// restart output.
    static constexpr Scalar hydrocarbonFloor = 1.0e-8;

    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;

    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using DiffusionIntensiveQuantities = ::Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = ::Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;

    FlashIntensiveQuantities() = default;

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

        const Scalar flashTolerance = Parameters::Get<Parameters::FlashTolerance<Scalar>>();
        const int flashVerbosity = Parameters::Get<Parameters::FlashVerbosity>();
        const auto ptFlashMethod =
            ptFlashMethodFromString(Parameters::Get<Parameters::FlashTwoPhaseMethod>());
        // TODO: the formulation here is still to begin with XMF and YMF values to derive ZMF value
        // TODO: we should check how we update ZMF in the newton update, since it is the primary variables.

        // extract the total molar densities of the components
        ComponentVector z(0.);
        {
            Evaluation lastZ = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                z[compIdx] = priVars.makeEvaluation(z0Idx + compIdx, timeIdx);
                lastZ -= z[compIdx];
            }
            z[numComponents - 1] = lastZ;

            Evaluation sumz = 0.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                // Clamp only the value; preserve derivatives. Replacing the Evaluation with
                // max() when the bound applies removes composition derivatives from a
                // vanished component's conservation equation and makes the cell Jacobian
                // block singular.
                if (z[compIdx] < 1e-8) {
                    z[compIdx].setValue(1e-8);
                }
                sumz += z[compIdx];
            }
            z /= sumz;
        }

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluidState_.setMoleFraction(compIdx, z[compIdx]);
        }

        Evaluation p = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setPressure(phaseIdx, p);
        }

        // Get initial K and L from storage initially (if enabled)
        const auto* hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        if (hint) {
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation& Ktmp = hint->fluidState().K(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = hint->fluidState().L();
             fluidState_.setLvalue(Ltmp);
        }
        else if (timeIdx == 0 && elemCtx.thermodynamicHint(dofIdx, 1)) {
             // checking the storage cache
             const auto& hint2 = elemCtx.thermodynamicHint(dofIdx, 1);
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation& Ktmp = hint2->fluidState().K(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = hint2->fluidState().L();
             fluidState_.setLvalue(Ltmp);
        }
        else {
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation Ktmp = fluidState_.wilsonK_(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = -1.0;
             fluidState_.setLvalue(Ltmp);
         }

        /////////////
        // Compute the phase compositions and densities
        /////////////
        if (flashVerbosity >= 1) {
            OpmLog::debug(fmt::format("Updating the intensive quantities for cell {}",
                                      elemCtx.globalSpaceIndex(dofIdx, timeIdx)));
        }
        const auto& eos_type = problem.getEosType();
        FlashSolver::solve(fluidState_, ptFlashMethod, flashTolerance, eos_type, flashVerbosity);

        if (flashVerbosity >= 5) {
            std::string phaseCompositions;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fmt::format_to(
                    std::back_inserter(phaseCompositions),
                    "  component {}: x = {}, y = {}\n",
                    compIdx,
                    getValue(fluidState_.moleFraction(FluidSystem::oilPhaseIdx, compIdx)),
                    getValue(fluidState_.moleFraction(FluidSystem::gasPhaseIdx, compIdx)));
            }
            OpmLog::debug(fmt::format("After the flash for cell {}: liquid fraction = {}\n{}",
                                      elemCtx.globalSpaceIndex(dofIdx, timeIdx),
                                      getValue(fluidState_.L()),
                                      phaseCompositions));
        }

        // Update phases
        typename FluidSystem::template ParameterCache<Evaluation> paramCache(eos_type);
        paramCache.updatePhase(fluidState_, FluidSystem::oilPhaseIdx);
        paramCache.updatePhase(fluidState_, FluidSystem::gasPhaseIdx);

        // Update saturation
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            Sw = priVars.makeEvaluation(water0Idx, timeIdx);
        }
        const Evaluation L = fluidState_.L();
        const Evaluation Vm_L = paramCache.correctedMolarVolume(FluidSystem::oilPhaseIdx);
        const Evaluation Vm_V = paramCache.correctedMolarVolume(FluidSystem::gasPhaseIdx);

        // Every component storage term is proportional to the pore-space share
        // occupied by hydrocarbons. If a cell holds only water, a zero share
        // removes the composition from the residual: the component equations
        // then depend on no composition variable, and the cell's diagonal
        // Jacobian block is singular. Every equilibrated cell below the
        // water-oil contact starts in that state. Floor the value while retaining
        // its derivatives, as for the composition above.
        Evaluation hydrocarbon = 1 - Sw;
        if (hydrocarbon < hydrocarbonFloor) {
            hydrocarbon.setValue(hydrocarbonFloor);
        }
        Evaluation So = max(hydrocarbon * (L * Vm_L / ( L * Vm_L + (1 - L) * Vm_V)), 0.0);
        Evaluation Sg = max(hydrocarbon - So, 0.0);
        const Scalar sumS = getValue(So) + getValue(Sg) + getValue(Sw);
        So /= sumS;
        Sg /= sumS;

        fluidState_.setSaturation(FluidSystem::oilPhaseIdx, So);
        fluidState_.setSaturation(FluidSystem::gasPhaseIdx, Sg);
        if constexpr (waterEnabled) {
            Sw /= sumS;
            fluidState_.setSaturation(FluidSystem::waterPhaseIdx, Sw);
        }

        // The compressibility factor comes from the unshifted EOS root, unlike
        // the saturations above.
        const Scalar R = Opm::Constants<Scalar>::R;
        const Evaluation Z_L = (paramCache.molarVolume(FluidSystem::oilPhaseIdx) *
                                fluidState_.pressure(FluidSystem::oilPhaseIdx)) /
                               (R * fluidState_.temperature(FluidSystem::oilPhaseIdx));
        const Evaluation Z_V = (paramCache.molarVolume(FluidSystem::gasPhaseIdx) *
                                fluidState_.pressure(FluidSystem::gasPhaseIdx)) /
                               (R * fluidState_.temperature(FluidSystem::gasPhaseIdx));
        fluidState_.setCompressFactor(FluidSystem::oilPhaseIdx, Z_L);
        fluidState_.setCompressFactor(FluidSystem::gasPhaseIdx, Z_V);

        if (flashVerbosity >= 5) {
            OpmLog::debug(fmt::format("Flash phase properties for cell {}: "
                                      "oil saturation = {}, gas saturation = {}, "
                                      "oil molar volume = {}, gas molar volume = {}",
                                      elemCtx.globalSpaceIndex(dofIdx, timeIdx),
                                      getValue(So),
                                      getValue(Sg),
                                      getValue(Vm_L),
                                      getValue(Vm_V)));
        }

        /////////////
        // Compute rel. perm and viscosity and densities
        /////////////
        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);

        // calculate relative permeability
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // set the phase viscosity and density
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx == static_cast<unsigned int>(FluidSystem::oilPhaseIdx) ||
                phaseIdx == static_cast<unsigned int>(FluidSystem::gasPhaseIdx))
            {
                paramCache.updatePhase(fluidState_, phaseIdx);
            }

            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);

            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);

            const Evaluation& rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);
        }

        /////////////
        // Compute the remaining quantities
        /////////////

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

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
    { return mobility_[phaseIdx]; }

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
