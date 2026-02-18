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
 * \copydoc Opm::BlackOilIntensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH
#define EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/blackoil/blackoilbioeffectsmodules.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoildiffusionmodule.hh>
#include <opm/models/blackoil/blackoildispersionmodule.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoilextbomodules.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/common/directionalmobility.hh>

#include <opm/utility/CopyablePtr.hpp>

#include <array>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
    , public BlackOilDiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public BlackOilDispersionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDispersion>() >
    , public BlackOilSolventIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableSolvent>()>
    , public BlackOilExtboIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableExtbo>()>
    , public BlackOilPolymerIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnablePolymer>()>
    , public BlackOilFoamIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableFoam>()>
    , public BlackOilBrineIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableBrine>()>
    , public BlackOilEnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnergyModuleType>()>
    , public BlackOilBioeffectsIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableBioeffects>()>
    , public BlackOilConvectiveMixingIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableConvectiveMixing>()>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };
    enum { enableConvectiveMixing = getPropValue<TypeTag, Properties::EnableConvectiveMixing>() };
    enum { enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>() };
    enum { enableMICP = Indices::enableMICP };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;
    static constexpr bool waterEnabled = Indices::waterEnabled;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;

    using Toolbox = MathToolbox<Evaluation>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = BlackOilDiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using DispersionIntensiveQuantities = BlackOilDispersionIntensiveQuantities<TypeTag, enableDispersion>;

    using DirectionalMobilityPtr = Utility::CopyablePtr<DirectionalMobility<TypeTag>>;
    using BrineModule = BlackOilBrineModule<TypeTag>;
    using BrineIntQua = BlackOilBrineIntensiveQuantities<TypeTag, enableSaltPrecipitation>;
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag>;
    using BioeffectsIntQua = BlackOilBioeffectsIntensiveQuantities<TypeTag, enableBioeffects>;

public:
    using FluidState = BlackOilFluidState<Evaluation,
                                          FluidSystem,
                                          energyModuleType == EnergyModules::ConstantTemperature,
                                          (energyModuleType == EnergyModules::FullyImplicitThermal || energyModuleType == EnergyModules::SequentialImplicitThermal),
                                          compositionSwitchEnabled,
                                          enableVapwat,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          enableDisgasInWater,
                                          Indices::numPhases>;
    using ScalarFluidState = BlackOilFluidState<Scalar,
                                                FluidSystem,
                                                energyModuleType == EnergyModules::ConstantTemperature,
                                                (energyModuleType == EnergyModules::FullyImplicitThermal || 
                                                    energyModuleType == EnergyModules::SequentialImplicitThermal),
                                                compositionSwitchEnabled,
                                                enableVapwat,
                                                enableBrine,
                                                enableSaltPrecipitation,
                                                enableDisgasInWater,
                                                Indices::numPhases>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    BlackOilIntensiveQuantities()
    {
        if constexpr (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
        }
        if constexpr (enableVapwat) {
            fluidState_.setRvw(0.0);
        }
        if constexpr (enableDisgasInWater) {
            fluidState_.setRsw(0.0);
        }
    }
    BlackOilIntensiveQuantities(const BlackOilIntensiveQuantities& other) = default;

    BlackOilIntensiveQuantities& operator=(const BlackOilIntensiveQuantities& other) = default;

    void updateTempSalt(const Problem& problem,
                        const PrimaryVariables& priVars,
                        const unsigned globalSpaceIdx,
                        const unsigned timeIdx,
                        const LinearizationType& lintype)
    {
        asImp_().updateTemperature_(problem, priVars, globalSpaceIdx, timeIdx, lintype);
        if constexpr (enableBrine) {
            asImp_().updateSaltConcentration_(priVars, timeIdx, lintype);
        }
    }

    void updateSaturations(const PrimaryVariables& priVars,
                           const unsigned timeIdx,
                           [[maybe_unused]] const LinearizationType lintype)
    {
        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                assert(Indices::waterSwitchIdx >= 0);
                if constexpr (Indices::waterSwitchIdx >= 0) {
                    Sw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
                }
            }
            else if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rsw ||
                     priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Disabled)
            {
                // water is enabled but is not a primary variable i.e. one component/phase case
                // or two-phase water + gas with only water present
                Sw = 1.0;
            } // else i.e. for MeaningWater() = Rvw, Sw is still 0.0;
        }
        Evaluation Sg = 0.0;
        if constexpr (gasEnabled) {
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                assert(Indices::compositionSwitchIdx >= 0);
                if constexpr (compositionSwitchEnabled) {
                    Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
                }
            }
            else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                Sg = 1.0 - Sw;
            }
            else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Disabled) {
                if constexpr (waterEnabled) {
                    Sg = 1.0 - Sw; // two phase water + gas
                } else {
                    // one phase case
                    Sg = 1.0;
                }
            }
        }
        Valgrind::CheckDefined(Sg);
        Valgrind::CheckDefined(Sw);

        Evaluation So = 1.0 - Sw - Sg;

        // deal with solvent
        if constexpr (enableSolvent) {
            if (priVars.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
                }
                else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    Sg -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
                }
            }
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            fluidState_.setSaturation(waterPhaseIdx, Sw);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            fluidState_.setSaturation(gasPhaseIdx, Sg);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            fluidState_.setSaturation(oilPhaseIdx, So);
        }
    }

    template <class ...Args>
    void updateRelpermAndPressures(const Problem& problem,
                                   const PrimaryVariables& priVars,
                                   const unsigned globalSpaceIdx,
                                   const unsigned timeIdx,
                                   const LinearizationType& lintype)
    {

        // Solvent saturation manipulation:
        // After this, gas saturation will actually be (gas sat + solvent sat)
        // until set back to just gas saturation in the corresponding call to
        // solventPostSatFuncUpdate_() further down.
        if constexpr (enableSolvent) {
            asImp_().solventPreSatFuncUpdate_(priVars, timeIdx, lintype);
        }

        // Phase relperms.
        problem.template updateRelperms<FluidState, Args...>(mobility_, dirMob_, fluidState_, globalSpaceIdx);

        // now we compute all phase pressures
        using EvalArr = std::array<Evaluation, numPhases>;
        EvalArr pC;
        const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        MaterialLaw::template capillaryPressures<EvalArr, FluidState, Args...>(pC, materialParams, fluidState_);

        // scaling the capillary pressure due to porosity changes
        if constexpr (enableBrine) {
            if (BrineModule::hasPcfactTables() &&
                priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp)
            {
                const unsigned satnumRegionIdx = problem.satnumRegionIndex(globalSpaceIdx);
                const Evaluation Sp = priVars.makeEvaluation(Indices::saltConcentrationIdx, timeIdx);
                const Evaluation porosityFactor  = min(1.0 - Sp, 1.0); //phi/phi_0
                const auto& pcfactTable = BrineModule::pcfactTable(satnumRegionIdx);
                const Evaluation pcFactor = pcfactTable.eval(porosityFactor, /*extrapolation=*/true);
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        pC[phaseIdx] *= pcFactor;
                    }
                }
            }
        }
        else if constexpr (enableBioeffects) {
            if (BioeffectsModule::hasPcfactTables() && referencePorosity_ > 0) {
                unsigned satnumRegionIdx = problem.satnumRegionIndex(globalSpaceIdx);
                const Evaluation Sb = priVars.makeEvaluation(Indices::biofilmVolumeFractionIdx, timeIdx);
                const Evaluation porosityFactor  = min(1.0 - Sb/referencePorosity_, 1.0); //phi/phi_0
                const auto& pcfactTable = BioeffectsModule::pcfactTable(satnumRegionIdx);
                const Evaluation pcFactor = pcfactTable.eval(porosityFactor, /*extrapolation=*/true);
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        pC[phaseIdx] *= pcFactor;
                    }
                }
            }
        }

        // oil is the reference phase for pressure
        if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
                }
            }
        }
        else if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pw) {
            const Evaluation& pw = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fluidState_.setPressure(phaseIdx, pw + (pC[phaseIdx] - pC[waterPhaseIdx]));
                }
            }
        }
        else {
            assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
                }
            }
        }

        // Update the Saturation functions for the blackoil solvent module.
        // Including setting gas saturation back to hydrocarbon gas saturation.
        // Note that this depend on the pressures, so it must be called AFTER the pressures
        // have been updated.
        if constexpr (enableSolvent) {
            asImp_().solventPostSatFuncUpdate_(problem, priVars, globalSpaceIdx, timeIdx, lintype);
        }
    }

    void updateRsRvRsw(const Problem& problem, const PrimaryVariables& priVars, const unsigned globalSpaceIdx, const unsigned timeIdx)
    {
        const unsigned pvtRegionIdx = priVars.pvtRegionIndex();

        const Scalar RvMax = FluidSystem::enableVaporizedOil()
            ? problem.maxOilVaporizationFactor(timeIdx, globalSpaceIdx)
            : 0.0;
        const Scalar RsMax = FluidSystem::enableDissolvedGas()
            ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
            : 0.0;
        const Scalar RswMax = FluidSystem::enableDissolvedGasInWater()
            ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
            : 0.0;

        Evaluation SoMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            SoMax = max(fluidState_.saturation(oilPhaseIdx),
                        problem.maxOilSaturation(globalSpaceIdx));
        }

        // take the meaning of the switching primary variable into account for the gas
        // and oil phase compositions

        if constexpr (compositionSwitchEnabled) {
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
                const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
                fluidState_.setRs(Rs);
            }
            else {
                if (FluidSystem::enableDissolvedGas()) { // Add So > 0? i.e. if only water set rs = 0)
                    const Evaluation& RsSat = enableExtbo ? asImp_().rs() :
                        FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                                oilPhaseIdx,
                                                                pvtRegionIdx,
                                                                SoMax);
                    fluidState_.setRs(min(RsMax, RsSat));
                }
                else {
                    fluidState_.setRs(0.0);
                }
            }

            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
                fluidState_.setRv(Rv);
            }
            else {
                if (FluidSystem::enableVaporizedOil() ) { // Add Sg > 0? i.e. if only water set rv = 0)
                    const Evaluation& RvSat = enableExtbo ? asImp_().rv() :
                        FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                                gasPhaseIdx,
                                                                pvtRegionIdx,
                                                                SoMax);
                    fluidState_.setRv(min(RvMax, RvSat));
                }
                else {
                    fluidState_.setRv(0.0);
                }
            }
        }

        if constexpr (enableVapwat) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rvw) {
                const auto& Rvw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
                fluidState_.setRvw(Rvw);
            }
            else {
                if (FluidSystem::enableVaporizedWater()) { // Add Sg > 0? i.e. if only water set rv = 0)
                    const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                                                        gasPhaseIdx,
                                                                                        pvtRegionIdx);
                    fluidState_.setRvw(RvwSat);
                }
            }
        }

        if constexpr (enableDisgasInWater) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rsw) {
                const auto& Rsw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
                fluidState_.setRsw(Rsw);
            }
            else {
                if (FluidSystem::enableDissolvedGasInWater()) {
                    const Evaluation& RswSat = FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                                                       waterPhaseIdx,
                                                                                       pvtRegionIdx);
                    fluidState_.setRsw(min(RswMax, RswSat));
                }
            }
        }
    }

    void updateMobilityAndInvB()
    {
        OPM_TIMEBLOCK_LOCAL(updateMobilityAndInvB, Subsystem::PvtProps);
        const unsigned pvtRegionIdx = fluidState_.pvtRegionIndex();

        // compute the phase densities and transform the phase permeabilities into mobilities
        int nmobilities = 1;
        constexpr int max_nmobilities = 4;
        std::array<std::array<Evaluation, numPhases>*, max_nmobilities> mobilities = { &mobility_};
        if (dirMob_) {
            for (int i = 0; i < 3; ++i) {
                mobilities[nmobilities] = &(dirMob_->getArray(i));
                ++nmobilities;
            }
        }
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const auto [b, mu] = FluidSystem::inverseFormationVolumeFactorAndViscosity(fluidState_, phaseIdx, pvtRegionIdx);
            fluidState_.setInvB(phaseIdx, b);
            for (int i = 0; i < nmobilities; ++i) {
                if (enableExtbo && phaseIdx == oilPhaseIdx) {
                    (*mobilities[i])[phaseIdx] /= asImp_().oilViscosity();
                }
                else if (enableExtbo && phaseIdx == gasPhaseIdx) {
                    (*mobilities[i])[phaseIdx] /= asImp_().gasViscosity();
                }
                else {
                    (*mobilities[i])[phaseIdx] /= mu;
                }
            }
        }
        Valgrind::CheckDefined(mobility_);
    }

    void updatePhaseDensities()
    {
        const unsigned pvtRegionIdx = fluidState_.pvtRegionIndex();

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGasInWater()) {
                rho += fluidState_.invB(waterPhaseIdx) *
                       fluidState_.Rsw() *
                       FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            rho = fluidState_.invB(gasPhaseIdx);
            rho *= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableVaporizedOil()) {
                rho += fluidState_.invB(gasPhaseIdx) *
                       fluidState_.Rv() *
                       FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
            if (FluidSystem::enableVaporizedWater()) {
                rho += fluidState_.invB(gasPhaseIdx) *
                       fluidState_.Rvw() *
                       FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(gasPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {
                rho += fluidState_.invB(oilPhaseIdx) *
                       fluidState_.Rs() *
                       FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }
    }

    void updatePorosity(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        // Retrieve the reference porosity from the problem.
        referencePorosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        // Account for other effects.
        this->updatePorosityImpl(problem, priVars, globalSpaceIdx, timeIdx);
    }

    void updatePorosity(const Problem& problem, const PrimaryVariables& priVars, const unsigned globalSpaceIdx, const unsigned timeIdx)
    {
        // Retrieve the reference porosity from the problem.
        referencePorosity_ = problem.porosity(globalSpaceIdx, timeIdx);
        // Account for other effects.
        this->updatePorosityImpl(problem, priVars, globalSpaceIdx, timeIdx);
    }

    void updatePorosityImpl(const Problem& problem, const PrimaryVariables& priVars, const unsigned globalSpaceIdx, const unsigned timeIdx)
    {
        const auto& linearizationType = problem.model().linearizer().getLinearizationType();

        // Start from the reference porosity.
        porosity_ = referencePorosity_;

        // the porosity must be modified by the compressibility of the
        // rock...
        const Scalar rockCompressibility = problem.rockCompressibility(globalSpaceIdx);
        if (rockCompressibility > 0.0) {
            const Scalar rockRefPressure = problem.rockReferencePressure(globalSpaceIdx);
            Evaluation x;
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                x = rockCompressibility * (fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            }
            else if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                x = rockCompressibility * (fluidState_.pressure(waterPhaseIdx) - rockRefPressure);
            }
            else {
                x = rockCompressibility * (fluidState_.pressure(gasPhaseIdx) - rockRefPressure);
            }
            porosity_ *= 1.0 + x + 0.5 * x * x;
        }

        // deal with water induced rock compaction
        porosity_ *= problem.template rockCompPoroMultiplier<Evaluation>(*this, globalSpaceIdx);

        // deal with bioeffects (minimum porosity of 1e-8 to prevent numerical issues)
        if constexpr (enableBioeffects) {
            const Evaluation biofilm_ = priVars.makeEvaluation(Indices::biofilmVolumeFractionIdx,
                                                               timeIdx, linearizationType);
            Evaluation calcite_ = 0.0;
            if constexpr (enableMICP) {
                calcite_ = priVars.makeEvaluation(Indices::calciteVolumeFractionIdx, timeIdx, linearizationType);
            }
            porosity_ -= min(biofilm_ + calcite_, referencePorosity_ - 1e-8);
        }

        // deal with salt-precipitation
        if (enableSaltPrecipitation && priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
            const Evaluation Sp = priVars.makeEvaluation(Indices::saltConcentrationIdx, timeIdx);
            porosity_ *= (1.0 - Sp);
        }
    }

    void assertFiniteMembers()
    {
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            assert(isfinite(fluidState_.density(phaseIdx)));
            assert(isfinite(fluidState_.saturation(phaseIdx)));
            assert(isfinite(fluidState_.temperature(phaseIdx)));
            assert(isfinite(fluidState_.pressure(phaseIdx)));
            assert(isfinite(fluidState_.invB(phaseIdx)));
        }
        assert(isfinite(fluidState_.Rs()));
        assert(isfinite(fluidState_.Rv()));
    }

    /*!
     * \copydoc IntensiveQuantities::update
     */
    template <class ...Args>
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);

        updateCommonPart<Args...>(problem, priVars, globalSpaceIdx, timeIdx);

        updatePorosity(elemCtx, dofIdx, timeIdx);

        // Below: things I want to move to elemCtx-less versions but have not done yet.

        if constexpr (enableSolvent) {
            asImp_().solventPvtUpdate_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (enableExtbo) {
            asImp_().zPvtUpdate_();
        }
        if constexpr (enablePolymer) {
            asImp_().polymerPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal ||
                      energyModuleType == EnergyModules::SequentialImplicitThermal) {
            asImp_().updateEnergyQuantities_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (enableFoam) {
            asImp_().foamPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (enableBioeffects) {
            asImp_().bioeffectsPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (enableBrine) {
            asImp_().saltPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        }
        if constexpr (enableConvectiveMixing) {
            // The ifs are here is to avoid extra calculations for
            // cases with dry runs and without CO2STORE and DRSDTCON.
            if (!problem.simulator().vanguard().eclState().getIOConfig().initOnly()) {
                if (problem.simulator().vanguard().eclState().runspec().co2Storage()) {
                    if (problem.drsdtconIsActive(globalSpaceIdx, problem.simulator().episodeIndex())) {
                        asImp_().updateSaturatedDissolutionFactor_();
                    }
                }
            }
        }

        // update the quantities which are required by the chosen
        // velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        if constexpr (enableDiffusion) {
            DiffusionIntensiveQuantities::update_(fluidState_, priVars.pvtRegionIndex(), elemCtx, dofIdx, timeIdx);
        }

        // update the dispersion specific quantities of the intensive quantities
        if constexpr (enableDispersion) {
            DispersionIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);
        }
    }

    template <class ...Args>
    void update(const Problem& problem, const PrimaryVariables& priVars, const unsigned globalSpaceIdx, const unsigned timeIdx)
    {
        // This is the version of update() that does not use any ElementContext.
        // It is limited by some modules that are not yet adapted to that.
        static_assert(!enableSolvent);
        static_assert(!enableExtbo);
        static_assert(!enablePolymer);
        static_assert(!enableFoam);
        static_assert(!enableMICP);
        static_assert(!enableBrine);
        static_assert(!enableDiffusion);
        static_assert(!enableDispersion);

        this->extrusionFactor_ = 1.0;// to avoid fixing parent update
        updateCommonPart<Args...>(problem, priVars, globalSpaceIdx, timeIdx);
        // Porosity requires separate calls so this can be instantiated with ReservoirProblem from the examples/ directory.
        updatePorosity(problem, priVars, globalSpaceIdx, timeIdx);

        // TODO: Here we should do the parts for solvent etc. at the bottom of the other update() function.
    }

    // This function updated the parts that are common to the IntensiveQuantities regardless of extensions used.
    template <class ...Args>
    void updateCommonPart(const Problem& problem, const PrimaryVariables& priVars, const unsigned globalSpaceIdx, const unsigned timeIdx)
    {
        OPM_TIMEBLOCK_LOCAL(blackoilIntensiveQuanititiesUpdate, Subsystem::SatProps | Subsystem::PvtProps);

        const auto& linearizationType = problem.model().linearizer().getLinearizationType();
        const unsigned pvtRegionIdx = priVars.pvtRegionIndex();

        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        updateTempSalt(problem, priVars, globalSpaceIdx, timeIdx, linearizationType);
        updateSaturations(priVars, timeIdx, linearizationType);
        updateRelpermAndPressures<Args...>(problem, priVars, globalSpaceIdx, timeIdx, linearizationType);

        // update extBO parameters
        if constexpr (enableExtbo) {
            asImp_().zFractionUpdate_(priVars, timeIdx);
        }

        updateRsRvRsw(problem, priVars, globalSpaceIdx, timeIdx);
        updateMobilityAndInvB();
        updatePhaseDensities();

        rockCompTransMultiplier_ = problem.template rockCompTransMultiplier<Evaluation>(*this, globalSpaceIdx);

#ifndef NDEBUG
        assertFiniteMembers();
#endif
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    const Evaluation& mobility(unsigned phaseIdx, FaceDir::DirEnum facedir) const
    {
        using Dir = FaceDir::DirEnum;
        if (dirMob_) {
            switch (facedir) {
                case Dir::XMinus:
                case Dir::XPlus:
                    return dirMob_->getArray(0)[phaseIdx];
                case Dir::YMinus:
                case Dir::YPlus:
                    return dirMob_->getArray(1)[phaseIdx];
                case Dir::ZMinus:
                case Dir::ZPlus:
                    return dirMob_->getArray(2)[phaseIdx];
                default:
                    throw std::runtime_error("Unexpected face direction");
            }
        }
        else {
            return mobility_[phaseIdx];
        }
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * The pressure-dependent transmissibility multiplier due to rock compressibility.
     */
    const Evaluation& rockCompTransMultiplier() const
    { return rockCompTransMultiplier_; }

    /*!
     * \brief Returns the index of the PVT region used to calculate the thermodynamic
     *        quantities.
     *
     * This allows to specify different Pressure-Volume-Temperature (PVT) relations in
     * different parts of the spatial domain.
     */
    auto pvtRegionIndex() const -> decltype(std::declval<FluidState>().pvtRegionIndex())
    { return fluidState_.pvtRegionIndex(); }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Evaluation relativePermeability(unsigned phaseIdx) const
    {
        // warning: slow
        return fluidState_.viscosity(phaseIdx) * mobility(phaseIdx);
    }

    /*!
     * \brief Returns the porosity of the rock at reference conditions.
     *
     * I.e., the porosity of rock which is not perturbed by pressure and temperature
     * changes.
     */
    Scalar referencePorosity() const
    { return referencePorosity_; }

    const Evaluation& permFactor() const
    {
        if constexpr (enableBioeffects) {
            return BioeffectsIntQua::permFactor();
        }
        else if constexpr (enableSaltPrecipitation) {
            return BrineIntQua::permFactor();
        }
        else {
            throw std::logic_error("permFactor() called but salt precipitation or bioeffects are disabled");
        }
    }

    /*!
     * \brief Returns the fluid system used by this intensive quantities.
     */
    OPM_HOST_DEVICE const auto& getFluidSystem() const
    {
        return fluidState_.fluidSystem();
    }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag, enableSolvent>;
    friend BlackOilExtboIntensiveQuantities<TypeTag, enableExtbo>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag, enablePolymer>;
    friend BlackOilEnergyIntensiveQuantities<TypeTag, energyModuleType>;
    friend BlackOilFoamIntensiveQuantities<TypeTag, enableFoam>;
    friend BlackOilBrineIntensiveQuantities<TypeTag, enableBrine>;
    friend BlackOilBioeffectsIntensiveQuantities<TypeTag, enableBioeffects>;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Scalar referencePorosity_;
    Evaluation porosity_;
    Evaluation rockCompTransMultiplier_;
    std::array<Evaluation, numPhases> mobility_;

    // Instead of writing a custom copy constructor and a custom assignment operator just to handle
    // the dirMob_ unique ptr member variable when copying BlackOilIntensiveQuantites (see for example
    // updateIntensitiveQuantities_() in fvbaseelementcontext.hh for a copy example) we write the below
    // custom wrapper class CopyablePtr which wraps the unique ptr and makes it copyable.
    //
    // The advantage of this approach is that we avoid having to call all the base class copy constructors and
    // assignment operators explicitly (which is needed when writing the custom copy constructor and assignment
    // operators) which could become a maintenance burden. For example, when adding a new base class (if that should
    // be needed sometime in the future) to BlackOilIntensiveQuantites we could forget to update the copy
    // constructor and assignment operators.
    //
    // We want each copy of the BlackOilIntensiveQuantites to be unique, (TODO: why?) so we have to make a copy
    // of the unique_ptr each time we copy construct or assign to it from another BlackOilIntensiveQuantites.
    // (On the other hand, if a copy could share the ptr with the original, a shared_ptr could be used instead and the
    // wrapper would not be needed)
    DirectionalMobilityPtr dirMob_;
};

} // namespace Opm

#endif
