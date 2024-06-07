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

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilextbomodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilfoammodules.hh"
#include "blackoilbrinemodules.hh"
#include "blackoilenergymodules.hh"
#include "blackoildiffusionmodule.hh"
#include "blackoildispersionmodule.hh"
#include "blackoilmicpmodules.hh"

#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/directionalmobility.hh>

#include <opm/utility/CopyablePtr.hpp>

#include <dune/common/fmatrix.hh>

#include <cstring>
#include <utility>

#include <fmt/format.h>

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
    , public BlackOilSolventIntensiveQuantities<TypeTag>
    , public BlackOilExtboIntensiveQuantities<TypeTag>
    , public BlackOilPolymerIntensiveQuantities<TypeTag>
    , public BlackOilFoamIntensiveQuantities<TypeTag>
    , public BlackOilBrineIntensiveQuantities<TypeTag>
    , public BlackOilEnergyIntensiveQuantities<TypeTag>
    , public BlackOilMICPIntensiveQuantities<TypeTag>
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
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;
    static constexpr bool waterEnabled = Indices::waterEnabled;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;

    using Toolbox = MathToolbox<Evaluation>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = BlackOilDiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using DispersionIntensiveQuantities = BlackOilDispersionIntensiveQuantities<TypeTag, enableDispersion>;

    using DirectionalMobilityPtr = Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;
    using BrineModule = BlackOilBrineModule<TypeTag>;


public:
    using FluidState = BlackOilFluidState<Evaluation,
                                          FluidSystem,
                                          enableTemperature,
                                          enableEnergy,
                                          compositionSwitchEnabled,
                                          enableVapwat,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          has_disgas_in_water,
                                          Indices::numPhases>;
    using ScalarFluidState = BlackOilFluidState<Scalar,
                                                FluidSystem,
                                                enableTemperature,
                                                enableEnergy,
                                                compositionSwitchEnabled,
                                                enableVapwat,
                                                enableBrine,
                                                enableSaltPrecipitation,
                                                has_disgas_in_water,
                                                Indices::numPhases>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    BlackOilIntensiveQuantities()
    {
        if (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
        }        
        if (enableVapwat) { 
            fluidState_.setRvw(0.0);
        }
        if (has_disgas_in_water) {
            fluidState_.setRsw(0.0);
        }
    }
    BlackOilIntensiveQuantities(const BlackOilIntensiveQuantities& other) = default;

    BlackOilIntensiveQuantities& operator=(const BlackOilIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        OPM_TIMEBLOCK_LOCAL(blackoilIntensiveQuanititiesUpdate);
        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& linearizationType = problem.model().linearizer().getLinearizationType();
        unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        Scalar RvMax = FluidSystem::enableVaporizedOil()
            ? problem.maxOilVaporizationFactor(timeIdx, globalSpaceIdx)
            : 0.0;
        Scalar RsMax = FluidSystem::enableDissolvedGas()
            ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
            : 0.0;
        Scalar RswMax = FluidSystem::enableDissolvedGasInWater()
            ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
            : 0.0;

        asImp_().updateTemperature_(elemCtx, dofIdx, timeIdx);

        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        asImp_().updateSaltConcentration_(elemCtx, dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                Sw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            } else if(priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rsw ||
                      priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Disabled) {
                      // water is enabled but is not a primary variable i.e. one component/phase case
                      // or two-phase water + gas with only water present
                Sw = 1.0;
            } // else i.e. for MeaningWater() = Rvw, Sw is still 0.0;
        }
        Evaluation Sg = 0.0;
        if constexpr (gasEnabled) {
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                Sg = 1.0 - Sw;
            } else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Disabled) {
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
            if(priVars.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
                } else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    Sg -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
                }
            }
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            fluidState_.setSaturation(waterPhaseIdx, Sw);

        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            fluidState_.setSaturation(gasPhaseIdx, Sg);

        if (FluidSystem::phaseIsActive(oilPhaseIdx))
            fluidState_.setSaturation(oilPhaseIdx, So);

        asImp_().solventPreSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        // now we compute all phase pressures
        std::array<Evaluation, numPhases> pC;
        const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        problem.updateRelperms(mobility_, dirMob_, fluidState_, globalSpaceIdx);

        // scaling the capillary pressure due to salt precipitation 
        if (BrineModule::hasPcfactTables() && priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
            unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);
            const Evaluation Sp = priVars.makeEvaluation(Indices::saltConcentrationIdx, timeIdx);
            const Evaluation porosityFactor  = min(1.0 - Sp, 1.0); //phi/phi_0
            const auto& pcfactTable = BrineModule::pcfactTable(satnumRegionIdx);
            const Evaluation pcFactor = pcfactTable.eval(porosityFactor, /*extrapolation=*/true);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    pC[phaseIdx] *= pcFactor;
                }
        }

        // oil is the reference phase for pressure
        if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
        } else if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pw) {
            const Evaluation& pw = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pw + (pC[phaseIdx] - pC[waterPhaseIdx]));
        } else {
            assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
        }

        // update the Saturation functions for the blackoil solvent module.
        asImp_().solventPostSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        // update extBO parameters
        asImp_().zFractionUpdate_(elemCtx, dofIdx, timeIdx);

        Evaluation SoMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            SoMax = max(fluidState_.saturation(oilPhaseIdx),
                        problem.maxOilSaturation(globalSpaceIdx));
        
        
        // take the meaning of the switching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
            const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRs(Rs);
        } else {
            if (FluidSystem::enableDissolvedGas()) { // Add So > 0? i.e. if only water set rs = 0)
                const Evaluation& RsSat = enableExtbo ? asImp_().rs() :
                FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                    oilPhaseIdx,
                                                    pvtRegionIdx,
                                                    SoMax);
                fluidState_.setRs(min(RsMax, RsSat));
            }
            else if constexpr (compositionSwitchEnabled)
                fluidState_.setRs(0.0);
        }
        }
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
            const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRv(Rv);
        } else {
            if (FluidSystem::enableVaporizedOil() ) { // Add Sg > 0? i.e. if only water set rv = 0)
                const Evaluation& RvSat = enableExtbo ? asImp_().rv() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRv(min(RvMax, RvSat));
            }
            else if constexpr (compositionSwitchEnabled)
                fluidState_.setRv(0.0);
        }

        if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rvw) {
            const auto& Rvw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            fluidState_.setRvw(Rvw);
        } else {
            if (FluidSystem::enableVaporizedWater()) { // Add Sg > 0? i.e. if only water set rv = 0)
                const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRvw(RvwSat);
            }
        }

        if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rsw) {
            const auto& Rsw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            fluidState_.setRsw(Rsw);
        } else {
            if (FluidSystem::enableDissolvedGasInWater()) {
                const Evaluation& RswSat = FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            waterPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRsw(min(RswMax, RswSat));
            }
        }

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.setRegionIndex(pvtRegionIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            paramCache.setMaxOilSat(SoMax);
        }
        paramCache.updateAll(fluidState_);

        // compute the phase densities and transform the phase permeabilities into mobilities
        int nmobilities = 1;
        std::vector<std::array<Evaluation,numPhases>*> mobilities = {&mobility_};
        if (dirMob_) {
            for (int i=0; i<3; i++) {
                nmobilities += 1;
                mobilities.push_back(&(dirMob_->getArray(i)));
            }
        }
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;
            const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState_, phaseIdx, pvtRegionIdx);
            fluidState_.setInvB(phaseIdx, b);
            const auto& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            for (int i = 0; i<nmobilities; i++) {
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

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGasInWater()) {
                const auto& oilVaporizationControl = problem.simulator().vanguard().schedule()[problem.episodeIndex()].oilvap();
				if(!oilVaporizationControl.drsdtConvective()) {
                    rho +=
                        fluidState_.invB(waterPhaseIdx) *
                        fluidState_.Rsw() *
                        FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                }
            }
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            rho = fluidState_.invB(gasPhaseIdx);
            rho *= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableVaporizedOil()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rv() *
                    FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
            if (FluidSystem::enableVaporizedWater()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rvw() *
                    FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(gasPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {                
                const auto& oilVaporizationControl = problem.simulator().vanguard().schedule()[problem.episodeIndex()].oilvap();
				if(!oilVaporizationControl.drsdtConvective()) {
                    rho +=
                        fluidState_.invB(oilPhaseIdx) *
                        fluidState_.Rs() *
                        FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
                }
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }

        // retrieve the porosity from the problem
        referencePorosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        porosity_ = referencePorosity_;

        // the porosity must be modified by the compressibility of the
        // rock...
        Scalar rockCompressibility = problem.rockCompressibility(globalSpaceIdx);
        if (rockCompressibility > 0.0) {
            Scalar rockRefPressure = problem.rockReferencePressure(globalSpaceIdx);
            Evaluation x;
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                x = rockCompressibility*(fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            } else if (FluidSystem::phaseIsActive(waterPhaseIdx)){
                x = rockCompressibility*(fluidState_.pressure(waterPhaseIdx) - rockRefPressure);
            } else {
                x = rockCompressibility*(fluidState_.pressure(gasPhaseIdx) - rockRefPressure);
            }
            porosity_ *= 1.0 + x + 0.5*x*x;
        }

        // deal with water induced rock compaction
        porosity_ *= problem.template rockCompPoroMultiplier<Evaluation>(*this, globalSpaceIdx);

        // the MICP processes change the porosity
        if constexpr (enableMICP){
          Evaluation biofilm_ = priVars.makeEvaluation(Indices::biofilmConcentrationIdx, timeIdx, linearizationType);
          Evaluation calcite_ = priVars.makeEvaluation(Indices::calciteConcentrationIdx, timeIdx, linearizationType);
          porosity_ += - biofilm_ - calcite_;
        }

        // deal with salt-precipitation
        if (enableSaltPrecipitation && priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
            Evaluation Sp = priVars.makeEvaluation(Indices::saltConcentrationIdx, timeIdx);
            porosity_ *= (1.0 - Sp);
        }

        rockCompTransMultiplier_ = problem.template rockCompTransMultiplier<Evaluation>(*this, globalSpaceIdx);

        asImp_().solventPvtUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().zPvtUpdate_();
        asImp_().polymerPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().updateEnergyQuantities_(elemCtx, dofIdx, timeIdx, paramCache);
        asImp_().foamPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().MICPPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().saltPropertiesUpdate_(elemCtx, dofIdx, timeIdx);

        // update the quantities which are required by the chosen
        // velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the dispersion specific quantities of the intensive quantities
        DispersionIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

#ifndef NDEBUG
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            assert(isfinite(fluidState_.density(phaseIdx)));
            assert(isfinite(fluidState_.saturation(phaseIdx)));
            assert(isfinite(fluidState_.temperature(phaseIdx)));
            assert(isfinite(fluidState_.pressure(phaseIdx)));
            assert(isfinite(fluidState_.invB(phaseIdx)));
        }
        assert(isfinite(fluidState_.Rs()));
        assert(isfinite(fluidState_.Rv()));

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
            switch(facedir) {
                case Dir::XMinus:
                case Dir::XPlus:
                    return dirMob_->mobilityX_[phaseIdx];
                case Dir::YMinus:
                case Dir::YPlus:
                    return dirMob_->mobilityY_[phaseIdx];
                case Dir::ZMinus:
                case Dir::ZPlus:
                    return dirMob_->mobilityZ_[phaseIdx];
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
    auto pvtRegionIndex() const
        -> decltype(std::declval<FluidState>().pvtRegionIndex())
    { return fluidState_.pvtRegionIndex(); }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Evaluation relativePermeability(unsigned phaseIdx) const
    {
        // warning: slow
        return fluidState_.viscosity(phaseIdx)*mobility(phaseIdx);
    }

    /*!
     * \brief Returns the porosity of the rock at reference conditions.
     *
     * I.e., the porosity of rock which is not perturbed by pressure and temperature
     * changes.
     */
    Scalar referencePorosity() const
    { return referencePorosity_; }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag>;
    friend BlackOilExtboIntensiveQuantities<TypeTag>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag>;
    friend BlackOilEnergyIntensiveQuantities<TypeTag>;
    friend BlackOilFoamIntensiveQuantities<TypeTag>;
    friend BlackOilBrineIntensiveQuantities<TypeTag>;
    friend BlackOilMICPIntensiveQuantities<TypeTag>;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Scalar referencePorosity_;
    Evaluation porosity_;
    Evaluation rockCompTransMultiplier_;
    std::array<Evaluation,numPhases> mobility_;

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
