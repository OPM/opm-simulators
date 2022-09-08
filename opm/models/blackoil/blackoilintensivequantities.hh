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
#include "blackoilmicpmodules.hh"
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

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
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableEvaporation = getPropValue<TypeTag, Properties::EnableEvaporation>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
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

    struct DirectionalMobility {
        using array_type = std::array<Evaluation,numPhases>;
        DirectionalMobility(const array_type& mX, const array_type& mY, const array_type& mZ) 
            : mobilityX_{mX}, mobilityY_{mY}, mobilityZ_{mZ} {}
        DirectionalMobility() : mobilityX_{}, mobilityY_{}, mobilityZ_{} {}
        array_type& getArray(int index) {
            switch(index) {
                case 0:
                    return mobilityX_;
                case 1:
                    return mobilityY_;
                case 2:
                    return mobilityZ_;
                default:
                    throw std::runtime_error("Unexpected mobility array index");
            }
        }
        array_type mobilityX_;
        array_type mobilityY_;
        array_type mobilityZ_;
    };

public:
    using FluidState = BlackOilFluidState<Evaluation,
                                          FluidSystem,
                                          enableTemperature,
                                          enableEnergy,
                                          compositionSwitchEnabled,
                                          enableEvaporation,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          Indices::numPhases>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    BlackOilIntensiveQuantities()
    {
        if (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
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

        asImp_().updateTemperature_(elemCtx, dofIdx, timeIdx);

        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        asImp_().updateSaltConcentration_(elemCtx, dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            if (priVars.primaryVarsMeaning() == PrimaryVariables::OnePhase_p) {
                Sw = 1.0;
            } else if (priVars.primaryVarsMeaning() != PrimaryVariables::Rvw_po_Sg && priVars.primaryVarsMeaning() != PrimaryVariables::Rvw_pg_Rv) {
                Sw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx);
            }
        }
        Evaluation Sg = 0.0;
        if constexpr (compositionSwitchEnabled)
        {
            if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                // -> threephase case
                assert( priVars.primaryVarsMeaning() != PrimaryVariables::OnePhase_p );
                Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
                // -> gas-water case
                Sg = 1.0 - Sw;

                // deal with solvent
                if (enableSolvent)
                    Sg -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);

            } else if (priVars.primaryVarsMeaning() == PrimaryVariables::Rvw_po_Sg) {
                 // -> oil-gas case
                 Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaning() == PrimaryVariables::Rvw_pg_Rv) {
                 // -> gas case
                 Sg = 1.0;
            } else
            {
                assert(priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs);
                // -> oil-water case
                Sg = 0.0;
            }
        }
        if constexpr (gasEnabled && waterEnabled && !oilEnabled) {
            Sg = 1.0 - Sw;
        }

        Valgrind::CheckDefined(Sg);
        Valgrind::CheckDefined(Sw);

        Evaluation So = 1.0 - Sw - Sg;

        // deal with solvent
        if constexpr (enableSolvent)
            So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);

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
        updateRelperms(problem, materialParams, globalSpaceIdx);

        // oil is the reference phase for pressure
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv || priVars.primaryVarsMeaning() == PrimaryVariables::Rvw_pg_Rv) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
        }

        else {
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
        }

        // take the meaning of the switching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
            // in the threephase case, gas and oil phases are potentially present, i.e.,
            // we use the compositions of the gas-saturated oil and oil-saturated gas.
            if (FluidSystem::enableDissolvedGas()) {
                const Evaluation& RsSat = enableExtbo ? asImp_().rs() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRs(min(RsMax, RsSat));
            }
            else if constexpr (compositionSwitchEnabled)
                fluidState_.setRs(0.0);

            if (FluidSystem::enableVaporizedOil()) {
                const Evaluation& RvSat = enableExtbo ? asImp_().rv() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRv(min(RvMax, RvSat));
            }
            else if constexpr (compositionSwitchEnabled)
                fluidState_.setRv(0.0);

            if (FluidSystem::enableVaporizedWater()) {
                const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRvw(RvwSat);
            }
        }
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Rvw_po_Sg) {
            // The switching variable is the water-gas ratio Rvw
            const auto& Rvw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx);
            fluidState_.setRvw(Rvw);

            if (FluidSystem::enableVaporizedOil()) {
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
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Rvw_pg_Rv) {
            // The switching variable is the water-gas ratio Rvw
            const auto& Rvw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx);
            fluidState_.setRvw(Rvw);

            const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRv(Rv);

            if (FluidSystem::enableDissolvedGas()) {
                // the oil phase is not present, but we need to compute its "composition" for
                // the gravity correction anyway
                const auto& RsSat = enableExtbo ? asImp_().rs() :
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
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
            // if the switching variable is the mole fraction of the gas component in the
            // oil phase, we can directly set the composition of the oil phase
            const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRs(min(RsMax, Rs));

            if (FluidSystem::enableVaporizedOil()) {
                // the gas phase is not present, but we need to compute its "composition"
                // for the gravity correction anyway
                const auto& RvSat = enableExtbo ? asImp_().rv() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRv(min(RvMax, RvSat));
            }
            else
                fluidState_.setRv(0.0);

            if (FluidSystem::enableVaporizedWater()) {
                const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRvw(RvwSat);
            }
        }
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
            const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRv(Rv);

            if (FluidSystem::enableDissolvedGas()) {
                // the oil phase is not present, but we need to compute its "composition" for
                // the gravity correction anyway
                const auto& RsSat = enableExtbo ? asImp_().rs() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRs(min(RsMax, RsSat));
            } else {
                fluidState_.setRs(0.0);
            }

            if (FluidSystem::enableVaporizedWater()) {
                const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRvw(RvwSat);
            }
        } else {
            assert(priVars.primaryVarsMeaning() == PrimaryVariables::OnePhase_p);
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
                rho +=
                    fluidState_.invB(oilPhaseIdx) *
                    fluidState_.Rs() *
                    FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
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
        if (enableSaltPrecipitation && priVars.primaryVarsMeaningBrine() == PrimaryVariables::Sp) {
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
                case Dir::XPlus:
                    return dirMob_->mobilityX_[phaseIdx];
                case Dir::YPlus:
                    return dirMob_->mobilityY_[phaseIdx];
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
     * different parts of the spatial domain. Note that this concept should be seen as a
     * work-around of the fact that the black-oil model does not capture the
     * thermodynamics well enough. (Because there is, err, only a single real world with
     * in which all substances follow the same physical laws and hence the same
     * thermodynamics.) Anyway: Since the ECL file format uses multiple PVT regions, we
     * support it as well in our black-oil model. (Note that, if it is not explicitly
     * specified, the PVT region index is 0.)
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

    void updateRelperms(const Problem& problem, const MaterialLawParams& materialParams, unsigned globalSpaceIdx)
    {
        // calculate relative permeabilities. note that we store the result into the
        // mobility_ class attribute. the division by the phase viscosity happens later.
        MaterialLaw::relativePermeabilities(mobility_, materialParams, fluidState_);
        Valgrind::CheckDefined(mobility_);
        const auto& materialLawManager = problem.materialLawManager();
        if (materialLawManager->hasDirectionalRelperms()) {
            auto satnumIdx = materialLawManager->satnumRegionIdx(globalSpaceIdx);
            using Dir = FaceDir::DirEnum;
            constexpr int ndim = 3;
            dirMob_ = std::make_unique<DirectionalMobility>();
            Dir facedirs[ndim] = {Dir::XPlus, Dir::YPlus, Dir::ZPlus};
            for (int i = 0; i<ndim; i++) {
                auto krnumSatIdx = materialLawManager->getKrnumSatIdx(globalSpaceIdx, facedirs[i]);
                auto& mob_array = dirMob_->getArray(i);
                if (krnumSatIdx != satnumIdx) {
                    // This hack is also used by StandardWell_impl.hpp:getMobilityEval() to temporarily use a different
                    // satnum index for a cell
                    const auto& paramsCell = materialLawManager->connectionMaterialLawParams(krnumSatIdx, globalSpaceIdx);
                    MaterialLaw::relativePermeabilities(mob_array, paramsCell, fluidState_);
                    // reset the cell's satnum index back to the original
                    materialLawManager->connectionMaterialLawParams(satnumIdx, globalSpaceIdx);
                }
                else {
                    // Copy the default (non-directional dependent) mobility
                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        mob_array[phaseIdx] = mobility_[phaseIdx];
                    }
                }
            }
        }
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Scalar referencePorosity_;
    Evaluation porosity_;
    Evaluation rockCompTransMultiplier_;
    std::array<Evaluation,numPhases> mobility_;
    std::shared_ptr<DirectionalMobility> dirMob_;
};

} // namespace Opm

#endif
