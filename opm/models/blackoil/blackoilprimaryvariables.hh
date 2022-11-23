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
 * \copydoc Opm::BlackOilPrimaryVariables
 */
#ifndef EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH
#define EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilextbomodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilenergymodules.hh"
#include "blackoilfoammodules.hh"
#include "blackoilbrinemodules.hh"
#include "blackoilmicpmodules.hh"

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>

#include <dune/common/fvector.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/common/Valgrind.hpp>

namespace Opm {
template <class TypeTag, bool enableSolvent>
class BlackOilSolventModule;

template <class TypeTag, bool enableExtbo>
class BlackOilExtboModule;

template <class TypeTag, bool enablePolymer>
class BlackOilPolymerModule;

template <class TypeTag, bool enableBrine>
class BlackOilBrineModule;

/*!
 * \ingroup BlackOilModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class BlackOilPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    // number of equations
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };

    // primary variable indices
    enum { waterSaturationIdx = Indices::waterSaturationIdx };
    enum { pressureSwitchIdx = Indices::pressureSwitchIdx };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };
    enum { saltConcentrationIdx  = Indices::saltConcentrationIdx };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;
    static constexpr bool waterEnabled = Indices::waterEnabled;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;

    // phase indices from the fluid system
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    // component indices from the fluid system
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableEvaporation = getPropValue<TypeTag, Properties::EnableEvaporation>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

    using Toolbox = MathToolbox<Evaluation>;
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using SolventModule = BlackOilSolventModule<TypeTag, enableSolvent>;
    using ExtboModule = BlackOilExtboModule<TypeTag, enableExtbo>;
    using PolymerModule = BlackOilPolymerModule<TypeTag, enablePolymer>;
    using EnergyModule = BlackOilEnergyModule<TypeTag, enableEnergy>;
    using FoamModule = BlackOilFoamModule<TypeTag, enableFoam>;
    using BrineModule = BlackOilBrineModule<TypeTag, enableBrine>;
    using MICPModule = BlackOilMICPModule<TypeTag, enableMICP>;

    static_assert(numPhases == 3, "The black-oil model assumes three phases!");
    static_assert(numComponents == 3, "The black-oil model assumes three components!");

public:
    enum PrimaryVarsMeaningWater {
        Sw,  // water saturation
        Rvw, // vaporized water
        W_disabled,
    };

    enum PrimaryVarsMeaningPressure {
        Po, // oil pressure
        Pg, // gas pressure
        Pw, // water pressure
    };
    enum PrimaryVarsMeaningGas {
        Sg, // gas saturation
        Rs, // dissolved gas in oil
        //Rsw, // dissolved gas in water
        Rv, // vapporized oil
        G_disabled,
    };

    enum PrimaryVarsMeaningBrine {
        Cs, // salt concentration
        Sp, // (precipitated) salt saturation
        Brine_disabled,
    };

    BlackOilPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
        pvtRegionIdx_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::SetUndefined(primaryVarsMeaningWater_);
        Valgrind::SetUndefined(primaryVarsMeaningGas_);
        Valgrind::SetUndefined(primaryVarsMeaningPressure_);
        Valgrind::SetUndefined(primaryVarsMeaningBrine_);

        pvtRegionIdx_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables& )
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables& value) = default;

    /*!
     * \brief Set the index of the region which should be used for PVT properties.
     *
     * The concept of PVT regions is a hack to work around the fact that the
     * pseudo-components used by the black oil model (i.e., oil, gas and water) change
     * their composition within the spatial domain. We implement them because, the ECL
     * file format mandates them.
     */
    void setPvtRegionIndex(unsigned value)
    { pvtRegionIdx_ = static_cast<unsigned short>(value); }

    /*!
     * \brief Return the index of the region which should be used for PVT properties.
     */
    unsigned pvtRegionIndex() const
    { return pvtRegionIdx_; }

    /*!
     * \brief Return the interpretation which should be applied to the switching primary
     *        variables.
     */
    PrimaryVarsMeaningWater primaryVarsMeaningWater() const
    { return primaryVarsMeaningWater_; }

    /*!
     * \brief Set the interpretation which should be applied to the switching primary
     *        variables.
     */
    void setPrimaryVarsMeaningWater(PrimaryVarsMeaningWater newMeaning)
    { primaryVarsMeaningWater_ = newMeaning; }

     /*!
     * \brief Return the interpretation which should be applied to the switching primary
     *        variables.
     */
    PrimaryVarsMeaningPressure primaryVarsMeaningPressure() const
    { return primaryVarsMeaningPressure_; }

    /*!
     * \brief Set the interpretation which should be applied to the switching primary
     *        variables.
     */
    void setPrimaryVarsMeaningPressure(PrimaryVarsMeaningPressure newMeaning)
    { primaryVarsMeaningPressure_ = newMeaning; }

     /*!
     * \brief Return the interpretation which should be applied to the switching primary
     *        variables.
     */
    PrimaryVarsMeaningGas primaryVarsMeaningGas() const
    { return primaryVarsMeaningGas_; }

    /*!
     * \brief Set the interpretation which should be applied to the switching primary
     *        variables.
     */
    void setPrimaryVarsMeaningGas(PrimaryVarsMeaningGas newMeaning)
    { primaryVarsMeaningGas_ = newMeaning; }

    PrimaryVarsMeaningBrine primaryVarsMeaningBrine() const
    { return primaryVarsMeaningBrine_; }

    /*!
     * \brief Set the interpretation which should be applied to the switching primary
     *        variables.
     */

    void setPrimaryVarsMeaningBrine(PrimaryVarsMeaningBrine newMeaning)
    { primaryVarsMeaningBrine_ = newMeaning; }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium = false)
    {
        using ConstEvaluation = typename std::remove_reference<typename FluidState::Scalar>::type;
        using FsEvaluation = typename std::remove_const<ConstEvaluation>::type;
        using FsToolbox = MathToolbox<FsEvaluation>;

#ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            Valgrind::CheckDefined(fluidState.temperature(0));
            Valgrind::CheckDefined(fluidState.temperature(phaseIdx));

            assert(fluidState.temperature(0) == fluidState.temperature(phaseIdx));
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // If your compiler bails out here, you're probably not using a suitable black
        // oil fluid system.
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.setRegionIndex(pvtRegionIdx_);
        paramCache.setMaxOilSat(FsToolbox::value(fluidState.saturation(oilPhaseIdx)));

        // create a mutable fluid state with well defined densities based on the input
        using NcpFlash = NcpFlash<Scalar, FluidSystem>;
        using FlashFluidState = CompositionalFluidState<Scalar, FluidSystem>;
        FlashFluidState fsFlash;
        fsFlash.setTemperature(FsToolbox::value(fluidState.temperature(/*phaseIdx=*/0)));
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fsFlash.setPressure(phaseIdx, FsToolbox::value(fluidState.pressure(phaseIdx)));
            fsFlash.setSaturation(phaseIdx, FsToolbox::value(fluidState.saturation(phaseIdx)));
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                fsFlash.setMoleFraction(phaseIdx, compIdx, FsToolbox::value(fluidState.moleFraction(phaseIdx, compIdx)));
        }

        paramCache.updateAll(fsFlash);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            Scalar rho = FluidSystem::template density<FlashFluidState, Scalar>(fsFlash, paramCache, phaseIdx);
            fsFlash.setDensity(phaseIdx, rho);
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                globalMolarities[compIdx] +=
                    fsFlash.saturation(phaseIdx) * fsFlash.molarity(phaseIdx, compIdx);
            }
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium

        // run the flash calculation
        NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        using ConstEvaluation = typename std::remove_reference<typename FluidState::Scalar>::type;
        using FsEvaluation = typename std::remove_const<ConstEvaluation>::type;
        using FsToolbox = MathToolbox<FsEvaluation>;

        bool gasPresent = FluidSystem::phaseIsActive(gasPhaseIdx)?(fluidState.saturation(gasPhaseIdx) > 0.0):false;
        bool oilPresent = FluidSystem::phaseIsActive(oilPhaseIdx)?(fluidState.saturation(oilPhaseIdx) > 0.0):false;
        bool waterPresent = FluidSystem::phaseIsActive(waterPhaseIdx)?(fluidState.saturation(waterPhaseIdx) > 0.0):false;
        const auto& saltSaturation = BlackOil::getSaltSaturation_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
        bool precipitatedSaltPresent = enableSaltPrecipitation?(saltSaturation > 0.0):false;
        bool oneActivePhases = FluidSystem::numActivePhases() == 1;
        // deal with the primary variables for the energy extension
        EnergyModule::assignPrimaryVars(*this, fluidState);

        // determine the meaning of the primary variables
        if (gasPresent && FluidSystem::enableVaporizedOil() && !oilPresent){
            primaryVarsMeaningPressure_ = Pg;
        } else if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            primaryVarsMeaningPressure_ = Po;
        } else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            primaryVarsMeaningPressure_ = Pg;
        } else {
            assert(FluidSystem::phaseIsActive(waterPhaseIdx));
            primaryVarsMeaningPressure_ = Pw;
        }

        // determine the meaning of the primary variables
        if ( waterPresent && gasPresent ){
            primaryVarsMeaningWater_ = Sw;
        } else if (gasPresent && FluidSystem::enableVaporizedWater()) {
            primaryVarsMeaningWater_ = Rvw;
        } else if (FluidSystem::phaseIsActive(waterPhaseIdx) && !oneActivePhases) {
            primaryVarsMeaningWater_ = Sw;
        } else {
            primaryVarsMeaningWater_ = W_disabled;
        }

        // determine the meaning of the primary variables
        if ( gasPresent && oilPresent ) {
            primaryVarsMeaningGas_ = Sg;
        } else if (oilPresent && FluidSystem::enableDissolvedGas()) {
            primaryVarsMeaningGas_ = Rs;
        } else if (gasPresent && FluidSystem::enableVaporizedOil()){
            primaryVarsMeaningGas_ = Rv;
        //} else if (waterPresent && FluidSystem::enableDissolvedGasInWater()) {
        //    primaryVarsMeaningGas_ = Rsw;
        } else if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx)) {
            primaryVarsMeaningGas_ = Sg;
        } else {
            primaryVarsMeaningGas_ = G_disabled;
        }

        if constexpr (enableSaltPrecipitation){
            if (precipitatedSaltPresent)
                primaryVarsMeaningBrine_ = Sp;
            else
                primaryVarsMeaningBrine_ = Cs;
        } else {
            primaryVarsMeaningBrine_ = Brine_disabled;
        }

        // assign the actual primary variables
        switch(primaryVarsMeaningPressure()) {
            case Po:
                (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
                break;
            case Pg:
                (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(gasPhaseIdx));
                break;
            case Pw:
                (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(waterPhaseIdx));
                break;
            default:
                throw std::logic_error("No valid primary variable selected for pressure");
        }
        switch(primaryVarsMeaningWater()) {
            case Sw:
            {
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
                break;
            }
            case Rvw:
            {
                const auto& rvw = BlackOil::getRvw_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
                (*this)[waterSaturationIdx] = rvw;
                break;
            }
            case W_disabled:
            {
                break;
            }
            default:
                throw std::logic_error("No valid primary variable selected for water");
        }
        switch(primaryVarsMeaningGas()) {
            case Sg:
            {
                (*this)[compositionSwitchIdx] = FsToolbox::value(fluidState.saturation(gasPhaseIdx));
                break;
            }
            case Rs:
            {
                const auto& rs = BlackOil::getRs_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
                (*this)[compositionSwitchIdx] = rs;
                break;
            }
            case Rv:
            {
                const auto& rv = BlackOil::getRv_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
                (*this)[compositionSwitchIdx] = rv;
                break;
            }
            //case Rsw:
            //{
                //const auto& Rsw = BlackOil::getRsw_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
                //(*this)[waterSaturationIdx] = Rsw;
            //    break;
            //}
            case G_disabled:
            {
                break;
            }
            default:
                throw std::logic_error("No valid primary variable selected for composision");
        }
        checkDefined();
    }

    /*!
     * \brief Adapt the interpretation of the switching variables to be physically
     *        meaningful.
     *
     * If the meaning of the primary variables changes, their values are also adapted in a
     * meaningful manner. (e.g. if the gas phase appears and the composition switching
     * variable changes its meaning from the gas dissolution factor Rs to the gas
     * saturation Sg, the value for this variable is set to zero.)
     * A Scalar eps can be passed to make the switching condition more strict.
     * Useful for avoiding ocsilation in the primaryVarsMeaning.
     *
     * \return true Iff the interpretation of one of the switching variables was changed
     */
    bool adaptPrimaryVariables(const Problem& problem, unsigned globalDofIdx, Scalar eps = 0.0)
    {
        static const Scalar thresholdWaterFilledCell = 1.0 - eps;

        // this function accesses quite a few black-oil specific low-level functions
        // directly for better performance (instead of going the canonical way through
        // the IntensiveQuantities). The reason is that most intensive quantities are not
        // required to be able to decide if the primary variables needs to be switched or
        // not, so it would be a waste to compute them.
        if (primaryVarsMeaningWater() == W_disabled && primaryVarsMeaningGas() == G_disabled){
            return false;
        }
        Scalar sw = 0.0;
        Scalar sg = 0.0;
        Scalar saltConcentration = 0.0;
        const Scalar& T = asImp_().temperature_();
        if (primaryVarsMeaningWater() == Sw)
            sw = (*this)[waterSaturationIdx];
        if (primaryVarsMeaningGas() == Sg)
            sg = (*this)[compositionSwitchIdx];

        if (primaryVarsMeaningGas() == G_disabled && gasEnabled)
            sg = 1.0 - sw; // water + gas case

        if constexpr (enableSaltPrecipitation) {
            Scalar saltSolubility = BrineModule::saltSol(pvtRegionIndex());
            if (primaryVarsMeaningBrine() == Sp) {
                saltConcentration = saltSolubility;
                Scalar saltSat = (*this)[saltConcentrationIdx];
                if (saltSat < -eps){ //precipitated salt dissappears
                    setPrimaryVarsMeaningBrine(Cs);
                    (*this)[saltConcentrationIdx] = saltSolubility; //set salt concentration to solubility limit
                }
            }
            else if (primaryVarsMeaningBrine() == Cs) {
                saltConcentration = (*this)[saltConcentrationIdx];
                if (saltConcentration > saltSolubility + eps){ //salt concentration exceeds solubility limit
                    setPrimaryVarsMeaningBrine(Sp);
                    (*this)[saltConcentrationIdx] = 0.0;
                }
            }
        }

        bool changed = false;
        // special case for cells with almost only water
        if (sw >= thresholdWaterFilledCell) {

            // make sure water saturations does not exceed 1.0
            if constexpr (waterEnabled)
                (*this)[Indices::waterSaturationIdx] = 1.0;
            // the hydrocarbon gas saturation is set to 0.0
            if constexpr (compositionSwitchEnabled)
                (*this)[Indices::compositionSwitchIdx] = 0.0;

            //const Scalar& po = (*this)[pressureSwitchIdx];
            changed = primaryVarsMeaningWater() != Sw || primaryVarsMeaningGas() != Sg;
            if(changed) {
                if constexpr (waterEnabled)
                    setPrimaryVarsMeaningWater(Sw);
                if constexpr (compositionSwitchEnabled)
                    setPrimaryVarsMeaningGas(Sg);

                //setPrimaryVarsMeaningPressure(Po);
                // use water pressure?
            }
            return changed;
        }

        switch(primaryVarsMeaningWater()) {
            case Sw:
            {
                if(sw < -eps && sg > eps && FluidSystem::enableVaporizedWater()) {
                    Scalar p = (*this)[pressureSwitchIdx];
                    if(primaryVarsMeaningPressure() == Po) {
                        std::array<Scalar, numPhases> pC = { 0.0 };
                        const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                        Scalar so = 1.0 - sg - solventSaturation_();
                        computeCapillaryPressures_(pC, so, sg + solventSaturation_(), /*sw=*/ 0.0, matParams);
                        p += (pC[gasPhaseIdx] - pC[oilPhaseIdx]);
                    }
                    Scalar rvwSat = FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_,
                                                                                   T,
                                                                                   p,
                                                                                   saltConcentration);
                    setPrimaryVarsMeaningWater(Rvw);
                    (*this)[Indices::waterSaturationIdx] = rvwSat; //primary variable becomes Rvw
                    changed = true;
                    break;
                }
                //if(sg < -eps && sw > eps && FluidSystem::enableDissolvedGasInWater()) {
                //    const Scalar& po = (*this)[pressureSwitchIdx];
                //    std::array<Scalar, numPhases> pC = { 0.0 };
                //    const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                //    Scalar so = 1.0 - sw - solventSaturation_();
                //    computeCapillaryPressures_(pC, so,  /*sg=*/ 0.0, sw, matParams);
                //    Scalar pw = po + (pC[waterPhaseIdx] - pC[oilPhaseIdx]);
                //    Scalar rswSat = FluidSystem::waterPvt().saturatedWaterDissolutionFactor(pvtRegionIdx_,
                //                                                                   T,
                //                                                                   pw,
                //                                                                   saltConcentration);
                //    setPrimaryVarsMeaningWater(Rsw);
                //    (*this)[Indices::waterSaturationIdx] = rswSat; //primary variable becomes Rsw
                //    changed = true;
                //    break;
                //}
                break;
            }
            case Rvw:
            {
                const Scalar& rvw = (*this)[waterSaturationIdx];
                Scalar p = (*this)[pressureSwitchIdx];
                if(primaryVarsMeaningPressure() == Po) {
                    std::array<Scalar, numPhases> pC = { 0.0 };
                    const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                    Scalar so = 1.0 - sg - solventSaturation_();
                    computeCapillaryPressures_(pC, so, sg + solventSaturation_(), /*sw=*/ 0.0, matParams);
                    p += (pC[gasPhaseIdx] - pC[oilPhaseIdx]);
                }
                Scalar rvwSat = FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_,
                                                                                   T,
                                                                                   p,
                                                                                   saltConcentration);
                if (rvw > rvwSat*(1.0 + eps)) {
                    // water phase appears
                    setPrimaryVarsMeaningWater(Sw);
                    (*this)[Indices::waterSaturationIdx] = 0.0; // water saturation
                    changed = true;
                }
                break;
            }
            case W_disabled:
            {
                break;
            }
            default:
                throw std::logic_error("No valid primary variable selected for water");
        }

        switch(primaryVarsMeaningGas()) {
            case Sg:
            {
                Scalar s = 1.0 - sw - solventSaturation_();
                if (sg < -eps && s > 0.0 && FluidSystem::enableDissolvedGas()) {
                    const Scalar& po = (*this)[pressureSwitchIdx];
                    setPrimaryVarsMeaningGas(Rs);
                    Scalar soMax = problem.maxOilSaturation(globalDofIdx);
                    Scalar rsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
                    Scalar rsSat = enableExtbo ? ExtboModule::rs(pvtRegionIndex(),
                                                                 po,
                                                                 zFraction_())
                                 : FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                                       T,
                                                                                       po,
                                                                                       s,
                                                                                       soMax);
                    (*this)[Indices::compositionSwitchIdx] = std::min(rsMax, rsSat);
                    changed = true;
                }
                Scalar so = 1.0 - sw - solventSaturation_() - sg;
                if (so < -eps && sg > 0.0 && FluidSystem::enableVaporizedOil()) {
                    // the oil phase disappeared and some hydrocarbon gas phase is still
                    // present, i.e., switch the primary variables to Rv.
                    // we only have the oil pressure readily available, but we need the gas
                    // pressure, i.e. we must determine capillary pressure
                    const Scalar& po = (*this)[pressureSwitchIdx];
                    std::array<Scalar, numPhases> pC = { 0.0 };
                    const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                    computeCapillaryPressures_(pC, /*so=*/0.0, sg + solventSaturation_(), sw, matParams);
                    Scalar pg = po + (pC[gasPhaseIdx] - pC[oilPhaseIdx]);

                    // we start at the Rv value that corresponds to that of oil-saturated
                    // hydrocarbon gas
                    setPrimaryVarsMeaningPressure(Pg);
                    (*this)[Indices::pressureSwitchIdx] = pg;
                    Scalar soMax = problem.maxOilSaturation(globalDofIdx);
                    Scalar rvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
                    Scalar rvSat = enableExtbo ? ExtboModule::rv(pvtRegionIndex(),
                                                                 pg,
                                                                 zFraction_())
                                 : FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                                        T,
                                                                                        pg,
                                                                                        Scalar(0),
                                                                                        soMax);
                    setPrimaryVarsMeaningGas(Rv);
                    (*this)[Indices::compositionSwitchIdx] = std::min(rvMax, rvSat);
                    changed = true;
                }
                break;
            }
            case Rs:
            {
                // Gas phase not present. The hydrocarbon gas phase
                // appears as soon as more of the gas component is present in the oil phase
                // than what saturated oil can hold.
                const Scalar& po = (*this)[pressureSwitchIdx];
                Scalar so = 1.0 - sw - solventSaturation_();
                Scalar soMax = std::max(so, problem.maxOilSaturation(globalDofIdx));
                Scalar rsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar rsSat = enableExtbo ? ExtboModule::rs(pvtRegionIndex(),
                                                         po,
                                                         zFraction_())
                         : FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                               T,
                                                                               po,
                                                                               so,
                                                                               soMax);

                Scalar rs = (*this)[Indices::compositionSwitchIdx];
                if (rs > std::min(rsMax, rsSat*(1.0 + eps))) {
                    // the gas phase appears, i.e., switch the primary variables to Sg
                    setPrimaryVarsMeaningGas(Sg);
                    (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation
                    changed = true;
                }
                break;
            }
            case Rv:
            {
                // The oil phase appears as
                // soon as more of the oil component is present in the hydrocarbon gas phase
                // than what saturated gas contains. Note that we use the blackoil specific
                // low-level PVT objects here for performance reasons.
                const Scalar& pg = (*this)[pressureSwitchIdx];
                Scalar soMax = problem.maxOilSaturation(globalDofIdx);
                Scalar rvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar rvSat = enableExtbo ? ExtboModule::rv(pvtRegionIndex(),
                                                            pg,
                                                            zFraction_())
                            : FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                                    T,
                                                                                    pg,
                                                                                    /*so=*/Scalar(0.0),
                                                                                    soMax);

                Scalar rv = (*this)[Indices::compositionSwitchIdx];
                if (rv > std::min(rvMax, rvSat*(1.0 + eps))) {
                    // switch to phase equilibrium mode because the oil phase appears. here
                    // we also need the capillary pressures to calculate the oil phase
                    // pressure using the gas phase pressure
                    Scalar sg2 = 1.0 - sw - solventSaturation_();
                    std::array<Scalar, numPhases> pC = { 0.0 };
                    const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                    computeCapillaryPressures_(pC,
                                            /*so=*/0.0,
                                            /*sg=*/sg2 + solventSaturation_(),
                                            sw,
                                            matParams);
                    Scalar po = pg + (pC[oilPhaseIdx] - pC[gasPhaseIdx]);

                    setPrimaryVarsMeaningGas(Sg);
                    setPrimaryVarsMeaningPressure(Po);
                    (*this)[Indices::pressureSwitchIdx] = po;
                    (*this)[Indices::compositionSwitchIdx] = sg2; // hydrocarbon gas saturation
                    changed = true;
                }
                break;
            }
            //case Rsw:
            //{
                //TODO
            //    break;
            //}
            case G_disabled:
            {
                break;
            }
            default:
                throw std::logic_error("No valid primary variable selected for water");
        }
        return changed;
    }

    bool chopAndNormalizeSaturations(){
        if (primaryVarsMeaningWater() == W_disabled && primaryVarsMeaningGas() == G_disabled){
            return false;
        }
        Scalar sw = 0.0;
        if (primaryVarsMeaningWater() == Sw)
            sw = (*this)[Indices::waterSaturationIdx];
        Scalar sg = 0.0;
        if (primaryVarsMeaningGas() == Sg)
            sg = (*this)[Indices::compositionSwitchIdx];

        Scalar ssol = 0.0;
        if constexpr (enableSolvent)
            ssol =(*this) [Indices::solventSaturationIdx];

        Scalar so = 1.0 - sw - sg - ssol;
        sw = std::min(std::max(sw,0.0),1.0);
        so = std::min(std::max(so,0.0),1.0);
        sg = std::min(std::max(sg,0.0),1.0);
        ssol = std::min(std::max(ssol,0.0),1.0);
        Scalar st = sw + so + sg + ssol;
        sw = sw/st;
        sg = sg/st;
        ssol = ssol/st;
        assert(st>0.5);
        if (primaryVarsMeaningWater() == Sw)
            (*this)[Indices::waterSaturationIdx] = sw;
        if (primaryVarsMeaningGas() == Sg)
            (*this)[Indices::compositionSwitchIdx] = sg;
        if constexpr (enableSolvent)
            (*this) [Indices::solventSaturationIdx] = ssol;

        return !(st==1);
    }

    BlackOilPrimaryVariables& operator=(const BlackOilPrimaryVariables& other) = default;
    BlackOilPrimaryVariables& operator=(Scalar value)
    {
        for (unsigned i = 0; i < numEq; ++i)
            (*this)[i] = value;

        return *this;
    }

    /*!
     * \brief Instruct valgrind to check the definedness of all attributes of this class.
     *
     * We cannot simply check the definedness of the whole object because there might be
     * "alignedness holes" in the memory layout which are caused by the pseudo primary
     * variables.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        // check the "real" primary variables
        for (unsigned i = 0; i < this->size(); ++i)
            Valgrind::CheckDefined((*this)[i]);

        // check the "pseudo" primary variables
        Valgrind::CheckDefined(primaryVarsMeaningWater_);
        Valgrind::CheckDefined(primaryVarsMeaningGas_);
        Valgrind::CheckDefined(primaryVarsMeaningPressure_);
        Valgrind::CheckDefined(primaryVarsMeaningBrine_);

        Valgrind::CheckDefined(pvtRegionIdx_);
#endif // NDEBUG
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar solventSaturation_() const
    {
        if constexpr (enableSolvent)
            return (*this)[Indices::solventSaturationIdx];
        else
            return 0.0;
    }

    Scalar zFraction_() const
    {
        if constexpr (enableExtbo)
            return (*this)[Indices::zFractionIdx];
        else
            return 0.0;
    }

    Scalar polymerConcentration_() const
    {
        if constexpr (enablePolymer)
            return (*this)[Indices::polymerConcentrationIdx];
        else
            return 0.0;
    }

    Scalar foamConcentration_() const
    {
        if constexpr (enableFoam)
            return (*this)[Indices::foamConcentrationIdx];
        else
            return 0.0;
    }

    Scalar saltConcentration_() const
    {
        if constexpr (enableBrine)
            return (*this)[Indices::saltConcentrationIdx];
        else
            return 0.0;
    }

    Scalar temperature_() const
    {
        if constexpr (enableEnergy)
            return (*this)[Indices::temperatureIdx];
        else
            return FluidSystem::reservoirTemperature();
    }

    Scalar microbialConcentration_() const
    {
        if constexpr (enableMICP)
            return (*this)[Indices::microbialConcentrationIdx];
        else
            return 0.0;
    }

    Scalar oxygenConcentration_() const
    {
        if constexpr (enableMICP)
            return (*this)[Indices::oxygenConcentrationIdx];
        else
            return 0.0;
    }

    Scalar ureaConcentration_() const
    {
        if constexpr (enableMICP)
            return (*this)[Indices::ureaConcentrationIdx];
        else
            return 0.0;
    }

    Scalar biofilmConcentration_() const
    {
        if constexpr (enableMICP)
            return (*this)[Indices::biofilmConcentrationIdx];
        else
            return 0.0;
    }

    Scalar calciteConcentration_() const
    {
        if constexpr (enableMICP)
            return (*this)[Indices::calciteConcentrationIdx];
        else
            return 0.0;
    }

    template <class Container>
    void computeCapillaryPressures_(Container& result,
                                    Scalar so,
                                    Scalar sg,
                                    Scalar sw,
                                    const MaterialLawParams& matParams) const
    {
        using SatOnlyFluidState = SimpleModularFluidState<Scalar,
                                                          numPhases,
                                                          numComponents,
                                                          FluidSystem,
                                                          /*storePressure=*/false,
                                                          /*storeTemperature=*/false,
                                                          /*storeComposition=*/false,
                                                          /*storeFugacity=*/false,
                                                          /*storeSaturation=*/true,
                                                          /*storeDensity=*/false,
                                                          /*storeViscosity=*/false,
                                                          /*storeEnthalpy=*/false>;
        SatOnlyFluidState fluidState;
        fluidState.setSaturation(waterPhaseIdx, sw);
        fluidState.setSaturation(oilPhaseIdx, so);
        fluidState.setSaturation(gasPhaseIdx, sg);

        MaterialLaw::capillaryPressures(result, matParams, fluidState);
    }

    PrimaryVarsMeaningWater primaryVarsMeaningWater_;
    PrimaryVarsMeaningPressure primaryVarsMeaningPressure_;
    PrimaryVarsMeaningGas primaryVarsMeaningGas_;
    PrimaryVarsMeaningBrine primaryVarsMeaningBrine_;
    unsigned short pvtRegionIdx_;
};


} // namespace Opm

#endif
