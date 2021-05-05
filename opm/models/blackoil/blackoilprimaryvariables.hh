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

    static const bool compositionSwitchEnabled = Indices::gasEnabled;
    static const bool waterEnabled = Indices::waterEnabled;

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
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
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

    static_assert(numPhases == 3, "The black-oil model assumes three phases!");
    static_assert(numComponents == 3, "The black-oil model assumes three components!");

public:
    enum PrimaryVarsMeaning {
        Sw_po_Sg, // threephase case
        Sw_po_Rs, // water + oil case
        Sw_pg_Rv, // water + gas case
        OnePhase_p, // onephase case
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
        Valgrind::SetUndefined(primaryVarsMeaning_);
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
    PrimaryVarsMeaning primaryVarsMeaning() const
    { return primaryVarsMeaning_; }

    /*!
     * \brief Set the interpretation which should be applied to the switching primary
     *        variables.
     */
    void setPrimaryVarsMeaning(PrimaryVarsMeaning newMeaning)
    { primaryVarsMeaning_ = newMeaning; }

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
        static const Scalar thresholdWaterFilledCell = 1.0 - 1e-6;
        bool onlyWater = FluidSystem::phaseIsActive(waterPhaseIdx)?(fluidState.saturation(waterPhaseIdx) > thresholdWaterFilledCell):false;

        // deal with the primary variables for the energy extension
        EnergyModule::assignPrimaryVars(*this, fluidState);

        // determine the meaning of the primary variables
        if (FluidSystem::numActivePhases() == 1) {
            primaryVarsMeaning_ = OnePhase_p;
        }
        else if ((gasPresent && oilPresent) || (onlyWater && FluidSystem::phaseIsActive(oilPhaseIdx))) {
            // gas and oil: both hydrocarbon phases are in equilibrium (i.e., saturated
            // with the "protagonist" component of the other phase.)
            primaryVarsMeaning_ = Sw_po_Sg;
        }
        else if (oilPresent) {
            // only oil: if dissolved gas is enabled, we need to consider the oil phase
            // composition, if it is disabled, the gas component must stick to its phase
            if (FluidSystem::enableDissolvedGas())
                primaryVarsMeaning_ = Sw_po_Rs;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }
        else {
            assert(gasPresent);
            // only gas: if vaporized oil is enabled, we need to consider the gas phase
            // composition, if it is disabled, the oil component must stick to its phase
            if (FluidSystem::enableVaporizedOil())
                primaryVarsMeaning_ = Sw_pg_Rv;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }

        // assign the actual primary variables
        if (primaryVarsMeaning() == OnePhase_p) {
            if (waterEnabled) {
                (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(waterPhaseIdx));
            } else {
                throw std::logic_error("For single-phase runs, only pure water is presently allowed.");
            }
            
        }
        else if (primaryVarsMeaning() == Sw_po_Sg) {
            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = FsToolbox::value(fluidState.saturation(gasPhaseIdx));
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            const auto& Rs = BlackOil::getRs_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);

            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = Rs;
        }
        else {
            assert(primaryVarsMeaning() == Sw_pg_Rv);

            const auto& Rv = BlackOil::getRv_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));

            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(gasPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = Rv;
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
        if (primaryVarsMeaning() == OnePhase_p){
            return false;
        }
        Scalar Sw = 0.0;
        if (waterEnabled)
            Sw = (*this)[Indices::waterSaturationIdx];

        if (primaryVarsMeaning() == Sw_po_Sg) {

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {

                // make sure water saturations does not exceed 1.0
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0;
                // the hydrocarbon gas saturation is set to 0.0
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] = 0.0;

                return false;
            }

            // phase equilibrium, i.e., both hydrocarbon phases are present.
            Scalar Sg = 0.0;
            if (compositionSwitchEnabled)
                Sg = (*this)[Indices::compositionSwitchIdx];

            Scalar So = 1.0 - Sw - Sg - solventSaturation_();

            Scalar So2 = 1.0 - Sw - solventSaturation_();
            if (Sg < -eps && So2 > 0.0 && FluidSystem::enableDissolvedGas()) {
                // the hydrocarbon gas phase disappeared and some oil phase is left,
                // i.e., switch the primary variables to { Sw, po, Rs }.
                //
                // by a "lucky" coincidence the pressure switching variable already
                // represents the oil phase pressure, so we do not need to change
                // this. For the gas dissolution factor, we use the low-level blackoil
                // PVT objects to calculate the mole fraction of gas saturated oil.
                Scalar po = (*this)[Indices::pressureSwitchIdx];
                Scalar T = asImp_().temperature_();
                Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
                Scalar RsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar RsSat = enableExtbo ? ExtboModule::rs(pvtRegionIndex(),
                                                             po,
                                                             zFraction_())
                             : FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                                   T,
                                                                                   po,
                                                                                   So2,
                                                                                   SoMax);

                setPrimaryVarsMeaning(Sw_po_Rs);
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] =
                        std::min(RsMax, RsSat);

                return true;
            }

            Scalar Sg2 = 1.0 - Sw - solventSaturation_();
            if (So < -eps && Sg2 > 0.0 && FluidSystem::enableVaporizedOil()) {
                // the oil phase disappeared and some hydrocarbon gas phase is still
                // present, i.e., switch the primary variables to { Sw, pg, Rv }.
                Scalar po = (*this)[Indices::pressureSwitchIdx];

                // we only have the oil pressure readily available, but we need the gas
                // pressure, i.e. we must determine capillary pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC, /*So=*/0.0, Sg2 + solventSaturation_(), Sw, matParams);
                Scalar pg = po + (pC[gasPhaseIdx] - pC[oilPhaseIdx]);

                // we start at the Rv value that corresponds to that of oil-saturated
                // hydrocarbon gas
                Scalar T = asImp_().temperature_();
                Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
                Scalar RvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar RvSat = enableExtbo ? ExtboModule::rv(pvtRegionIndex(),
                                                             pg,
                                                             zFraction_())
                             : FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                                    T,
                                                                                    pg,
                                                                                    Scalar(0),
                                                                                    SoMax);
                setPrimaryVarsMeaning(Sw_pg_Rv);
                (*this)[Indices::pressureSwitchIdx] = pg;
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] = std::min(RvMax, RvSat);

                return true;
            }

            return false;
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            assert(compositionSwitchEnabled);

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {
                // switch back to phase equilibrium mode if the oil phase vanishes (i.e.,
                // the water-only case)
                setPrimaryVarsMeaning(Sw_po_Sg);
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0; // water saturation

                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            // Only the oil and the water phases are present. The hydrocarbon gas phase
            // appears as soon as more of the gas component is present in the oil phase
            // than what saturated oil can hold.
            Scalar T = asImp_().temperature_();
            Scalar po = (*this)[Indices::pressureSwitchIdx];
            Scalar So = 1.0 - Sw - solventSaturation_();
            Scalar SoMax = std::max(So, problem.maxOilSaturation(globalDofIdx));
            Scalar RsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
            Scalar RsSat = enableExtbo ? ExtboModule::rs(pvtRegionIndex(),
                                                         po,
                                                         zFraction_())
                         : FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                               T,
                                                                               po,
                                                                               So,
                                                                               SoMax);

            Scalar Rs = (*this)[Indices::compositionSwitchIdx];
            if (Rs > std::min(RsMax, RsSat*(1.0 + eps))) {
                // the gas phase appears, i.e., switch the primary variables to { Sw, po,
                // Sg }.
                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            return false;
        }
        else {
            assert(primaryVarsMeaning() == Sw_pg_Rv);
            assert(compositionSwitchEnabled);

            Scalar pg = (*this)[Indices::pressureSwitchIdx];
            Scalar Sg = 1.0 - Sw - solventSaturation_();

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {
                // switch to phase equilibrium mode because the hydrocarbon gas phase
                // disappears. here we need the capillary pressures to calculate the oil
                // phase pressure using the gas phase pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC,
                                           /*So=*/0.0,
                                           /*Sg=*/Sg + solventSaturation_(),
                                           Sw,
                                           matParams);
                Scalar po = pg + (pC[oilPhaseIdx] - pC[gasPhaseIdx]);

                setPrimaryVarsMeaning(Sw_po_Sg);
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0;

                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            // Only the gas and the water phases are present. The oil phase appears as
            // soon as more of the oil component is present in the hydrocarbon gas phase
            // than what saturated gas contains. Note that we use the blackoil specific
            // low-level PVT objects here for performance reasons.
            Scalar T = asImp_().temperature_();
            Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
            Scalar RvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
            Scalar RvSat = enableExtbo ? ExtboModule::rv(pvtRegionIndex(),
                                                         pg,
                                                         zFraction_())
                         : FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                                T,
                                                                                pg,
                                                                                /*So=*/Scalar(0.0),
                                                                                SoMax);

            Scalar Rv = (*this)[Indices::compositionSwitchIdx];
            if (Rv > std::min(RvMax, RvSat*(1.0 + eps))) {
                // switch to phase equilibrium mode because the oil phase appears. here
                // we also need the capillary pressures to calculate the oil phase
                // pressure using the gas phase pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC,
                                           /*So=*/0.0,
                                           /*Sg=*/Sg + solventSaturation_(),
                                           Sw,
                                           matParams);
                Scalar po = pg + (pC[oilPhaseIdx] - pC[gasPhaseIdx]);

                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = Sg; // hydrocarbon gas saturation

                return true;
            }

            return false;
        }

        assert(false);
        return false;
    }

    bool chopAndNormalizeSaturations(){
        if (primaryVarsMeaning() == OnePhase_p){
            return false;
        }
        Scalar Sw = 0.0;
        if (waterEnabled)
            Sw = (*this)[Indices::waterSaturationIdx];

        if (primaryVarsMeaning() == Sw_po_Sg) {
            Scalar Sg = 0.0;
            if (compositionSwitchEnabled)
                Sg = (*this)[Indices::compositionSwitchIdx];
            Scalar Ssol = 0.0;
            if (enableSolvent)
                Ssol =(*this) [Indices::solventSaturationIdx];
            
            Scalar So = 1.0 - Sw - Sg - Ssol;
            Sw = std::min(std::max(Sw,0.0),1.0);
            So = std::min(std::max(So,0.0),1.0);
            Sg = std::min(std::max(Sg,0.0),1.0);
            Ssol = std::min(std::max(Ssol,0.0),1.0);
            
            Scalar St = Sw + So + Sg + Ssol;
            Sw = Sw/St;
            Sg = Sg/St;
            Ssol = Ssol/St;
            assert(St>0.5);
            if (waterEnabled)
                (*this)[Indices::waterSaturationIdx]= Sw;
            if (compositionSwitchEnabled)
                (*this)[Indices::compositionSwitchIdx] = Sg;
            if (enableSolvent)
                (*this)[Indices::solventSaturationIdx] = Ssol;
            return not(St==1);
            
        }else if (primaryVarsMeaning() == Sw_po_Rs) { 
            Scalar Ssol = 0.0;
            if (enableSolvent)
                Ssol = (*this)[Indices::solventSaturationIdx];
            Scalar So = 1.0 - Sw  - Ssol;
            Sw = std::min(std::max(Sw,0.0),1.0);
            So = std::min(std::max(So,0.0),1.0);
            Ssol = std::min(std::max(Ssol,0.0),1.0);
            //Sg = 0.0;
            Scalar St = Sw + So + Ssol;
            assert(St>0.5);
            Sw=Sw/St;
            Ssol=Ssol/St;
            if (waterEnabled)
                (*this)[Indices::waterSaturationIdx]= Sw;
            if (enableSolvent)
                (*this)[Indices::solventSaturationIdx] = Ssol;            
            return not(St==1);
        }else{
            assert(primaryVarsMeaning() == Sw_pg_Rv);
            assert(compositionSwitchEnabled);
            Scalar Ssol=0.0;
            if (enableSolvent)
                Ssol = (*this)[Indices::solventSaturationIdx];
            Scalar Sg = 1.0 - Sw - Ssol;
            //So = 0.0;
            Sw = std::min(std::max(Sw,0.0),1.0);
            Sg = std::min(std::max(Sg,0.0),1.0);
            Ssol = std::min(std::max(Ssol,0.0),1.0);
            Scalar St = Sw + Sg + Ssol;
            assert(St>0.5);
            Sw=Sw/St;
            Ssol=Ssol/St;
            if (waterEnabled)
                (*this)[Indices::waterSaturationIdx]= Sw;
            if (enableSolvent)
                (*this)[Indices::solventSaturationIdx] = Ssol;
            return not(St==1);            
        }
        assert(false); 
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
        Valgrind::CheckDefined(primaryVarsMeaning_);
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
        if (!enableSolvent)
            return 0.0;

        return (*this)[Indices::solventSaturationIdx];
    }

    Scalar zFraction_() const
    {
        if (!enableExtbo)
            return 0.0;

        return (*this)[Indices::zFractionIdx];
    }

    Scalar polymerConcentration_() const
    {
        if (!enablePolymer)
            return 0.0;

        return (*this)[Indices::polymerConcentrationIdx];
    }

    Scalar foamConcentration_() const
    {
        if (!enableFoam)
            return 0.0;

        return (*this)[Indices::foamConcentrationIdx];
    }

    Scalar saltConcentration_() const
    {
        if (!enableBrine)
            return 0.0;

        return (*this)[Indices::saltConcentrationIdx];
    }

    Scalar temperature_() const
    {
        if (!enableEnergy)
            return FluidSystem::reservoirTemperature();

        return (*this)[Indices::temperatureIdx];
    }

    template <class Container>
    void computeCapillaryPressures_(Container& result,
                                    Scalar So,
                                    Scalar Sg,
                                    Scalar Sw,
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
        fluidState.setSaturation(waterPhaseIdx, Sw);
        fluidState.setSaturation(oilPhaseIdx, So);
        fluidState.setSaturation(gasPhaseIdx, Sg);

        MaterialLaw::capillaryPressures(result, matParams, fluidState);
    }

    PrimaryVarsMeaning primaryVarsMeaning_;
    unsigned short pvtRegionIdx_;
};


} // namespace Opm

#endif
