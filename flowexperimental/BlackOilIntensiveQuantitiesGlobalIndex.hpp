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
 * \copydoc Opm::BlackOilIntensiveQuantitiesGlobalIndex
 */
#ifndef OPM_BLACK_OIL_INTENSIVE_QUANTITIES_GLOBAL_INDEX_HPP
#define OPM_BLACK_OIL_INTENSIVE_QUANTITIES_GLOBAL_INDEX_HPP

#include "BlackOilEnergyIntensiveQuantitiesGlobalIndex.hpp"

#include <dune/common/fmatrix.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/blackoil/blackoilextbomodules.hh>
#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoildiffusionmodule.hh>
#include <opm/models/blackoil/blackoilbioeffectsmodules.hh>
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/common/directionalmobility.hh>

#include <opm/utility/CopyablePtr.hpp>

#include <stdexcept>
#include <utility>

#include <fmt/format.h>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model using global indices
 */
template <class TypeTag>
class BlackOilIntensiveQuantitiesGlobalIndex
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
    , public BlackOilDiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>()>
    , public BlackOilSolventIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableSolvent>()>
    , public BlackOilExtboIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableExtbo>()>
    , public BlackOilPolymerIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnablePolymer>()>
    , public BlackOilFoamIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableFoam>()>
    , public BlackOilBrineIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableBrine>()>
    , public BlackOilEnergyIntensiveQuantitiesGlobalIndex<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>
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
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableConvectiveMixing = getPropValue<TypeTag, Properties::EnableConvectiveMixing>() };
    enum { enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>() };
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

    using DirectionalMobilityPtr = Utility::CopyablePtr<DirectionalMobility<TypeTag>>;

public:
    using FluidState = BlackOilFluidState<Evaluation,
                                          FluidSystem,
                                          enableTemperature,
                                          enableEnergy,
                                          compositionSwitchEnabled,
                                          enableVapwat,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          false,
                                          Indices::numPhases>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    BlackOilIntensiveQuantitiesGlobalIndex()
    {
        if constexpr (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
        }
        if constexpr (enableVapwat) {
            fluidState_.setRvw(0.0);
        }
    }

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        this->update(problem,priVars,globalSpaceIdx,timeIdx);
    }

    void update(const Problem& problem,
                const PrimaryVariables& priVars,
                unsigned globalSpaceIdx,
                unsigned timeIdx)
    {

        this->extrusionFactor_ = 1.0;// to avoid fixing parent update
        OPM_TIMEBLOCK_LOCAL(UpdateIntensiveQuantities, Subsystem::PvtProps | Subsystem::SatProps);
        Scalar RvMax = FluidSystem::enableVaporizedOil()
            ? problem.maxOilVaporizationFactor(timeIdx, globalSpaceIdx)
            : 0.0;
        Scalar RsMax = FluidSystem::enableDissolvedGas()
            ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
            : 0.0;

        asImp_().updateTemperature_(problem, priVars, globalSpaceIdx, timeIdx);

        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        //asImp_().updateSaltConcentration_(elemCtx, dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                Sw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Disabled) {
                // water is enabled but is not a primary variable i.e. one phase case
                Sw = 1.0;
            }
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
            So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
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

        std::array<Evaluation, numPhases> pC{};
        computeRelpermAndPC(mobility_, pC, problem, fluidState_, globalSpaceIdx);
        // oil is the reference phase for pressure
        if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
                }
            }
        } else {
            const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
                }
            }
        }

        Evaluation SoMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            SoMax = max(fluidState_.saturation(oilPhaseIdx),
                        problem.maxOilSaturation(globalSpaceIdx));
        }

        // take the meaning of the switching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
            const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRs(Rs);
        } else {
            if (FluidSystem::enableDissolvedGas()) { // Add So > 0? i.e. if only water set rs = 0)
                OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRs, Subsystem::PvtProps);
                const Evaluation& RsSat = enableExtbo ? asImp_().rs() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRs(min(RsMax, RsSat));
            }
            else if constexpr (compositionSwitchEnabled) {
                fluidState_.setRs(0.0);
            }
        }

        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
            const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRv(Rv);
        } else {
            if (FluidSystem::enableVaporizedOil() ) { // Add Sg > 0? i.e. if only water set rv = 0)
                OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRv, Subsystem::PvtProps);
                //NB! should save the indexing for later evalustion
                const Evaluation& RvSat = enableExtbo ? asImp_().rv() :
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRv(min(RvMax, RvSat));
            }
            else if constexpr (compositionSwitchEnabled) {
                fluidState_.setRv(0.0);
            }
        }

        if constexpr (enableVapwat) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rvw) {
                const auto& Rvw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
                fluidState_.setRvw(Rvw);
            } else {
                //NB! should save the indexing for later evaluation
                if (FluidSystem::enableVaporizedWater()) { // Add Sg > 0? i.e. if only water set rv = 0)
                    OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRv, Subsystem::PvtProps);
                    const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                                                        gasPhaseIdx,
                                                                                        pvtRegionIdx);
                    fluidState_.setRvw(RvwSat);
                }
            }
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            computeInverseFormationVolumeFactorAndViscosity(fluidState_, phaseIdx, pvtRegionIdx);
        }
        Valgrind::CheckDefined(mobility_);

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            OPM_TIMEBLOCK_LOCAL(UpdateWDensity, Subsystem::PvtProps);
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            OPM_TIMEBLOCK_LOCAL(UpdateGDensity, Subsystem::PvtProps);
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
            OPM_TIMEBLOCK_LOCAL(UpdateODensity, Subsystem::PvtProps);
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {
                rho += fluidState_.invB(oilPhaseIdx) *
                       fluidState_.Rs() *
                       FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }

        // retrieve the porosity from the problem
        referencePorosity_ = problem.porosity(globalSpaceIdx,timeIdx);
        porosity_ = referencePorosity_;
        // the porosity must be modified by the compressibility of the
        // rock...
        Scalar rockCompressibility = problem.rockCompressibility(globalSpaceIdx);
        if (rockCompressibility > 0.0) {
            OPM_TIMEBLOCK_LOCAL(UpdateRockCompressibility, Subsystem::PvtProps);
            Scalar rockRefPressure = problem.rockReferencePressure(globalSpaceIdx);
            Evaluation x;
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                x = rockCompressibility*(fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            } else if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                x = rockCompressibility*(fluidState_.pressure(waterPhaseIdx) - rockRefPressure);
            } else {
                x = rockCompressibility*(fluidState_.pressure(gasPhaseIdx) - rockRefPressure);
            }
            porosity_ *= 1.0 + x + 0.5 * x * x;
        }

        // deal with water induced rock compaction
        porosity_ *= problem.template rockCompPoroMultiplier<Evaluation>(*this, globalSpaceIdx);

        rockCompTransMultiplier_ = problem.template rockCompTransMultiplier<Evaluation>(*this, globalSpaceIdx);

        if constexpr (enableConvectiveMixing) {
            asImp_().updateSaturatedDissolutionFactor_();
        }

#ifndef NDEBUG
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
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
            case Dir::XPlus:
                return dirMob_->getArray(0)[phaseIdx];
            case Dir::YPlus:
                return dirMob_->getArray(1)[phaseIdx];
            case Dir::ZPlus:
                return dirMob_->getArray(2)[phaseIdx];
            default:
                throw std::runtime_error("Unexpected face direction");
            }
        } else {
            return mobility_[phaseIdx];
        }
    }

    void computeInverseFormationVolumeFactorAndViscosity(FluidState& fluidState,
                                                         unsigned phaseIdx,
                                                         unsigned pvtRegionIdx) {
        OPM_TIMEBLOCK_LOCAL(UpdateInverseFormationFactorAndViscosity, Subsystem::PvtProps);
        {
            OPM_TIMEBLOCK_LOCAL(UpdateFormationFactor, Subsystem::PvtProps);
            const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
            fluidState_.setInvB(phaseIdx, b);
        }
        {
            OPM_TIMEBLOCK_LOCAL(UpdateViscosity, Subsystem::PvtProps);
            typename FluidSystem::template ParameterCache<Evaluation> paramCache;
            paramCache.setRegionIndex(pvtRegionIdx);
            paramCache.updateAll(fluidState_);

            const auto& mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            mobility_[phaseIdx] /= mu;
        }
    }

    void computeRelpermAndPC(std::array<Evaluation,numPhases>& mobility,
                             std::array<Evaluation, numPhases>& pC,
                             const Problem& problem,
                             const FluidState& fluidState,
                             unsigned globalSpaceIdx){
        OPM_TIMEBLOCK_LOCAL(UpdateRelperm, Subsystem::SatProps);
        const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        problem.updateRelperms(mobility, dirMob_, fluidState, globalSpaceIdx);
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
    auto pvtRegionIndex() const -> decltype(std::declval<FluidState>().pvtRegionIndex())
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

    const Evaluation& permFactor() const
    {
        throw std::logic_error("permFactor() is not yet implemented for compositional modeling"); 
    }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag, enableSolvent>;
    friend BlackOilExtboIntensiveQuantities<TypeTag, enableExtbo>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag, enablePolymer>;
    friend BlackOilEnergyIntensiveQuantitiesGlobalIndex<TypeTag, enableEnergy>;
    friend BlackOilFoamIntensiveQuantities<TypeTag, enableFoam>;
    friend BlackOilBrineIntensiveQuantities<TypeTag, enableBrine>;
    friend BlackOilBioeffectsIntensiveQuantities<TypeTag, enableBioeffects>;

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
