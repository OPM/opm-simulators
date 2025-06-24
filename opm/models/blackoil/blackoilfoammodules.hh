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
 * \brief Contains the classes required to extend the black-oil model to include the effects of foam.
 */
#ifndef EWOMS_BLACK_OIL_FOAM_MODULE_HH
#define EWOMS_BLACK_OIL_FOAM_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/models/blackoil/blackoilfoamparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <cassert>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace Opm {

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model to include the effects of foam.
 */
template <class TypeTag, bool enableFoamV = getPropValue<TypeTag, Properties::EnableFoam>()>
class BlackOilFoamModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilFoamParams<Scalar>::TabulatedFunction;

    static constexpr unsigned foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned contiFoamEqIdx = Indices::contiFoamEqIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableFoam = enableFoamV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;

    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };

public:
    //! \brief Set parameters.
    static void setParams(BlackOilFoamParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil foam module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Register all foam specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model&,
                                      Simulator&)
    {
        if constexpr (enableFoam) {
            if (Parameters::Get<Parameters::EnableVtkOutput>()) {
                OpmLog::warning("VTK output requested, currently unsupported by the foam module.");
            }
        }
        //model.addOutputModule(new VtkBlackOilFoamModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if constexpr (enableFoam) {
            return pvIdx == foamConcentrationIdx;
        }
        else {
            return false;
        }
    }

    static std::string primaryVarName([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));
        return "foam_concentration";
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
       assert(primaryVarApplies(pvIdx));

       // TODO: it may be beneficial to chose this differently.
       return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableFoam) {
            return eqIdx == contiFoamEqIdx;
        }
        else {
            return false;
        }
    }

    static std::string eqName([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^foam";
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
       assert(eqApplies(eqIdx));

       // TODO: it may be beneficial to chose this differently.
       return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if constexpr (enableFoam) {
            const auto& fs = intQuants.fluidState();

            LhsEval surfaceVolume = Toolbox::template decay<LhsEval>(intQuants.porosity());
            if (params_.transport_phase_ == Phase::WATER) {
                surfaceVolume *= Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx)) *
                                 Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx));
            } else if (params_.transport_phase_ == Phase::GAS) {
                surfaceVolume *= Toolbox::template decay<LhsEval>(fs.saturation(gasPhaseIdx)) *
                                 Toolbox::template decay<LhsEval>(fs.invB(gasPhaseIdx));
            } else if (params_.transport_phase_ == Phase::SOLVENT) {
                if constexpr (enableSolvent) {
                    surfaceVolume *= Toolbox::template decay<LhsEval>( intQuants.solventSaturation()) *
                                     Toolbox::template decay<LhsEval>(intQuants.solventInverseFormationVolumeFactor());
                }
            } else {
                throw std::runtime_error("Transport phase is GAS/WATER/SOLVENT");
            }

            // Avoid singular matrix if no gas is present.
            surfaceVolume = max(surfaceVolume, 1e-10);

            // Foam/surfactant in free phase.
            const LhsEval freeFoam = surfaceVolume *
                                     Toolbox::template decay<LhsEval>(intQuants.foamConcentration());

            // Adsorbed foam/surfactant.
            const LhsEval adsorbedFoam =
                Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity()) *
                Toolbox::template decay<LhsEval>(intQuants.foamRockDensity()) *
                Toolbox::template decay<LhsEval>(intQuants.foamAdsorbed());

            const LhsEval accumulationFoam = freeFoam + adsorbedFoam;
            storage[contiFoamEqIdx] += accumulationFoam;
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableFoam) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
            const unsigned inIdx = extQuants.interiorIndex();

            // The effect of the mobility reduction factor is
            // incorporated in the mobility for the relevant phase,
            // so fluxes do not need modification here.
            switch (transportPhase()) {
                case Phase::WATER: {
                    const unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
                    const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
                    if (upIdx == inIdx) {
                        flux[contiFoamEqIdx] =
                            extQuants.volumeFlux(waterPhaseIdx) *
                            up.fluidState().invB(waterPhaseIdx) *
                            up.foamConcentration();
                    } else {
                        flux[contiFoamEqIdx] =
                            extQuants.volumeFlux(waterPhaseIdx) *
                            decay<Scalar>(up.fluidState().invB(waterPhaseIdx)) *
                            decay<Scalar>(up.foamConcentration());
                    }
                    break;
                }
                case Phase::GAS: {
                    const unsigned upIdx = extQuants.upstreamIndex(gasPhaseIdx);
                    const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
                    if (upIdx == inIdx) {
                        flux[contiFoamEqIdx] =
                            extQuants.volumeFlux(gasPhaseIdx) *
                            up.fluidState().invB(gasPhaseIdx) *
                            up.foamConcentration();
                    } else {
                        flux[contiFoamEqIdx] =
                            extQuants.volumeFlux(gasPhaseIdx) *
                            decay<Scalar>(up.fluidState().invB(gasPhaseIdx)) *
                            decay<Scalar>(up.foamConcentration());
                    }
                    break;
                }
                case Phase::SOLVENT:
                    if constexpr (enableSolvent) {
                        const unsigned upIdx = extQuants.solventUpstreamIndex();
                        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
                        if (upIdx == inIdx) {
                            flux[contiFoamEqIdx] =
                                extQuants.solventVolumeFlux() *
                                up.solventInverseFormationVolumeFactor() *
                                up.foamConcentration();
                        } else {
                            flux[contiFoamEqIdx] =
                                extQuants.solventVolumeFlux() *
                                decay<Scalar>(up.solventInverseFormationVolumeFactor()) *
                                decay<Scalar>(up.foamConcentration());
                        }
                    } else {
                        throw std::runtime_error("Foam transport phase is SOLVENT but SOLVENT is not activated.");
                    }
                    break;
                default:
                    throw std::runtime_error("Foam transport phase must be GAS/WATER/SOLVENT.");
            }
        }
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables&,
                                     const EqVector&)
    {
        // do not consider the change of foam primary variables for convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity([[maybe_unused]] const Model& model,
                                [[maybe_unused]] std::ostream& outstream,
                                [[maybe_unused]] const DofEntity& dof)
    {
        if constexpr (enableFoam) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
            outstream << priVars[foamConcentrationIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity([[maybe_unused]] Model& model,
                                  [[maybe_unused]] std::istream& instream,
                                  [[maybe_unused]] const DofEntity& dof)
    {
        if constexpr (enableFoam) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
            PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

            instream >> priVars0[foamConcentrationIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1[foamConcentrationIdx] = priVars0[foamConcentrationIdx];
        }
    }

    static const Scalar foamRockDensity(const ElementContext& elemCtx,
                                        unsigned scvIdx,
                                        unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.foamRockDensity_[satnumRegionIdx];
    }

    static bool foamAllowDesorption(const ElementContext& elemCtx,
                                    unsigned scvIdx,
                                    unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.foamAllowDesorption_[satnumRegionIdx];
    }

    static const TabulatedFunction& adsorbedFoamTable(const ElementContext& elemCtx,
                                                      unsigned scvIdx,
                                                      unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.adsorbedFoamTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& gasMobilityMultiplierTable(const ElementContext& elemCtx,
                                                               unsigned scvIdx,
                                                               unsigned timeIdx)
    {
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.gasMobilityMultiplierTable_[pvtnumRegionIdx];
    }

    static const typename BlackOilFoamParams<Scalar>::FoamCoefficients&
    foamCoefficients(const ElementContext& elemCtx,
                     const unsigned scvIdx,
                     const unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.foamCoefficients_[satnumRegionIdx];
    }

    static Phase transportPhase()
    { return params_.transport_phase_; }

private:
    static BlackOilFoamParams<Scalar> params_;
};

template <class TypeTag, bool enableFoam>
BlackOilFoamParams<typename BlackOilFoamModule<TypeTag, enableFoam>::Scalar>
BlackOilFoamModule<TypeTag, enableFoam>::params_;

template <class TypeTag, bool enableFoam>
class BlackOilFoamIntensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilFoamIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag>
class BlackOilFoamIntensiveQuantities<TypeTag, /*enableFoam=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using FoamModule = BlackOilFoamModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };

    static constexpr int foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void foamPropertiesUpdate_(const ElementContext& elemCtx,
                               unsigned dofIdx,
                               unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        foamConcentration_ = priVars.makeEvaluation(foamConcentrationIdx, timeIdx);
        const auto& fs = asImp_().fluidState_;

        // Compute gas mobility reduction factor
        Evaluation mobilityReductionFactor = 1.0;
        if constexpr (false) {
            // The functional model is used.
            // TODO: allow this model.
            // In order to do this we must allow transport to be in the water phase, not just the gas phase.
            const auto& foamCoefficients = FoamModule::foamCoefficients(elemCtx, dofIdx, timeIdx);

            const Scalar fm_mob = foamCoefficients.fm_mob;

            const Scalar fm_surf = foamCoefficients.fm_surf;
            const Scalar ep_surf = foamCoefficients.ep_surf;

            const Scalar fm_oil = foamCoefficients.fm_oil;
            const Scalar fl_oil = foamCoefficients.fl_oil;
            const Scalar ep_oil = foamCoefficients.ep_oil;

            const Scalar fm_dry = foamCoefficients.fm_dry;
            const Scalar ep_dry = foamCoefficients.ep_dry;

            const Scalar fm_cap = foamCoefficients.fm_cap;
            const Scalar ep_cap = foamCoefficients.ep_cap;

            const Evaluation C_surf = foamConcentration_;
            const Evaluation Ca = 1e10; // TODO: replace with proper capillary number.
            const Evaluation S_o = fs.saturation(oilPhaseIdx);
            const Evaluation S_w = fs.saturation(waterPhaseIdx);

            const Evaluation F1 = pow(C_surf / fm_surf, ep_surf);
            const Evaluation F2 = pow((fm_oil - S_o) / (fm_oil - fl_oil), ep_oil);
            const Evaluation F3 = pow(fm_cap / Ca, ep_cap);
            const Evaluation F7 = 0.5 + atan(ep_dry * (S_w - fm_dry)) / M_PI;

            mobilityReductionFactor = 1. / (1. + fm_mob * F1 * F2 * F3 * F7);
        } else {
            // The tabular model is used.
            // Note that the current implementation only includes the effect of foam concentration (FOAMMOB),
            // and not the optional pressure dependence (FOAMMOBP) or shear dependence (FOAMMOBS).
            const auto& gasMobilityMultiplier = FoamModule::gasMobilityMultiplierTable(elemCtx, dofIdx, timeIdx);
            mobilityReductionFactor = gasMobilityMultiplier.eval(foamConcentration_, /* extrapolate = */ true);
        }

        // adjust mobility
        switch (FoamModule::transportPhase()) {
            case Phase::WATER:
                asImp_().mobility_[waterPhaseIdx] *= mobilityReductionFactor;
                break;
            case Phase::GAS:
                asImp_().mobility_[gasPhaseIdx] *= mobilityReductionFactor;
                break;
            case Phase::SOLVENT:
                if constexpr (enableSolvent) {
                    asImp_().solventMobility_ *= mobilityReductionFactor;
                } else {
                    throw std::runtime_error("Foam transport phase is SOLVENT but SOLVENT is not activated.");
                }
                break;
            default:
                throw std::runtime_error("Foam transport phase must be GAS/WATER/SOLVENT.");
        }

        foamRockDensity_ = FoamModule::foamRockDensity(elemCtx, dofIdx, timeIdx);

        const auto& adsorbedFoamTable = FoamModule::adsorbedFoamTable(elemCtx, dofIdx, timeIdx);
        foamAdsorbed_ = adsorbedFoamTable.eval(foamConcentration_, /*extrapolate=*/true);
        if (!FoamModule::foamAllowDesorption(elemCtx, dofIdx, timeIdx)) {
            throw std::runtime_error("Foam module does not support the 'no desorption' option.");
        }
    }

    const Evaluation& foamConcentration() const
    { return foamConcentration_; }

    Scalar foamRockDensity() const
    { return foamRockDensity_; }

    const Evaluation& foamAdsorbed() const
    { return foamAdsorbed_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation foamConcentration_;
    Scalar foamRockDensity_;
    Evaluation foamAdsorbed_;
};

template <class TypeTag>
class BlackOilFoamIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void foamPropertiesUpdate_(const ElementContext&,
                               unsigned,
                               unsigned)
    {}

    const Evaluation& foamConcentration() const
    { throw std::runtime_error("foamConcentration() called but foam is disabled"); }

    Scalar foamRockDensity() const
    { throw std::runtime_error("foamRockDensity() called but foam is disabled"); }

    Scalar foamAdsorbed() const
    { throw std::runtime_error("foamAdsorbed() called but foam is disabled"); }
};

} // namespace Opm

#endif
