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
 * \brief Contains the classes required to extend the black-oil model by bioeffects.
 */
#ifndef OPM_BLACK_OIL_BIOEFFECTS_MODULE_HH
#define OPM_BLACK_OIL_BIOEFFECTS_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/models/blackoil/blackoilbioeffectsparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/io/vtkblackoilbioeffectsmodule.hpp>

#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by bioeffects.
 * 
 * The two implemented model extensions are MICP and biofilm effects in
 * underground storage. For details on the mathematical models, refer to the 
 * OPM Flow manual.
 * 
 * I) MICP (Microbially Induced Calcite Precipitation)
 * MICP is a novel and sustainable technology that leverages biochemical
 * processes to form barriers through calcium carbonate cementation. This
 * approach shows promise for sealing leakage zones in geological formations.
 * 
 * The conceptual model includes the following key mechanisms:
 * - Suspended microbes attach to pore walls, forming biofilm.
 * - A growth solution is introduced to stimulate biofilm development.
 * - The biofilm utilizes a cementation solution to produce calcite.
 * - Calcite precipitates reduce pore space, thereby decreasing rock
 *   permeability.
 * 
 * This implementation considers a single-phase (water) system with the
 * following primary variables:
 * - Pressure of the water phase
 * - Concentration of suspended microbes
 * - Concentration of oxygen
 * - Concentration of urea
 * - Volume fraction of biofilm
 * - Volume fraction of calcite
 * 
 * II) Biofilm effects in underground applications (e.g., hydrogen storage)
 * 
 * Biofilm-related effects in subsurface applications such as hydrogen
 * storage include reduced injectivity and hydrogen loss. The conceptual
 * model includes the following mechanisms:
 * 
 * - Biofilm is present in the storage site prior to injection.
 * - The biofilm consumes injected hydrogen/CO2, leading to clogging
 *   effects.
 * 
 * This implementation considers a two-phase (gas + water) system with the
 * following primary variables:
 * - Pressure of the gas phase
 * - Water saturation / dissolved gas in water
 * - Concentration of suspended microbes
 * - Volume fraction of biofilm
 */
template <class TypeTag, bool enableBioeffectsV = getPropValue<TypeTag, Properties::EnableBioeffects>()>
class BlackOilBioeffectsModule
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
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilBioeffectsParams<Scalar>::TabulatedFunction;

    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };

    static constexpr unsigned microbialConcentrationIdx = Indices::microbialConcentrationIdx;
    static constexpr unsigned oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
    static constexpr unsigned ureaConcentrationIdx = Indices::ureaConcentrationIdx;
    static constexpr unsigned biofilmVolumeFractionIdx = Indices::biofilmVolumeFractionIdx;
    static constexpr unsigned calciteVolumeFractionIdx = Indices::calciteVolumeFractionIdx;
    static constexpr unsigned contiMicrobialEqIdx = Indices::contiMicrobialEqIdx;
    static constexpr unsigned contiOxygenEqIdx = Indices::contiOxygenEqIdx;
    static constexpr unsigned contiUreaEqIdx = Indices::contiUreaEqIdx;
    static constexpr unsigned contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
    static constexpr unsigned contiCalciteEqIdx = Indices::contiCalciteEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableBioeffects = enableBioeffectsV;
    static constexpr bool enableMICP = Indices::enableMICP;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:
    //! \brief Set parameters.
    static void setParams(BlackOilBioeffectsParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil bioeffects module.
     */
    static void registerParameters()
    {
        if constexpr (enableBioeffects)
            VtkBlackOilBioeffectsModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all bioeffects specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if constexpr (enableBioeffects)
            model.addOutputModule(std::make_unique<VtkBlackOilBioeffectsModule<TypeTag>>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableBioeffects)
            if constexpr (enableMICP)
                return eqIdx == contiMicrobialEqIdx || eqIdx == contiOxygenEqIdx || eqIdx == contiUreaEqIdx 
                    || eqIdx == contiBiofilmEqIdx || eqIdx == contiCalciteEqIdx;
            else
               return eqIdx == contiMicrobialEqIdx || eqIdx == contiBiofilmEqIdx;
        else
            return false;
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
        if constexpr (enableBioeffects) {
            const auto& fs = intQuants.fluidState();
            LhsEval surfaceVolumeWater = Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx)) *
                                         Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx)) *
                                         Toolbox::template decay<LhsEval>(intQuants.porosity());
            // avoid singular matrix if no water is present
            surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);
            // suspended microbes in water phase
            const LhsEval accumulationMicrobes = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.microbialConcentration());
            storage[contiMicrobialEqIdx] += accumulationMicrobes;
            // biofilm
            const LhsEval accumulationBiofilm = Toolbox::template decay<LhsEval>(intQuants.biofilmVolumeFraction());
            storage[contiBiofilmEqIdx] += accumulationBiofilm;
            if constexpr (enableMICP) {
                // oxygen in water phase
                const LhsEval accumulationOxygen = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.oxygenConcentration());
                storage[contiOxygenEqIdx] += accumulationOxygen;
                // urea in water phase (applying the scaling factor for the urea equation)
                const LhsEval accumulationUrea = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.ureaConcentration());
                storage[contiUreaEqIdx] += accumulationUrea;
                storage[contiUreaEqIdx] *= getPropValue<TypeTag, Properties::BlackOilUreaScalingFactor>();
                // calcite
                const LhsEval accumulationCalcite = Toolbox::template decay<LhsEval>(intQuants.calciteVolumeFraction());
                storage[contiCalciteEqIdx] += accumulationCalcite;
            }
        }
    }

    template <class UpEval, class Eval, class IntensiveQuantities>
    static void addBioeffectsFluxes_(RateVector& flux,
                                     unsigned phaseIdx,
                                     const Eval& volumeFlux,
                                     const IntensiveQuantities& upFs)
    {
        if (phaseIdx == waterPhaseIdx) {
            if constexpr (enableBioeffects) {
                flux[contiMicrobialEqIdx] =
                    decay<UpEval>(upFs.microbialConcentration())
                    * volumeFlux;
                if constexpr (enableMICP) {
                    flux[contiOxygenEqIdx] =
                        decay<UpEval>(upFs.oxygenConcentration())
                        * volumeFlux;
                    flux[contiUreaEqIdx] =
                        decay<UpEval>(upFs.ureaConcentration())
                        * volumeFlux;
                }
            }
        }
    }

    // since the urea concentration can be much larger than 1, then we apply a scaling factor
    static void applyScaling(RateVector& flux)
    {
        if constexpr (enableMICP) {
            flux[contiUreaEqIdx] *= getPropValue<TypeTag, Properties::BlackOilUreaScalingFactor>();
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableBioeffects) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
            unsigned focusIdx = elemCtx.focusDofIndex();
            unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
            flux[contiMicrobialEqIdx] = 0.0;
            if constexpr (enableMICP) {
                flux[contiOxygenEqIdx] = 0.0;
                flux[contiUreaEqIdx] = 0.0;
            }
            if (upIdx == focusIdx)
                addBioeffectsFluxes_<Evaluation>(flux, elemCtx, scvfIdx, timeIdx);
            else
                addBioeffectsFluxes_<Scalar>(flux, elemCtx, scvfIdx, timeIdx);
        }
    }

    template <class UpstreamEval>
    static void addBioeffectsFluxes_(RateVector& flux,
                               const ElementContext& elemCtx,
                               unsigned scvfIdx,
                               unsigned timeIdx)
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = up.fluidState();
        const auto& volFlux = extQuants.volumeFlux(waterPhaseIdx);
        addBioeffectsFluxes_<UpstreamEval>(flux, volFlux, fs);
    }

    static void addSource(RateVector& source,
                          const Problem& problem,
                          const IntensiveQuantities& intQuants,
                          unsigned globalSpaceIdex)
    {
        if constexpr (enableBioeffects) {
            const auto b = intQuants.fluidState().invB(waterPhaseIdx);
            unsigned satnumIdx = problem.satnumRegionIndex(globalSpaceIdex);
            Scalar rho_b = densityBiofilm(satnumIdx);
            Scalar k_d = microbialDeathRate(satnumIdx);
            Scalar mu = maximumGrowthRate(satnumIdx);
            Scalar k_n = halfVelocityGrowth(satnumIdx);
            Scalar Y = yieldGrowthCoefficient(satnumIdx);
            Scalar k_a = microbialAttachmentRate(satnumIdx);
            Scalar k_str = detachmentRate(satnumIdx);
            Scalar eta = detachmentExponent(satnumIdx);
            const auto& velocityInf = problem.model().linearizer().getVelocityInfo();
            auto velocityInfos = velocityInf[globalSpaceIdex];
            Scalar normVelocityCell =
                std::accumulate(velocityInfos.begin(), velocityInfos.end(), 0.0,
                                [](const auto acc, const auto& info)
                                { return max(acc, std::abs(info.velocity[waterPhaseIdx])); });
            if constexpr (enableMICP) {
                Scalar rho_c = densityCalcite(satnumIdx);
                Scalar k_u = halfVelocityUrea(satnumIdx);
                Scalar mu_u = maximumUreaUtilization(satnumIdx);
                Scalar F = oxygenConsumptionFactor(satnumIdx);
                Scalar Y_uc = yieldUreaToCalciteCoefficient(satnumIdx);

                // compute Monod terms (the negative region is replaced by a straight line)
                // Schäfer et al (1998) https://doi.org/10.1016/S0169-7722(97)00060-0
                Evaluation k_g = mu * intQuants.oxygenConcentration() / (k_n + intQuants.oxygenConcentration());
                Evaluation k_c = mu_u * intQuants.ureaConcentration() / (k_u + intQuants.ureaConcentration());
                if (intQuants.oxygenConcentration() < 0) {
                    k_g = mu * intQuants.oxygenConcentration() / k_n;
                }
                if (intQuants.ureaConcentration() < 0) {
                    k_c = mu_u * intQuants.ureaConcentration() / k_u;
                }

                // compute the processes
                // see https://doi.org/10.1016/j.ijggc.2021.103256 for the MICP processes in the model
                source[Indices::contiMicrobialEqIdx] += intQuants.microbialConcentration() * intQuants.porosity() *
                                                        b * (Y * k_g - k_d - k_a) +
                                                        rho_b * intQuants.biofilmVolumeFraction() * k_str * pow(normVelocityCell / intQuants.porosity(), eta);

                source[Indices::contiOxygenEqIdx] -= (intQuants.microbialConcentration() * intQuants.porosity() * 
                                                    b + rho_b * intQuants.biofilmVolumeFraction()) * F * k_g;

                source[Indices::contiUreaEqIdx] -= rho_b * intQuants.biofilmVolumeFraction() * k_c;

                source[Indices::contiBiofilmEqIdx] += intQuants.biofilmVolumeFraction() * (Y * k_g - k_d -
                                                    k_str * pow(normVelocityCell / intQuants.porosity(), eta) - Y_uc * (rho_b / rho_c) * 
                                                    intQuants.biofilmVolumeFraction() * k_c / (intQuants.porosity() + 
                                                    intQuants.biofilmVolumeFraction())) + k_a * intQuants.microbialConcentration() *
                                                    intQuants.porosity() * b / rho_b;

                source[Indices::contiCalciteEqIdx] += (rho_b / rho_c) * intQuants.biofilmVolumeFraction() * Y_uc * k_c;
                
                // since the urea concentration can be much larger than 1, then we apply a scaling factor
                source[Indices::contiUreaEqIdx] *= getPropValue<TypeTag, Properties::BlackOilUreaScalingFactor>();
            }
            else {                
                const Scalar normVelocityCellG =
                std::accumulate(velocityInfos.begin(), velocityInfos.end(), 0.0,
                                [](const auto acc, const auto& info)
                                { return max(acc, std::abs(info.velocity[1])); });
                normVelocityCell = max(normVelocityCellG, normVelocityCell);
                // convert Rsw to concentration to use in source term
                const auto& fs = intQuants.fluidState();
                const auto& Sw = fs.saturation(waterPhaseIdx);
                const auto& Rsw = fs.Rsw();
                const auto& rhow = fs.density(waterPhaseIdx);
                unsigned pvtRegionIndex = fs.pvtRegionIndex();

                const auto& xG = RswToMassFraction(pvtRegionIndex, Rsw);

                // get the porosity and and gas density for convenience
                const Evaluation& poro = intQuants.porosity();
                Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);

                // calculate biofilm growth rate
                Evaluation k_g = mu * (xG * rhow * poro * Sw / (xG * rhow * poro * Sw + k_n));
                if (xG < 0) {
                    k_g = mu * (xG * rhow * poro * Sw / k_n);
                }

                // compute source terms
                // decay, detachment, and attachment rate of suspended microbes
                source[contiMicrobialEqIdx] += Sw * intQuants.microbialConcentration() * intQuants.porosity() * b
                                               * (- k_a - k_d)
                                               + intQuants.biofilmVolumeFraction() * rho_b * k_str
                                               * pow(normVelocityCell / intQuants.porosity(), eta);
                // biofilm growth and decay rate
                source[contiBiofilmEqIdx] += (k_g - k_d - k_str * pow(normVelocityCell / intQuants.porosity(), eta))
                                             * intQuants.biofilmVolumeFraction()
                                             + k_a * Sw * intQuants.microbialConcentration() * intQuants.porosity() * b / rho_b;

                // biofilm consumption of dissolved gas is proportional to biofilm growth rate
                unsigned activeGasCompIdx = FluidSystem::canonicalToActiveCompIdx(gasCompIdx);
                source[activeGasCompIdx] -= intQuants.biofilmVolumeFraction() * rho_b * k_g / (Y * rho_gRef);
            }
        }
    }

    static void addSource([[maybe_unused]] RateVector& source,
                          [[maybe_unused]] const ElementContext& elemCtx,
                          [[maybe_unused]] unsigned dofIdx,
                          [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableMICP) {
            const auto& problem = elemCtx.problem();
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
            addSource(source, problem, intQuants, dofIdx);
        }
    }

    static const Scalar densityBiofilm(unsigned satnumRegionIdx)
    {
        return params_.densityBiofilm_[satnumRegionIdx];
    }

    static const Scalar densityCalcite(unsigned satnumRegionIdx)
    {
        return params_.densityCalcite_[satnumRegionIdx];
    }

    static const Scalar detachmentRate(unsigned satnumRegionIdx)
    {
        return params_.detachmentRate_[satnumRegionIdx]; 
    }

    static const Scalar detachmentExponent(unsigned satnumRegionIdx)
    {
        return params_.detachmentExponent_[satnumRegionIdx];
    }

    static const Scalar halfVelocityGrowth(unsigned satnumRegionIdx)
    {
        return params_.halfVelocityGrowth_[satnumRegionIdx];
    }

    static const Scalar halfVelocityUrea(unsigned satnumRegionIdx)
    {
        return params_.halfVelocityUrea_[satnumRegionIdx];
    }

    static const Scalar maximumGrowthRate(unsigned satnumRegionIdx)
    {
        return params_.maximumGrowthRate_[satnumRegionIdx];
    }

    static const Scalar maximumUreaUtilization(unsigned satnumRegionIdx)
    {
        return params_.maximumUreaUtilization_[satnumRegionIdx];
    }

    static const Scalar microbialAttachmentRate(unsigned satnumRegionIdx)
    {
        return params_.microbialAttachmentRate_[satnumRegionIdx];
    }

    static const Scalar microbialDeathRate(unsigned satnumRegionIdx)
    {
        return params_.microbialDeathRate_[satnumRegionIdx];
    }

    static const Scalar oxygenConsumptionFactor(unsigned satnumRegionIdx)
    {
        return params_.oxygenConsumptionFactor_[satnumRegionIdx];
    }

    static const Scalar yieldGrowthCoefficient(unsigned satnumRegionIdx)
    {
        return params_.yieldGrowthCoefficient_[satnumRegionIdx];
    }

    static const Scalar yieldUreaToCalciteCoefficient(unsigned satnumRegionIdx)
    {
        return params_.yieldUreaToCalciteCoefficient_[satnumRegionIdx];
    }

    static const Scalar bioDiffCoefficient(unsigned pvtRegionIdx, unsigned compIdx)
    {
        return params_.bioDiffCoefficient_[pvtRegionIdx][compIdx];
    }

    static const TabulatedFunction& permfactTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.permfactTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& permfactTable(unsigned satnumRegionIdx)
    {
        return params_.permfactTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& pcfactTable(unsigned satnumRegionIdx)
    {
        return params_.pcfactTable_[satnumRegionIdx];
    }

    static bool hasPcfactTables()
    {
        if constexpr (enableBioeffects && !enableMICP) {
            return !params_.pcfactTable_.empty();
        }
        else {
            return false;
        }
    }

private:
    static BlackOilBioeffectsParams<Scalar> params_;

    static Evaluation RswToMassFraction(unsigned regionIdx, const Evaluation& Rsw) {
        Scalar rho_wRef = FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, regionIdx);
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);

        const Evaluation rho_oG = Rsw * rho_gRef;

        return rho_oG / (rho_wRef + rho_oG);
    }
};


template <class TypeTag, bool enableBioeffectsV>
BlackOilBioeffectsParams<typename BlackOilBioeffectsModule<TypeTag, enableBioeffectsV>::Scalar>
BlackOilBioeffectsModule<TypeTag, enableBioeffectsV>::params_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilBioeffectsIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        bioeffects extension of the black-oil model.
 */
template <class TypeTag, bool enableBioeffectsV = getPropValue<TypeTag, Properties::EnableBioeffects>()>
class BlackOilBioeffectsIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag>;

    static constexpr int microbialConcentrationIdx = Indices::microbialConcentrationIdx;
    static constexpr int oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
    static constexpr int ureaConcentrationIdx = Indices::ureaConcentrationIdx;
    static constexpr int biofilmVolumeFractionIdx = Indices::biofilmVolumeFractionIdx;
    static constexpr int calciteVolumeFractionIdx = Indices::calciteVolumeFractionIdx;
    static constexpr int temperatureIdx = Indices::temperatureIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool enableMICP = Indices::enableMICP;

public:

    /*!
     * \brief Update the intensive properties needed to handle bioeffects from the
     *        primary variables
     *
     */
    void bioeffectsPropertiesUpdate_(const ElementContext& elemCtx,
                                     unsigned dofIdx,
                                     unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const Scalar referencePorosity_ = elemCtx.problem().referencePorosity(dofIdx, timeIdx);
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);

        microbialConcentration_ = priVars.makeEvaluation(microbialConcentrationIdx, timeIdx, linearizationType);
        biofilmVolumeFraction_ = priVars.makeEvaluation(biofilmVolumeFractionIdx, timeIdx, linearizationType);
        biofilmMass_ = biofilmVolumeFraction_ * BioeffectsModule::densityBiofilm(satnumRegionIdx);
        calciteVolumeFraction_ = 0.0;
        if constexpr (enableMICP) {
            oxygenConcentration_ = priVars.makeEvaluation(oxygenConcentrationIdx, timeIdx, linearizationType);
            ureaConcentration_ = priVars.makeEvaluation(ureaConcentrationIdx, timeIdx, linearizationType);
            calciteVolumeFraction_ = priVars.makeEvaluation(calciteVolumeFractionIdx, timeIdx, linearizationType);
            calciteMass_ = calciteVolumeFraction_ * BioeffectsModule::densityCalcite(satnumRegionIdx);
        }
        const Evaluation poroFact = min(1.0 - (biofilmVolumeFraction_ + calciteVolumeFraction_) /
                                              (referencePorosity_), 1.0); //phi/phi_0

        const auto& permfactTable = BioeffectsModule::permfactTable(satnumRegionIdx);
        permFactor_ = permfactTable.eval(poroFact, /*extrapolation=*/true);
    }

    const Evaluation& microbialConcentration() const
    { return microbialConcentration_; }

    const Evaluation& oxygenConcentration() const
    { return oxygenConcentration_; }

    const Evaluation& ureaConcentration() const
    { return ureaConcentration_; }

    const Evaluation& biofilmVolumeFraction() const
    { return biofilmVolumeFraction_; }

    const Evaluation& calciteVolumeFraction() const
    { return calciteVolumeFraction_; }

    const Evaluation biofilmMass() const
    { return biofilmMass_; }

    const Evaluation calciteMass() const
    { return calciteMass_; }

    const Evaluation& permFactor() const
    { return permFactor_; }

protected:
    Evaluation microbialConcentration_;
    Evaluation oxygenConcentration_;
    Evaluation ureaConcentration_;
    Evaluation biofilmVolumeFraction_;
    Evaluation calciteVolumeFraction_;
    Evaluation biofilmMass_;
    Evaluation calciteMass_;
    Evaluation permFactor_;
    Evaluation pcFactor_;

};

template <class TypeTag>
class BlackOilBioeffectsIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void bioeffectsPropertiesUpdate_(const ElementContext&,
                                     unsigned,
                                     unsigned)
    {}

    const Evaluation& microbialConcentration() const
    { throw std::logic_error("microbialConcentration() called but MICP is disabled"); }

    const Evaluation& oxygenConcentration() const
    { throw std::logic_error("oxygenConcentration() called but MICP is disabled"); }

    const Evaluation& ureaConcentration() const
    { throw std::logic_error("ureaConcentration() called but MICP is disabled"); }

    const Evaluation& biofilmVolumeFraction() const
    { throw std::logic_error("biofilmVolumeFraction() called but biofilm/MICP is disabled"); }

    const Evaluation& calciteVolumeFraction() const
    { throw std::logic_error("calciteVolumeFraction() called but MICP is disabled"); }

    const Evaluation& biofilmMass() const
    { throw std::logic_error("biofilmMass() called but biofilm/MICP is disabled"); }

    const Evaluation& calciteMass() const
    { throw std::logic_error("calciteMass() called but MICP is disabled"); }

    const Evaluation& permFactor() const
    { throw std::logic_error("permFactor() called but biofilm/MICP is disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilBioeffectsExtensiveQuantities
 *
 * \brief Provides the bioeffects specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableBioeffectsV = getPropValue<TypeTag, Properties::EnableBioeffects>()>
class BlackOilBioeffectsExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
};

template <class TypeTag>
class BlackOilBioeffectsExtensiveQuantities<TypeTag, false>{};

} // namespace Opm

#endif
