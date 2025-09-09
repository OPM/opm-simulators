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
 * \brief Contains the classes required to extend the black-oil model by polymer.
 */
#ifndef EWOMS_BLACK_OIL_POLYMER_MODULE_HH
#define EWOMS_BLACK_OIL_POLYMER_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/blackoil/blackoilpolymerparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/io/vtkblackoilpolymermodule.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by polymer.
 */
template <class TypeTag, bool enablePolymerV = getPropValue<TypeTag, Properties::EnablePolymer>()>
class BlackOilPolymerModule
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

    using TabulatedFunction = typename BlackOilPolymerParams<Scalar>::TabulatedFunction;
    using TabulatedTwoDFunction = typename BlackOilPolymerParams<Scalar>::TabulatedTwoDFunction;

    static constexpr unsigned polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr unsigned polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
    static constexpr unsigned contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr unsigned contiPolymerMolarWeightEqIdx = Indices::contiPolymerMWEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enablePolymer = enablePolymerV;
    static constexpr bool enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>();

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    //! \brief Set parameters.
    static void setParams(BlackOilPolymerParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
    * \brief get the PLYMWINJ table
    */
    static TabulatedTwoDFunction& getPlymwinjTable(const int tableNumber)
    {
        const auto iterTable = params_.plymwinjTables_.find(tableNumber);
        if (iterTable != params_.plymwinjTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the PLYMWINJ table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRWAT table
    */
    static TabulatedTwoDFunction& getSkprwatTable(const int tableNumber)
    {
        const auto iterTable = params_.skprwatTables_.find(tableNumber);
        if (iterTable != params_.skprwatTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRWAT table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRPOLY table
    */
    static typename BlackOilPolymerParams<Scalar>::SkprpolyTable&
    getSkprpolyTable(const int tableNumber)
    {
        const auto iterTable = params_.skprpolyTables_.find(tableNumber);
        if (iterTable != params_.skprpolyTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRPOLY table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
     * \brief Register all run-time parameters for the black-oil polymer module.
     */
    static void registerParameters()
    {
        if constexpr (enablePolymer) {
            VtkBlackOilPolymerModule<TypeTag>::registerParameters();
        }
    }

    /*!
     * \brief Register all polymer specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if constexpr (enablePolymer) {
            model.addOutputModule(std::make_unique<VtkBlackOilPolymerModule<TypeTag>>(simulator));
        }
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if constexpr (enablePolymer) {
            if constexpr (enablePolymerMolarWeight) {
                return pvIdx == polymerConcentrationIdx || pvIdx == polymerMoleWeightIdx;
            }
            else {
                return pvIdx == polymerConcentrationIdx;
            }
        }
        else {
            return false;
        }
    }

    static std::string primaryVarName(unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        if (pvIdx == polymerConcentrationIdx) {
            return "polymer_waterconcentration";
        }
        else {
            return "polymer_molecularweight";
        }
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enablePolymer) {
            if constexpr (enablePolymerMolarWeight) {
                return eqIdx == contiPolymerEqIdx || eqIdx == contiPolymerMolarWeightEqIdx;
            }
            else {
                return eqIdx == contiPolymerEqIdx;
            }
        }
        else {
            return false;
        }
    }

    static std::string eqName(unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        if (eqIdx == contiPolymerEqIdx) {
            return "conti^polymer";
        }
        else {
            return "conti^polymer_molecularweight";
        }
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
        if constexpr (enablePolymer) {
            const auto& fs = intQuants.fluidState();

            // avoid singular matrix if no water is present.
            const LhsEval surfaceVolumeWater =
                    max(Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(intQuants.porosity()),
                        1e-10);

            // polymer in water phase
            const LhsEval massPolymer =
                surfaceVolumeWater *
                Toolbox::template decay<LhsEval>(intQuants.polymerConcentration()) *
                (1.0 - Toolbox::template decay<LhsEval>(intQuants.polymerDeadPoreVolume()));

            // polymer in solid phase
            const LhsEval adsorptionPolymer =
                    Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity()) *
                    Toolbox::template decay<LhsEval>(intQuants.polymerRockDensity()) *
                    Toolbox::template decay<LhsEval>(intQuants.polymerAdsorption());

            LhsEval accumulationPolymer = massPolymer + adsorptionPolymer;

            storage[contiPolymerEqIdx] += accumulationPolymer;

            // tracking the polymer molecular weight
            if constexpr (enablePolymerMolarWeight) {
                accumulationPolymer = max(accumulationPolymer, 1e-10);

                storage[contiPolymerMolarWeightEqIdx]  +=
                    accumulationPolymer * Toolbox::template decay<LhsEval>(intQuants.polymerMoleWeight());
            }
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enablePolymer) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

            const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
            const unsigned inIdx = extQuants.interiorIndex();
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            const unsigned contiWaterEqIdx =
                Indices::conti0EqIdx + FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx);

            if (upIdx == inIdx) {
                flux[contiPolymerEqIdx] =
                        extQuants.volumeFlux(waterPhaseIdx) *
                        up.fluidState().invB(waterPhaseIdx) *
                        up.polymerViscosityCorrection() /
                        extQuants.polymerShearFactor() *
                        up.polymerConcentration();

                // modify water
                flux[contiWaterEqIdx] /= extQuants.waterShearFactor();
            }
            else {
                flux[contiPolymerEqIdx] =
                        extQuants.volumeFlux(waterPhaseIdx) *
                        decay<Scalar>(up.fluidState().invB(waterPhaseIdx)) *
                        decay<Scalar>(up.polymerViscosityCorrection()) /
                        decay<Scalar>(extQuants.polymerShearFactor()) *
                        decay<Scalar>(up.polymerConcentration());

                // modify water
                flux[contiWaterEqIdx] /= decay<Scalar>(extQuants.waterShearFactor());
            }

            // flux related to transport of polymer molecular weight
            if constexpr (enablePolymerMolarWeight) {
                if (upIdx == inIdx) {
                    flux[contiPolymerMolarWeightEqIdx] =
                        flux[contiPolymerEqIdx] * up.polymerMoleWeight();
                }
                else {
                    flux[contiPolymerMolarWeightEqIdx] =
                        flux[contiPolymerEqIdx] * decay<Scalar>(up.polymerMoleWeight());
                }
            }
        }
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables&,
                                     const EqVector&)
    {
        // do not consider consider the cange of polymer primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if constexpr (enablePolymer) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
            outstream << priVars[polymerConcentrationIdx];
            outstream << priVars[polymerMoleWeightIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if constexpr (enablePolymer) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
            PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

            instream >> priVars0[polymerConcentrationIdx];
            instream >> priVars0[polymerMoleWeightIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1[polymerConcentrationIdx] = priVars0[polymerConcentrationIdx];
            priVars1[polymerMoleWeightIdx] = priVars0[polymerMoleWeightIdx];
        }
    }

    static Scalar plyrockDeadPoreVolume(const ElementContext& elemCtx,
                                        unsigned scvIdx,
                                        unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyrockDeadPoreVolume_[satnumRegionIdx];
    }

    static Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyrockResidualResistanceFactor_[satnumRegionIdx];
    }

    static Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
                                           unsigned scvIdx,
                                           unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyrockRockDensityFactor_[satnumRegionIdx];
    }

    static Scalar plyrockAdsorbtionIndex(const ElementContext& elemCtx,
                                         unsigned scvIdx,
                                         unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyrockAdsorbtionIndex_[satnumRegionIdx];
    }

    static Scalar plyrockMaxAdsorbtion(const ElementContext& elemCtx,
                                       unsigned scvIdx,
                                       unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyrockMaxAdsorbtion_[satnumRegionIdx];
    }

    static const TabulatedFunction&
    plyadsAdsorbedPolymer(const ElementContext& elemCtx,
                          unsigned scvIdx,
                          unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyadsAdsorbedPolymer_[satnumRegionIdx];
    }

    static const TabulatedFunction&
    plyviscViscosityMultiplierTable(const ElementContext& elemCtx,
                                    unsigned scvIdx,
                                    unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx =
            elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const TabulatedFunction&
    plyviscViscosityMultiplierTable(unsigned pvtnumRegionIdx)
    { return params_.plyviscViscosityMultiplierTable_[pvtnumRegionIdx]; }

    static Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
                                         unsigned scvIdx,
                                         unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx =
            elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plymaxMaxConcentration_[polymerMixRegionIdx];
    }

    static Scalar plymixparToddLongstaff(const ElementContext& elemCtx,
                                         unsigned scvIdx,
                                         unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx =
            elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plymixparToddLongstaff_[polymerMixRegionIdx];
    }

    static const typename BlackOilPolymerParams<Scalar>::PlyvmhCoefficients&
    plyvmhCoefficients(const ElementContext& elemCtx,
                       const unsigned scvIdx,
                       const unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx =
            elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.plyvmhCoefficients_[polymerMixRegionIdx];
    }

    static bool hasPlyshlog()
    { return params_.hasPlyshlog_; }

    static bool hasShrate()
    { return params_.hasShrate_; }

    static Scalar shrate(unsigned pvtnumRegionIdx)
    { return params_.shrate_[pvtnumRegionIdx]; }

    /*!
     * \brief Computes the shear factor
     *
     * Input is polymer concentration and either the water velocity or the shrate if hasShrate_ is true.
     * The pvtnumRegionIdx is needed to make sure the right table is used.
     */
    template <class Evaluation>
    static Evaluation computeShearFactor(const Evaluation& polymerConcentration,
                                         unsigned pvtnumRegionIdx,
                                         const Evaluation& v0)
    {
        using ToolboxLocal = MathToolbox<Evaluation>;

        const auto& viscosityMultiplierTable = params_.plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
        const Scalar viscosityMultiplier =
            viscosityMultiplierTable.eval(scalarValue(polymerConcentration), /*extrapolate=*/true);

        const Scalar eps = 1e-14;
        // return 1.0 if the polymer has no effect on the water.
        if (std::abs((viscosityMultiplier - 1.0)) < eps) {
            return ToolboxLocal::createConstant(v0, 1.0);
        }

        const std::vector<Scalar>& shearEffectRefLogVelocity =
            params_.plyshlogShearEffectRefLogVelocity_[pvtnumRegionIdx];
        const auto v0AbsLog = log(abs(v0));
        // return 1.0 if the velocity /sharte is smaller than the first velocity entry.
        if (v0AbsLog < shearEffectRefLogVelocity[0]) {
            return ToolboxLocal::createConstant(v0, 1.0);
        }

        // compute shear factor from input
        // Z = (1 + (P - 1) * M(v)) / P
        // where M(v) is computed from user input
        // and P = viscosityMultiplier
        const std::vector<Scalar>& shearEffectRefMultiplier =
            params_.plyshlogShearEffectRefMultiplier_[pvtnumRegionIdx];
        const std::size_t numTableEntries = shearEffectRefLogVelocity.size();
        assert(shearEffectRefMultiplier.size() == numTableEntries);

        std::vector<Scalar> shearEffectMultiplier(numTableEntries, 1.0);
        for (std::size_t i = 0; i < numTableEntries; ++i) {
            shearEffectMultiplier[i] = (1.0 + (viscosityMultiplier - 1.0) *
                                               shearEffectRefMultiplier[i]) / viscosityMultiplier;
            shearEffectMultiplier[i] = log(shearEffectMultiplier[i]);
        }
        // store the logarithmic velocity and logarithmic multipliers in a table for easy look up and
        // linear interpolation in the logarithmic space.
        const TabulatedFunction logShearEffectMultiplier =
            TabulatedFunction(numTableEntries, shearEffectRefLogVelocity,
                              shearEffectMultiplier, /*bool sortInputs =*/ false);

        // Find sheared velocity (v) that satisfies
        // F = log(v) + log (Z) - log(v0) = 0;

        // Set up the function
        // u = log(v)
        auto F = [&logShearEffectMultiplier, &v0AbsLog](const Evaluation& u) {
            return u + logShearEffectMultiplier.eval(u, true) - v0AbsLog;
        };
        // and its derivative
        auto dF = [&logShearEffectMultiplier](const Evaluation& u) {
            return 1 + logShearEffectMultiplier.evalDerivative(u, true);
        };

        // Solve F = 0 using Newton
        // Use log(v0) as initial value for u
        auto u = v0AbsLog;
        bool converged = false;
        // TODO make this into parameters
        for (int i = 0; i < 20; ++i) {
            const auto f = F(u);
            const auto df = dF(u);
            u -= f / df;
            if (std::abs(scalarValue(f)) < 1e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Not able to compute shear velocity. \n");
        }

        // return the shear factor
        return exp(logShearEffectMultiplier.eval(u, /*extrapolate=*/true));
    }

    Scalar molarMass() const
    { return 0.25; } // kg/mol

private:
    static BlackOilPolymerParams<Scalar> params_;
};

template <class TypeTag, bool enablePolymerV>
BlackOilPolymerParams<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::params_;

template <class TypeTag, bool enablePolymerV>
class BlackOilPolymerIntensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag>
class BlackOilPolymerIntensiveQuantities<TypeTag, /*enablePolymerV=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using PolymerModule = BlackOilPolymerModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;

public:
    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void polymerPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        polymerConcentration_ = priVars.makeEvaluation(polymerConcentrationIdx, timeIdx, linearizationType);
        if constexpr (enablePolymerMolarWeight) {
            polymerMoleWeight_ = priVars.makeEvaluation(polymerMoleWeightIdx, timeIdx, linearizationType);
        }

        // permeability reduction due to polymer
        const Scalar& maxAdsorbtion = PolymerModule::plyrockMaxAdsorbtion(elemCtx, dofIdx, timeIdx);
        const auto& plyadsAdsorbedPolymer = PolymerModule::plyadsAdsorbedPolymer(elemCtx, dofIdx, timeIdx);
        polymerAdsorption_ = plyadsAdsorbedPolymer.eval(polymerConcentration_, /*extrapolate=*/true);
        if (static_cast<int>(PolymerModule::plyrockAdsorbtionIndex(elemCtx, dofIdx, timeIdx)) ==
            BlackOilPolymerParams<Scalar>::NoDesorption)
        {
            const auto maxPolymerAdsorption =
                elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
            polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption), polymerAdsorption_);
        }

        // compute resitanceFactor
        const Scalar& residualResistanceFactor =
            PolymerModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) *
                                                   polymerAdsorption_ / maxAdsorbtion;

        // compute effective viscosities
        if constexpr (!enablePolymerMolarWeight) {
            const Scalar cmax = PolymerModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);
            const auto& fs = asImp_().fluidState_;
            const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
            const auto& viscosityMultiplier =
                PolymerModule::plyviscViscosityMultiplierTable(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityMixture =
                viscosityMultiplier.eval(polymerConcentration_, /*extrapolate=*/true) * muWater;

            // Do the Todd-Longstaff mixing
            const Scalar plymixparToddLongstaff = PolymerModule::plymixparToddLongstaff(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityPolymer = viscosityMultiplier.eval(cmax, /*extrapolate=*/true) * muWater;
            const Evaluation viscosityPolymerEffective =
                pow(viscosityMixture, plymixparToddLongstaff) * pow(viscosityPolymer, 1.0 - plymixparToddLongstaff);
            const Evaluation viscosityWaterEffective =
                pow(viscosityMixture, plymixparToddLongstaff) * pow(muWater, 1.0 - plymixparToddLongstaff);

            const Evaluation cbar = polymerConcentration_ / cmax;
            // waterViscosity / effectiveWaterViscosity
            waterViscosityCorrection_ = muWater * ((1.0 - cbar) / viscosityWaterEffective + cbar / viscosityPolymerEffective);
            // effectiveWaterViscosity / effectivePolymerViscosity
            polymerViscosityCorrection_ =  (muWater / waterViscosityCorrection_) / viscosityPolymerEffective;
        }
        else { // based on PLYVMH
            const auto& plyvmhCoefficients = PolymerModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
            const Scalar k_mh = plyvmhCoefficients.k_mh;
            const Scalar a_mh = plyvmhCoefficients.a_mh;
            const Scalar gamma = plyvmhCoefficients.gamma;
            const Scalar kappa = plyvmhCoefficients.kappa;

            // viscosity model based on Mark-Houwink equation and Huggins equation
            // 1000 is a emperical constant, most likely related to unit conversion
            const Evaluation intrinsicViscosity = k_mh * pow(polymerMoleWeight_ * 1000., a_mh);
            const Evaluation x = polymerConcentration_ * intrinsicViscosity;
            waterViscosityCorrection_ = 1.0 / (1.0 + gamma * (x + kappa * x * x));
            polymerViscosityCorrection_ = 1.0;
        }

        // adjust water mobility
        asImp_().mobility_[waterPhaseIdx] *= waterViscosityCorrection_ / resistanceFactor;

        // update rock properties
        polymerDeadPoreVolume_ = PolymerModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        polymerRockDensity_ = PolymerModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);
    }

    const Evaluation& polymerConcentration() const
    { return polymerConcentration_; }

    const Evaluation& polymerMoleWeight() const
    {
        if constexpr (enablePolymerMolarWeight) {
            return polymerMoleWeight_;
        }
        else {
            throw std::logic_error("polymerMoleWeight() is called but polymer milecular weight is disabled");
        }
    }

    Scalar polymerDeadPoreVolume() const
    { return polymerDeadPoreVolume_; }

    const Evaluation& polymerAdsorption() const
    { return polymerAdsorption_; }

    Scalar polymerRockDensity() const
    { return polymerRockDensity_; }

    // effectiveWaterViscosity / effectivePolymerViscosity
    const Evaluation& polymerViscosityCorrection() const
    { return polymerViscosityCorrection_; }

    // waterViscosity / effectiveWaterViscosity
    const Evaluation& waterViscosityCorrection() const
    { return waterViscosityCorrection_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerConcentration_;
    // polymer molecular weight
    Evaluation polymerMoleWeight_;
    Scalar polymerDeadPoreVolume_;
    Scalar polymerRockDensity_;
    Evaluation polymerAdsorption_;
    Evaluation polymerViscosityCorrection_;
    Evaluation waterViscosityCorrection_;
};

template <class TypeTag>
class BlackOilPolymerIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void polymerPropertiesUpdate_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    {}

    const Evaluation& polymerMoleWeight() const
    { throw std::logic_error("polymerMoleWeight() called but polymer molecular weight is disabled"); }

    const Evaluation& polymerConcentration() const
    { throw std::runtime_error("polymerConcentration() called but polymers are disabled"); }

    const Evaluation& polymerDeadPoreVolume() const
    { throw std::runtime_error("polymerDeadPoreVolume() called but polymers are disabled"); }

    const Evaluation& polymerAdsorption() const
    { throw std::runtime_error("polymerAdsorption() called but polymers are disabled"); }

    const Evaluation& polymerRockDensity() const
    { throw std::runtime_error("polymerRockDensity() called but polymers are disabled"); }

    const Evaluation& polymerViscosityCorrection() const
    { throw std::runtime_error("polymerViscosityCorrection() called but polymers are disabled"); }

    const Evaluation& waterViscosityCorrection() const
    { throw std::runtime_error("waterViscosityCorrection() called but polymers are disabled"); }
};

template <class TypeTag, bool enablePolymerV>
class BlackOilPolymerExtensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerExtensiveQuantities
 *
 * \brief Provides the polymer specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag>
class BlackOilPolymerExtensiveQuantities<TypeTag, /*enablePolymerV=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr unsigned waterPhaseIdx =  FluidSystem::waterPhaseIdx;

    using Toolbox = MathToolbox<Evaluation>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;

public:
    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using permeabilities, i.e., *not* via transmissibilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateShearMultipliersPerm(const ElementContext&,
                                    unsigned,
                                    unsigned)
    {
        throw std::runtime_error("The extension of the blackoil model for polymers is not yet "
                                 "implemented for problems specified using permeabilities.");
    }

    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using transmissibilities, i.e., *not* via permeabilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateShearMultipliers(const ElementContext& elemCtx,
                                unsigned scvfIdx,
                                unsigned timeIdx)
    {
        waterShearFactor_ = 1.0;
        polymerShearFactor_ = 1.0;

        if (!PolymerModule::hasPlyshlog()) {
            return;
        }

        const ExtensiveQuantities& extQuants = asImp_();
        const unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const unsigned interiorDofIdx = extQuants.interiorIndex();
        const unsigned exteriorDofIdx = extQuants.exteriorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        // compute water velocity from flux
        const Evaluation poroAvg = intQuantsIn.porosity() * 0.5 + intQuantsEx.porosity() * 0.5;
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvfIdx, timeIdx);
        const Evaluation& Sw = up.fluidState().saturation(waterPhaseIdx);
        const unsigned cellIdx = elemCtx.globalSpaceIndex(scvfIdx, timeIdx);
        const auto& materialLawManager = elemCtx.problem().materialLawManager();
        const auto& scaledDrainageInfo =
                materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);
        const Scalar& Swcr = scaledDrainageInfo.Swcr;

        // guard against zero porosity and no mobile water
        const Evaluation denom = max(poroAvg * (Sw - Swcr), 1e-12);
        Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / denom;

        // if shrate is specified. Compute shrate based on the water velocity
        if (PolymerModule::hasShrate()) {
            const Evaluation& relWater = up.relativePermeability(waterPhaseIdx);
            const Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
            if (trans > 0.0) {
                const Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
                const auto dist = elemCtx.pos(interiorDofIdx, timeIdx) -  elemCtx.pos(exteriorDofIdx, timeIdx);
                // compute permeability from transmissibility.
                const Scalar absPerm = trans / faceArea * dist.two_norm();
                waterVolumeVelocity *=
                    PolymerModule::shrate(pvtnumRegionIdx) * sqrt(poroAvg * Sw / (relWater * absPerm));
                assert(isfinite(waterVolumeVelocity));
            }
        }

        // compute share factors for water and polymer
        waterShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity);
        polymerShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity * up.polymerViscosityCorrection());
    }

    const Evaluation& polymerShearFactor() const
    { return polymerShearFactor_; }

    const Evaluation& waterShearFactor() const
    { return waterShearFactor_; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerShearFactor_;
    Evaluation waterShearFactor_;
};

template <class TypeTag>
class BlackOilPolymerExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    void updateShearMultipliers(const ElementContext&,
                                unsigned,
                                unsigned)
    {}

    void updateShearMultipliersPerm(const ElementContext&,
                                    unsigned,
                                    unsigned)
    {}

    const Evaluation& polymerShearFactor() const
    { throw std::runtime_error("polymerShearFactor() called but polymers are disabled"); }

    const Evaluation& waterShearFactor() const
    { throw std::runtime_error("waterShearFactor() called but polymers are disabled"); }
};

} // namespace Opm

#endif
