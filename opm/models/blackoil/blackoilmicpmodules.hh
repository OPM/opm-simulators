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
 * \brief Contains the classes required to extend the black-oil model by MICP.
 */
#ifndef EWOMS_BLACK_OIL_MICP_MODULE_HH
#define EWOMS_BLACK_OIL_MICP_MODULE_HH

#include "blackoilproperties.hh"

#include <opm/models/blackoil/blackoilmicpparams.hh>
#include <opm/models/io/vtkblackoilmicpmodule.hh>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/MICPpara.hpp>
#endif

#include <dune/common/fvector.hh>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by MICP.
 */
template <class TypeTag, bool enableMICPV = getPropValue<TypeTag, Properties::EnableMICP>()>
class BlackOilMICPModule
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

    static constexpr unsigned microbialConcentrationIdx = Indices::microbialConcentrationIdx;
    static constexpr unsigned oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
    static constexpr unsigned ureaConcentrationIdx = Indices::ureaConcentrationIdx;
    static constexpr unsigned biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;
    static constexpr unsigned calciteConcentrationIdx = Indices::calciteConcentrationIdx;
    static constexpr unsigned contiMicrobialEqIdx = Indices::contiMicrobialEqIdx;
    static constexpr unsigned contiOxygenEqIdx = Indices::contiOxygenEqIdx;
    static constexpr unsigned contiUreaEqIdx = Indices::contiUreaEqIdx;
    static constexpr unsigned contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
    static constexpr unsigned contiCalciteEqIdx = Indices::contiCalciteEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableMICP = enableMICPV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:

#if HAVE_ECL_INPUT
    //
     //* \brief Initialize all internal data structures needed by the MICP module
     //
    static void initFromState(const EclipseState& eclState)
    {
        // some sanity checks: if MICP is enabled, the MICP keyword must be
        // present, if MICP is disabled the keyword must not be present.
        if (enableMICP && !eclState.runspec().micp()) {
            throw std::runtime_error("Non-trivial MICP treatment requested at compile time, but "
                                     "the deck does not contain the MICP keyword");
        }
        else if (!enableMICP && eclState.runspec().micp()) {
            throw std::runtime_error("MICP treatment disabled at compile time, but the deck "
                                     "contains the MICP keyword");
        }

        if (!eclState.runspec().micp())
            return; // MICP treatment is supposed to be disabled*/

        // initialize the objects which deal with the MICPpara keyword
        const auto& MICPpara = eclState.getMICPpara();
                setMICPpara(MICPpara.getDensityBiofilm(),
                           MICPpara.getDensityCalcite(),
                           MICPpara.getDetachmentRate(),
                           MICPpara.getCriticalPorosity(),
                           MICPpara.getFittingFactor(),
                           MICPpara.getHalfVelocityOxygen(),
                           MICPpara.getHalfVelocityUrea(),
                           MICPpara.getMaximumGrowthRate(),
                           MICPpara.getMaximumUreaUtilization(),
                           MICPpara.getMicrobialAttachmentRate(),
                           MICPpara.getMicrobialDeathRate(),
                           MICPpara.getMinimumPermeability(),
                           MICPpara.getOxygenConsumptionFactor(),
                           MICPpara.getYieldGrowthCoefficient(),
                           MICPpara.getMaximumOxygenConcentration(),
                           MICPpara.getMaximumUreaConcentration(),
                           MICPpara.getToleranceBeforeClogging());
        // obtain the porosity for the clamp in the blackoilnewtonmethod
        params_.phi_ = eclState.fieldProps().get_double("PORO");
    }
#endif

    /*!
     * \brief The simulator stops if "clogging" has been (almost) reached in any of the cells.
     *
     * I.e., porosity - biofilm - calcite < tol_clgg, where tol_clgg is a given tolerance. In the
     * implemented model a permebaility-porosity relatonship is used where a minimum
     * permeability value is reached if porosity - biofilm - calcite < phi_crit.
     */
    static void checkCloggingMICP(const Model& model, const Scalar phi, unsigned dofIdx)
    {
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/1)[dofIdx];
        if (phi - priVars[biofilmConcentrationIdx] - priVars[calciteConcentrationIdx] < toleranceBeforeClogging())
            throw std::logic_error("Clogging has been (almost) reached in at least one cell\n");
    }

    /*!
     * \brief Specify the MICP properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setMICPpara(const Scalar& densityBiofilm,
                            const Scalar& densityCalcite,
                            const Scalar& detachmentRate,
                            const Scalar& criticalPorosity,
                            const Scalar& fittingFactor,
                            const Scalar& halfVelocityOxygen,
                            const Scalar& halfVelocityUrea,
                            const Scalar& maximumGrowthRate,
                            const Scalar& maximumUreaUtilization,
                            const Scalar& microbialAttachmentRate,
                            const Scalar& microbialDeathRate,
                            const Scalar& minimumPermeability,
                            const Scalar& oxygenConsumptionFactor,
                            const Scalar& yieldGrowthCoefficient,
                            const Scalar& maximumOxygenConcentration,
                            const Scalar& maximumUreaConcentration,
                            const Scalar& toleranceBeforeClogging)
    {
        params_.densityBiofilm_ = densityBiofilm;
        params_.densityCalcite_ = densityCalcite;
        params_.detachmentRate_ = detachmentRate;
        params_.criticalPorosity_ = criticalPorosity;
        params_.fittingFactor_ = fittingFactor;
        params_.halfVelocityOxygen_ = halfVelocityOxygen;
        params_.halfVelocityUrea_ = halfVelocityUrea;
        params_.maximumGrowthRate_ = maximumGrowthRate;
        params_.maximumUreaUtilization_ = maximumUreaUtilization;
        params_.microbialAttachmentRate_ = microbialAttachmentRate;
        params_.microbialDeathRate_ = microbialDeathRate;
        params_.minimumPermeability_ = minimumPermeability;
        params_.oxygenConsumptionFactor_ = oxygenConsumptionFactor;
        params_.yieldGrowthCoefficient_ = yieldGrowthCoefficient;
        params_.maximumOxygenConcentration_ = maximumOxygenConcentration;
        params_.maximumUreaConcentration_ = maximumUreaConcentration;
        params_.toleranceBeforeClogging_ = toleranceBeforeClogging;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil MICP module.
     */
    static void registerParameters()
    {
        if (!enableMICP)
            // MICP has been disabled at compile time
            return;

        VtkBlackOilMICPModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all MICP specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableMICP)
            // MICP has been disabled at compile time
            return;

        model.addOutputModule(new VtkBlackOilMICPModule<TypeTag>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableMICP)
            return false;

        // All MICP components are true here
        return eqIdx == contiMicrobialEqIdx || eqIdx == contiOxygenEqIdx || eqIdx == contiUreaEqIdx || eqIdx == contiBiofilmEqIdx || eqIdx == contiCalciteEqIdx;
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
        if (!enableMICP)
            return;

        LhsEval surfaceVolumeWater = Toolbox::template decay<LhsEval>(intQuants.porosity());
        // avoid singular matrix if no water is present.
        surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);
        // Suspended microbes in water phase
        const LhsEval massMicrobes = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.microbialConcentration());
        LhsEval accumulationMicrobes = massMicrobes;
        storage[contiMicrobialEqIdx] += accumulationMicrobes;
        // Oxygen in water phase
        const LhsEval massOxygen = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.oxygenConcentration());
        LhsEval accumulationOxygen = massOxygen;
        storage[contiOxygenEqIdx] += accumulationOxygen;
        // Urea in water phase
        const LhsEval massUrea = surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.ureaConcentration());
        LhsEval accumulationUrea = massUrea;
        storage[contiUreaEqIdx] += accumulationUrea;
        // Biofilm
        const LhsEval massBiofilm = Toolbox::template decay<LhsEval>(intQuants.biofilmConcentration());
        LhsEval accumulationBiofilm = massBiofilm;
        storage[contiBiofilmEqIdx] += accumulationBiofilm;
        // Calcite
        const LhsEval massCalcite = Toolbox::template decay<LhsEval>(intQuants.calciteConcentration());
        LhsEval accumulationCalcite = massCalcite;
        storage[contiCalciteEqIdx] += accumulationCalcite;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableMICP)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);


        if (upIdx == inIdx) {
            flux[contiMicrobialEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * up.microbialConcentration();
            flux[contiOxygenEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * up.oxygenConcentration();
            flux[contiUreaEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * up.ureaConcentration();
        }
        else {
            flux[contiMicrobialEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * decay<Scalar>(up.microbialConcentration());
            flux[contiOxygenEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * decay<Scalar>(up.oxygenConcentration());
            flux[contiUreaEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * decay<Scalar>(up.ureaConcentration());
        }
    }

    // See https://doi.org/10.1016/j.ijggc.2021.103256 for the micp processes in the model.
    static void addSource(RateVector& source,
                            const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)

    {
        if (!enableMICP)
            return;

        // compute dpW (max norm of the pressure gradient in the cell center)
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& K = elemCtx.problem().intrinsicPermeability(elemCtx, dofIdx, 0);
        size_t numInteriorFaces = elemCtx.numInteriorFaces(timeIdx);
        Evaluation dpW = 0;
        for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
          const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
          unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
          const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
          const Evaluation& mobWater = up.mobility(waterPhaseIdx);

          // compute water velocity from flux
          Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / (K[0][0] * mobWater);
          dpW = std::max(dpW, abs(waterVolumeVelocity));
        }

        // get the model parameters
        Scalar k_a = microbialAttachmentRate();
        Scalar k_d = microbialDeathRate();
        Scalar rho_b = densityBiofilm();
        Scalar rho_c = densityCalcite();
        Scalar k_str = detachmentRate();
        Scalar k_o = halfVelocityOxygen();
        Scalar k_u = halfVelocityUrea() / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp)
        Scalar mu = maximumGrowthRate();
        Scalar mu_u = maximumUreaUtilization() / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp)
        Scalar Y_sb = yieldGrowthCoefficient();
        Scalar F = oxygenConsumptionFactor();
        Scalar Y_uc = 1.67 * 10; //Multiplying by scaling factor 10 (see WellInterface_impl.hpp)

        // compute the processes
        source[Indices::contiMicrobialEqIdx] += intQuants.microbialConcentration() * intQuants.porosity() *
                                                      (Y_sb * mu * intQuants.oxygenConcentration() / (k_o + intQuants.oxygenConcentration()) - k_d - k_a)
                                                + rho_b * intQuants.biofilmConcentration() * k_str * pow(intQuants.porosity() * dpW, 0.58);

        source[Indices::contiOxygenEqIdx] -= (intQuants.microbialConcentration() * intQuants.porosity() + rho_b * intQuants.biofilmConcentration()) *
                                                    F * mu * intQuants.oxygenConcentration() / (k_o + intQuants.oxygenConcentration());

        source[Indices::contiUreaEqIdx] -= rho_b * intQuants.biofilmConcentration() * mu_u * intQuants.ureaConcentration() / (k_u + intQuants.ureaConcentration());

        source[Indices::contiBiofilmEqIdx] += intQuants.biofilmConcentration() * (Y_sb * mu * intQuants.oxygenConcentration() / (k_o + intQuants.oxygenConcentration()) - k_d
                                                              - k_str * pow(intQuants.porosity() * dpW, 0.58) - Y_uc * (rho_b / rho_c) * intQuants.biofilmConcentration() * mu_u *
                                                                  (intQuants.ureaConcentration() / (k_u + intQuants.ureaConcentration())) / (intQuants.porosity() + intQuants.biofilmConcentration()))
                                              + k_a * intQuants.microbialConcentration() * intQuants.porosity() / rho_b;

        source[Indices::contiCalciteEqIdx] += (rho_b / rho_c) * intQuants.biofilmConcentration() * Y_uc * mu_u * intQuants.ureaConcentration() / (k_u + intQuants.ureaConcentration());
    }

    static const Scalar densityBiofilm()
    {
        return params_.densityBiofilm_;
    }

    static const Scalar densityCalcite()
    {
        return params_.densityCalcite_;
    }

    static const Scalar detachmentRate()
    {
        return params_.detachmentRate_;
    }

    static const Scalar criticalPorosity()
    {
        return params_.criticalPorosity_;
    }

    static const Scalar fittingFactor()
    {
        return params_.fittingFactor_;
    }

    static const Scalar halfVelocityOxygen()
    {
        return params_.halfVelocityOxygen_;
    }

    static const Scalar halfVelocityUrea()
    {
        return params_.halfVelocityUrea_;
    }

    static const Scalar maximumGrowthRate()
    {
        return params_.maximumGrowthRate_;
    }

    static const Scalar maximumOxygenConcentration()
    {
        return params_.maximumOxygenConcentration_;
    }

    static const Scalar maximumUreaConcentration()
    {
        return params_.maximumUreaConcentration_ / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp);
    }

    static const Scalar maximumUreaUtilization()
    {
        return params_.maximumUreaUtilization_;
    }

    static const Scalar microbialAttachmentRate()
    {
        return params_.microbialAttachmentRate_;
    }

    static const Scalar microbialDeathRate()
    {
        return params_.microbialDeathRate_;
    }

    static const Scalar minimumPermeability()
    {
        return params_.minimumPermeability_;
    }

    static const Scalar oxygenConsumptionFactor()
    {
        return params_.oxygenConsumptionFactor_;
    }

    static const Scalar toleranceBeforeClogging()
    {
        return params_.toleranceBeforeClogging_;
    }

    static const Scalar yieldGrowthCoefficient()
    {
        return params_.yieldGrowthCoefficient_;
    }

    static const std::vector<Scalar> phi()
    {
        return params_.phi_;
    }

private:
    static BlackOilMICPParams<Scalar> params_;
};


template <class TypeTag, bool enableMICPV>
BlackOilMICPParams<typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar>
BlackOilMICPModule<TypeTag, enableMICPV>::params_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilMICPIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        MICP extension of the black-oil model.
 */
template <class TypeTag, bool enableMICPV = getPropValue<TypeTag, Properties::EnableMICP>()>
class BlackOilMICPIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using MICPModule = BlackOilMICPModule<TypeTag>;

    static constexpr int microbialConcentrationIdx = Indices::microbialConcentrationIdx;
    static constexpr int oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
    static constexpr int ureaConcentrationIdx = Indices::ureaConcentrationIdx;
    static constexpr int biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;
    static constexpr int calciteConcentrationIdx = Indices::calciteConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;


public:

    /*!
     * \brief Update the intensive properties needed to handle MICP from the
     *        primary variables
     *
     */
    void MICPPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& K = elemCtx.problem().intrinsicPermeability(elemCtx, dofIdx, timeIdx);
        Scalar referencePorosity_ = elemCtx.problem().porosity(elemCtx, dofIdx, timeIdx);
        Scalar eta = MICPModule::fittingFactor();
        Scalar k_min = MICPModule::minimumPermeability();
        Scalar phi_crit = MICPModule::criticalPorosity();

        microbialConcentration_ = priVars.makeEvaluation(microbialConcentrationIdx, timeIdx, linearizationType);
        oxygenConcentration_ = priVars.makeEvaluation(oxygenConcentrationIdx, timeIdx, linearizationType);
        ureaConcentration_ = priVars.makeEvaluation(ureaConcentrationIdx, timeIdx, linearizationType);
        biofilmConcentration_ = priVars.makeEvaluation(biofilmConcentrationIdx, timeIdx, linearizationType);
        calciteConcentration_ = priVars.makeEvaluation(calciteConcentrationIdx, timeIdx, linearizationType);

        // Permeability reduction due to MICP, by adjusting the water mobility
        asImp_().mobility_[waterPhaseIdx] *= max((pow((intQuants.porosity() - phi_crit) / (referencePorosity_ - phi_crit), eta) + k_min / K[0][0])/(1. + k_min / K[0][0]), k_min / K[0][0]);

    }

    const Evaluation& microbialConcentration() const
    { return microbialConcentration_; }

    const Evaluation& oxygenConcentration() const
    { return oxygenConcentration_; }

    const Evaluation& ureaConcentration() const
    { return ureaConcentration_; }

    const Evaluation& biofilmConcentration() const
    { return biofilmConcentration_; }

    const Evaluation& calciteConcentration() const
    { return calciteConcentration_; }


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation microbialConcentration_;
    Evaluation oxygenConcentration_;
    Evaluation ureaConcentration_;
    Evaluation biofilmConcentration_;
    Evaluation calciteConcentration_;

};

template <class TypeTag>
class BlackOilMICPIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void MICPPropertiesUpdate_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    { }

    const Evaluation& microbialConcentration() const
    { throw std::logic_error("microbialConcentration() called but MICP is disabled"); }

    const Evaluation& oxygenConcentration() const
    { throw std::logic_error("oxygenConcentration() called but MICP is disabled"); }

    const Evaluation& ureaConcentration() const
    { throw std::logic_error("ureaConcentration() called but MICP is disabled"); }

    const Evaluation& biofilmConcentration() const
    { throw std::logic_error("biofilmConcentration() called but MICP is disabled"); }

    const Evaluation& calciteConcentration() const
    { throw std::logic_error("calciteConcentration() called but MICP is disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilMICPExtensiveQuantities
 *
 * \brief Provides the MICP specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableMICPV = getPropValue<TypeTag, Properties::EnableMICP>()>
class BlackOilMICPExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilMICPExtensiveQuantities<TypeTag, false>{};

} // namespace Opm

#endif
