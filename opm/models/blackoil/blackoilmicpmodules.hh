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
        phi_ = eclState.fieldProps().get_double("PORO");
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
        if (phi - priVars[biofilmConcentrationIdx] - priVars[calciteConcentrationIdx] < MICPparaToleranceBeforeClogging())
            throw std::logic_error("Clogging has been (almost) reached in at least one cell\n");
    }

    /*!
     * \brief Specify the MICP properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setMICPpara(const Scalar& MICPparaDensityBiofilm,
                            const Scalar& MICPparaDensityCalcite,
                            const Scalar& MICPparaDetachmentRate,
                            const Scalar& MICPparaCriticalPorosity,
                            const Scalar& MICPparaFittingFactor,
                            const Scalar& MICPparaHalfVelocityOxygen,
                            const Scalar& MICPparaHalfVelocityUrea,
                            const Scalar& MICPparaMaximumGrowthRate,
                            const Scalar& MICPparaMaximumUreaUtilization,
                            const Scalar& MICPparaMicrobialAttachmentRate,
                            const Scalar& MICPparaMicrobialDeathRate,
                            const Scalar& MICPparaMinimumPermeability,
                            const Scalar& MICPparaOxygenConsumptionFactor,
                            const Scalar& MICPparaYieldGrowthCoefficient,
                            const Scalar& MICPparaMaximumOxygenConcentration,
                            const Scalar& MICPparaMaximumUreaConcentration,
                            const Scalar& MICPparaToleranceBeforeClogging)
    {
        MICPparaDensityBiofilm_ = MICPparaDensityBiofilm;
        MICPparaDensityCalcite_ = MICPparaDensityCalcite;
        MICPparaDetachmentRate_ = MICPparaDetachmentRate;
        MICPparaCriticalPorosity_ = MICPparaCriticalPorosity;
        MICPparaFittingFactor_ = MICPparaFittingFactor;
        MICPparaHalfVelocityOxygen_ = MICPparaHalfVelocityOxygen;
        MICPparaHalfVelocityUrea_ = MICPparaHalfVelocityUrea;
        MICPparaMaximumGrowthRate_ = MICPparaMaximumGrowthRate;
        MICPparaMaximumUreaUtilization_ = MICPparaMaximumUreaUtilization;
        MICPparaMicrobialAttachmentRate_ = MICPparaMicrobialAttachmentRate;
        MICPparaMicrobialDeathRate_ = MICPparaMicrobialDeathRate;
        MICPparaMinimumPermeability_ = MICPparaMinimumPermeability;
        MICPparaOxygenConsumptionFactor_ = MICPparaOxygenConsumptionFactor;
        MICPparaYieldGrowthCoefficient_ = MICPparaYieldGrowthCoefficient;
        MICPparaMaximumOxygenConcentration_ = MICPparaMaximumOxygenConcentration;
        MICPparaMaximumUreaConcentration_ = MICPparaMaximumUreaConcentration;
        MICPparaToleranceBeforeClogging_ = MICPparaToleranceBeforeClogging;
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
        Scalar k_a = MICPparaMicrobialAttachmentRate();
        Scalar k_d = MICPparaMicrobialDeathRate();
        Scalar rho_b = MICPparaDensityBiofilm();
        Scalar rho_c = MICPparaDensityCalcite();
        Scalar k_str = MICPparaDetachmentRate();
        Scalar k_o = MICPparaHalfVelocityOxygen();
        Scalar k_u = MICPparaHalfVelocityUrea() / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp)
        Scalar mu = MICPparaMaximumGrowthRate();
        Scalar mu_u = MICPparaMaximumUreaUtilization() / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp)
        Scalar Y_sb = MICPparaYieldGrowthCoefficient();
        Scalar F = MICPparaOxygenConsumptionFactor();
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

    static const Scalar MICPparaDensityBiofilm()
    {
        return MICPparaDensityBiofilm_;
    }

    static const Scalar MICPparaDensityCalcite()
    {
        return MICPparaDensityCalcite_;
    }

    static const Scalar MICPparaDetachmentRate()
    {
        return MICPparaDetachmentRate_;
    }

    static const Scalar MICPparaCriticalPorosity()
    {
        return MICPparaCriticalPorosity_;
    }

    static const Scalar MICPparaFittingFactor()
    {
        return MICPparaFittingFactor_;
    }

    static const Scalar MICPparaHalfVelocityOxygen()
    {
        return MICPparaHalfVelocityOxygen_;
    }

    static const Scalar MICPparaHalfVelocityUrea()
    {
        return MICPparaHalfVelocityUrea_;
    }

    static const Scalar MICPparaMaximumGrowthRate()
    {
        return MICPparaMaximumGrowthRate_;
    }

    static const Scalar MICPparaMaximumOxygenConcentration()
    {
        return MICPparaMaximumOxygenConcentration_;
    }

    static const Scalar MICPparaMaximumUreaConcentration()
    {
        return MICPparaMaximumUreaConcentration_ / 10.0;//Dividing by scaling factor 10 (see WellInterface_impl.hpp);
    }

    static const Scalar MICPparaMaximumUreaUtilization()
    {
        return MICPparaMaximumUreaUtilization_;
    }

    static const Scalar MICPparaMicrobialAttachmentRate()
    {
        return MICPparaMicrobialAttachmentRate_;
    }

    static const Scalar MICPparaMicrobialDeathRate()
    {
        return MICPparaMicrobialDeathRate_;
    }

    static const Scalar MICPparaMinimumPermeability()
    {
        return MICPparaMinimumPermeability_;
    }

    static const Scalar MICPparaOxygenConsumptionFactor()
    {
        return MICPparaOxygenConsumptionFactor_;
    }

    static const Scalar MICPparaToleranceBeforeClogging()
    {
        return MICPparaToleranceBeforeClogging_;
    }

    static const Scalar MICPparaYieldGrowthCoefficient()
    {
        return MICPparaYieldGrowthCoefficient_;
    }

    static const std::vector<Scalar> phi()
    {
        return phi_;
    }

private:
    static Scalar MICPparaDensityBiofilm_;
    static Scalar MICPparaDensityCalcite_;
    static Scalar MICPparaDetachmentRate_;
    static Scalar MICPparaCriticalPorosity_;
    static Scalar MICPparaFittingFactor_;
    static Scalar MICPparaHalfVelocityOxygen_;
    static Scalar MICPparaHalfVelocityUrea_;
    static Scalar MICPparaMaximumGrowthRate_;
    static Scalar MICPparaMaximumUreaUtilization_;
    static Scalar MICPparaMicrobialAttachmentRate_;
    static Scalar MICPparaMicrobialDeathRate_;
    static Scalar MICPparaMinimumPermeability_;
    static Scalar MICPparaOxygenConsumptionFactor_;
    static Scalar MICPparaYieldGrowthCoefficient_;
    static Scalar MICPparaMaximumOxygenConcentration_;
    static Scalar MICPparaMaximumUreaConcentration_;
    static Scalar MICPparaToleranceBeforeClogging_;
    static std::vector<Scalar> phi_;
};


template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaDensityBiofilm_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaDensityCalcite_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaDetachmentRate_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaCriticalPorosity_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaFittingFactor_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaHalfVelocityOxygen_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaHalfVelocityUrea_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMaximumGrowthRate_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMaximumUreaUtilization_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMicrobialAttachmentRate_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMicrobialDeathRate_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMinimumPermeability_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaOxygenConsumptionFactor_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaYieldGrowthCoefficient_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMaximumOxygenConcentration_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaMaximumUreaConcentration_;

template <class TypeTag, bool enableMICPV>
typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar
BlackOilMICPModule<TypeTag, enableMICPV>::MICPparaToleranceBeforeClogging_;

template <class TypeTag, bool enableMICPV>
std::vector<typename BlackOilMICPModule<TypeTag, enableMICPV>::Scalar>
BlackOilMICPModule<TypeTag, enableMICPV>::phi_;


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
        Scalar eta = MICPModule::MICPparaFittingFactor();
        Scalar k_min = MICPModule::MICPparaMinimumPermeability();
        Scalar phi_crit = MICPModule::MICPparaCriticalPorosity();

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
