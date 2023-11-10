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
 * \brief Contains the classes required to extend the black-oil model by solvents.
 */
#ifndef EWOMS_BLACK_OIL_SOLVENT_MODULE_HH
#define EWOMS_BLACK_OIL_SOLVENT_MODULE_HH

#include "blackoilproperties.hh"

#include <opm/common/Exceptions.hpp>

#include <opm/models/blackoil/blackoilsolventparams.hh>
#include <opm/models/io/vtkblackoilsolventmodule.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PmiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SorwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgcwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TlpmixpaTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by solvents.
 */
template <class TypeTag, bool enableSolventV = getPropValue<TypeTag, Properties::EnableSolvent>()>
class BlackOilSolventModule
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
    using SolventPvt = typename BlackOilSolventParams<Scalar>::SolventPvt;
    using Co2GasPvt = typename BlackOilSolventParams<Scalar>::Co2GasPvt;
    using H2GasPvt = typename BlackOilSolventParams<Scalar>::H2GasPvt;
    using BrineCo2Pvt = typename BlackOilSolventParams<Scalar>::BrineCo2Pvt;
    using BrineH2Pvt = typename BlackOilSolventParams<Scalar>::BrineH2Pvt;

    using TabulatedFunction = typename BlackOilSolventParams<Scalar>::TabulatedFunction;

    static constexpr unsigned solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr unsigned contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr unsigned enableSolvent = enableSolventV;
    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;
    static constexpr bool blackoilConserveSurfaceVolume = getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>();
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;


public:
#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the solvent module
     */
    static void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        // some sanity checks: if solvents are enabled, the SOLVENT keyword must be
        // present, if solvents are disabled the keyword must not be present.
        if (enableSolvent && !eclState.runspec().phases().active(Phase::SOLVENT))
            throw std::runtime_error("Non-trivial solvent treatment requested at compile "
                                     "time, but the deck does not contain the SOLVENT keyword");
        else if (!enableSolvent && eclState.runspec().phases().active(Phase::SOLVENT))
            throw std::runtime_error("Solvent treatment disabled at compile time, but the deck "
                                     "contains the SOLVENT keyword");

        if (!eclState.runspec().phases().active(Phase::SOLVENT))
            return; // solvent treatment is supposed to be disabled

        params_.co2sol_ = eclState.runspec().co2Sol();
        params_.h2sol_ = eclState.runspec().h2Sol();

        if (isCO2Sol() && isH2Sol()) {
             throw std::runtime_error("CO2SOL and H2SOL can not be used together");
        }

        if (isCO2Sol() || isH2Sol() ) {
            if (isCO2Sol()) {
                params_.co2GasPvt_.initFromState(eclState, schedule);
                params_.brineCo2Pvt_.initFromState(eclState, schedule);
            } else { 
                params_.h2GasPvt_.initFromState(eclState, schedule);
                params_.brineH2Pvt_.initFromState(eclState, schedule);
            }
            if (eclState.getSimulationConfig().hasDISGASW()) {
                params_.rsSolw_active_ = true;
            }
        } else 
            params_.solventPvt_.initFromState(eclState, schedule);

        const auto& tableManager = eclState.getTableManager();
        // initialize the objects which deal with the SSFN keyword
        const auto& ssfnTables = tableManager.getSsfnTables();
        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        params_.setNumSatRegions(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
            const auto& ssfnTable = ssfnTables.template getTable<SsfnTable>(satRegionIdx);
            params_.ssfnKrg_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                                           ssfnTable.getGasRelPermMultiplierColumn(),
                                                           /*sortInput=*/true);
            params_.ssfnKrs_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                                           ssfnTable.getSolventRelPermMultiplierColumn(),
                                                           /*sortInput=*/true);
        }

        // initialize the objects needed for miscible solvent and oil simulations
        params_.isMiscible_ = false;
        if (!eclState.getTableManager().getMiscTables().empty()) {
            params_.isMiscible_ = true;

            unsigned numMiscRegions = 1;

            // misicible hydrocabon relative permeability wrt water
            const auto& sof2Tables = tableManager.getSof2Tables();
            if (!sof2Tables.empty()) {
                // resize the attributes of the object
                params_.sof2Krn_.resize(numSatRegions);
                for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                    const auto& sof2Table = sof2Tables.template getTable<Sof2Table>(satRegionIdx);
                    params_.sof2Krn_[satRegionIdx].setXYContainers(sof2Table.getSoColumn(),
                                                                   sof2Table.getKroColumn(),
                                                                   /*sortInput=*/true);
                }
            }
            else
                throw std::runtime_error("SOF2 must be specified in MISCIBLE (SOLVENT) runs\n");

            const auto& miscTables = tableManager.getMiscTables();
            if (!miscTables.empty()) {
                assert(numMiscRegions == miscTables.size());

                // resize the attributes of the object
                params_.misc_.resize(numMiscRegions);
                for (unsigned miscRegionIdx = 0; miscRegionIdx < numMiscRegions; ++miscRegionIdx) {
                    const auto& miscTable = miscTables.template getTable<MiscTable>(miscRegionIdx);

                    // solventFraction = Ss / (Ss + Sg);
                    const auto& solventFraction = miscTable.getSolventFractionColumn();
                    const auto& misc = miscTable.getMiscibilityColumn();
                    params_.misc_[miscRegionIdx].setXYContainers(solventFraction, misc);
                }
            }
            else
                throw std::runtime_error("MISC must be specified in MISCIBLE (SOLVENT) runs\n");

            // resize the attributes of the object
            params_.pmisc_.resize(numMiscRegions);
            const auto& pmiscTables = tableManager.getPmiscTables();
            if (!pmiscTables.empty()) {
                assert(numMiscRegions == pmiscTables.size());

                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    const auto& pmiscTable = pmiscTables.template getTable<PmiscTable>(regionIdx);

                    // Copy data
                    const auto& po = pmiscTable.getOilPhasePressureColumn();
                    const auto& pmisc = pmiscTable.getMiscibilityColumn();

                    params_.pmisc_[regionIdx].setXYContainers(po, pmisc);
                }
            }
            else {
                std::vector<double> x = {0.0,1.0e20};
                std::vector<double> y = {1.0,1.0};
                TabulatedFunction constant = TabulatedFunction(2, x, y);
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    params_.pmisc_[regionIdx] = constant;
                }
            }

            // miscible relative permeability multipleiers
            params_.msfnKrsg_.resize(numSatRegions);
            params_.msfnKro_.resize(numSatRegions);
            const auto& msfnTables = tableManager.getMsfnTables();
            if (!msfnTables.empty()) {
                assert(numSatRegions == msfnTables.size());

                for (unsigned regionIdx = 0; regionIdx < numSatRegions; ++regionIdx) {
                    const MsfnTable& msfnTable = msfnTables.template getTable<MsfnTable>(regionIdx);

                    // Copy data
                    // Ssg = Ss + Sg;
                    const auto& Ssg = msfnTable.getGasPhaseFractionColumn();
                    const auto& krsg = msfnTable.getGasSolventRelpermMultiplierColumn();
                    const auto& kro = msfnTable.getOilRelpermMultiplierColumn();

                    params_.msfnKrsg_[regionIdx].setXYContainers(Ssg, krsg);
                    params_.msfnKro_[regionIdx].setXYContainers(Ssg, kro);
                }
            }
            else {
                std::vector<double> x = {0.0,1.0};
                std::vector<double> y = {1.0,0.0};
                TabulatedFunction unit = TabulatedFunction(2, x, x);
                TabulatedFunction invUnit = TabulatedFunction(2, x, y);

                for (unsigned regionIdx = 0; regionIdx < numSatRegions; ++regionIdx) {
                    params_.setMsfn(regionIdx, unit, invUnit);
                }
            }
            // resize the attributes of the object
            params_.sorwmis_.resize(numMiscRegions);
            const auto& sorwmisTables = tableManager.getSorwmisTables();
            if (!sorwmisTables.empty()) {
                assert(numMiscRegions == sorwmisTables.size());

                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    const auto& sorwmisTable = sorwmisTables.template getTable<SorwmisTable>(regionIdx);

                    // Copy data
                    const auto& sw = sorwmisTable.getWaterSaturationColumn();
                    const auto& sorwmis = sorwmisTable.getMiscibleResidualOilColumn();

                    params_.sorwmis_[regionIdx].setXYContainers(sw, sorwmis);
                }
            }
            else {
                // default
                std::vector<double> x = {0.0,1.0};
                std::vector<double> y = {0.0,0.0};
                TabulatedFunction zero = TabulatedFunction(2, x, y);
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    params_.sorwmis_[regionIdx] = zero;
                }
            }

            // resize the attributes of the object
            params_.sgcwmis_.resize(numMiscRegions);
            const auto& sgcwmisTables = tableManager.getSgcwmisTables();
            if (!sgcwmisTables.empty()) {
                assert(numMiscRegions == sgcwmisTables.size());

                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    const auto& sgcwmisTable = sgcwmisTables.template getTable<SgcwmisTable>(regionIdx);

                    // Copy data
                    const auto& sw = sgcwmisTable.getWaterSaturationColumn();
                    const auto& sgcwmis = sgcwmisTable.getMiscibleResidualGasColumn();

                    params_.sgcwmis_[regionIdx].setXYContainers(sw, sgcwmis);
                }
            }
            else {
                // default
                std::vector<double> x = {0.0,1.0};
                std::vector<double> y = {0.0,0.0};
                TabulatedFunction zero = TabulatedFunction(2, x, y);
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx)
                    params_.sgcwmis_[regionIdx] = zero;
            }

            const auto& tlmixpar = eclState.getTableManager().getTLMixpar();
            if (!tlmixpar.empty()) {
                // resize the attributes of the object
                params_.tlMixParamViscosity_.resize(numMiscRegions);
                params_.tlMixParamDensity_.resize(numMiscRegions);

                assert(numMiscRegions == tlmixpar.size());
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    const auto& tlp = tlmixpar[regionIdx];
                    params_.tlMixParamViscosity_[regionIdx] = tlp.viscosity_parameter;
                    params_.tlMixParamDensity_[regionIdx] = tlp.density_parameter;
                }
            }
            else
                throw std::runtime_error("TLMIXPAR must be specified in MISCIBLE (SOLVENT) runs\n");

            // resize the attributes of the object
            params_.tlPMixTable_.resize(numMiscRegions);
            if (!eclState.getTableManager().getTlpmixpaTables().empty()) {
                const auto& tlpmixparTables = tableManager.getTlpmixpaTables();
                if (!tlpmixparTables.empty()) {
                    assert(numMiscRegions == tlpmixparTables.size());
                    for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                        const auto& tlpmixparTable = tlpmixparTables.template getTable<TlpmixpaTable>(regionIdx);

                        // Copy data
                        const auto& po = tlpmixparTable.getOilPhasePressureColumn();
                        const auto& tlpmixpa = tlpmixparTable.getMiscibilityColumn();

                        params_.tlPMixTable_[regionIdx].setXYContainers(po, tlpmixpa);
                    }
                }
                else {
                    // if empty keyword. Try to use the pmisc table as default.
                    if (params_.pmisc_.size() > 0)
                        params_.tlPMixTable_ = params_.pmisc_;
                    else
                        throw std::invalid_argument("If the pressure dependent TL values in "
                                                    "TLPMIXPA is defaulted (no entries), then "
                                                    "the PMISC tables must be specified.");
                }
            }
            else {
                // default
                std::vector<double> x = {0.0,1.0e20};
                std::vector<double> y = {1.0,1.0};
                TabulatedFunction ones = TabulatedFunction(2, x, y);
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx)
                    params_.tlPMixTable_[regionIdx] = ones;
            }
        }
    }
#endif

    /*!
     * \brief Specify the solvent PVT of a all PVT regions.
     */
    static void setSolventPvt(const SolventPvt& value)
    { params_.solventPvt_ = value; }


    static void setIsMiscible(const bool isMiscible)
    { params_.isMiscible_ = isMiscible; }

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        if constexpr (enableSolvent)
            VtkBlackOilSolventModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if constexpr (enableSolvent)
            model.addOutputModule(new VtkBlackOilSolventModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if constexpr (enableSolvent)
            return pvIdx == solventSaturationIdx;
        else
            return false;
    }

    static std::string primaryVarName([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        return "saturation_solvent";
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableSolvent)
            return eqIdx == contiSolventEqIdx;
        else
            return false;
    }

    static std::string eqName([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if constexpr (enableSolvent) {
            if constexpr (blackoilConserveSurfaceVolume) {
                storage[contiSolventEqIdx] +=
                        Toolbox::template decay<LhsEval>(intQuants.porosity())
                        * Toolbox::template decay<LhsEval>(intQuants.solventSaturation())
                        * Toolbox::template decay<LhsEval>(intQuants.solventInverseFormationVolumeFactor());
                if (isSolubleInWater()) {
                    storage[contiSolventEqIdx] += Toolbox::template decay<LhsEval>(intQuants.porosity())
                        * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(waterPhaseIdx))
                        * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(waterPhaseIdx)) 
                        * Toolbox::template decay<LhsEval>(intQuants.rsSolw());
                }
            }
            else {
                storage[contiSolventEqIdx] +=
                        Toolbox::template decay<LhsEval>(intQuants.porosity())
                        * Toolbox::template decay<LhsEval>(intQuants.solventSaturation())
                        * Toolbox::template decay<LhsEval>(intQuants.solventDensity());
                if (isSolubleInWater()) {
                    storage[contiSolventEqIdx] += Toolbox::template decay<LhsEval>(intQuants.porosity())
                        * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(waterPhaseIdx))
                        * Toolbox::template decay<LhsEval>(intQuants.fluidState().density(waterPhaseIdx))
                        * Toolbox::template decay<LhsEval>(intQuants.rsSolw());
                }
                         
            }
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)

    {
        if constexpr (enableSolvent) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

            unsigned upIdx = extQuants.solventUpstreamIndex();
            unsigned inIdx = extQuants.interiorIndex();
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            if constexpr (blackoilConserveSurfaceVolume) {
                if (upIdx == inIdx) 
                    flux[contiSolventEqIdx] =
                            extQuants.solventVolumeFlux()
                            *up.solventInverseFormationVolumeFactor();
                else
                    flux[contiSolventEqIdx] =
                            extQuants.solventVolumeFlux()
                            *decay<Scalar>(up.solventInverseFormationVolumeFactor());


                if (isSolubleInWater()) {
                    if (upIdx == inIdx) 
                        flux[contiSolventEqIdx] =
                                extQuants.volumeFlux(waterPhaseIdx)
                                * up.fluidState().invB(waterPhaseIdx)
                                * up.rsSolw();
                    else
                        flux[contiSolventEqIdx] =
                                extQuants.volumeFlux(waterPhaseIdx)
                                *decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                                *decay<Scalar>(up.rsSolw());
                }
            }
            else {
                if (upIdx == inIdx)
                    flux[contiSolventEqIdx] =
                            extQuants.solventVolumeFlux()
                            *up.solventDensity();
                else
                    flux[contiSolventEqIdx] =
                            extQuants.solventVolumeFlux()
                            *decay<Scalar>(up.solventDensity());


                if (isSolubleInWater()) {
                    if (upIdx == inIdx) 
                        flux[contiSolventEqIdx] =
                                extQuants.volumeFlux(waterPhaseIdx)
                                * up.fluidState().density(waterPhaseIdx)
                                * up.rsSolw();
                    else
                        flux[contiSolventEqIdx] =
                                extQuants.volumeFlux(waterPhaseIdx)
                                *decay<Scalar>(up.fluidState().density(waterPhaseIdx))
                                *decay<Scalar>(up.rsSolw());
                }
            }
        }
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar solventSaturation,
                                  Scalar solventRsw)
    {
        if constexpr (!enableSolvent) {
            priVars.setPrimaryVarsMeaningSolvent(PrimaryVariables::SolventMeaning::Disabled);
            return;
        }
        // Determine the meaning of the solvent primary variables
        if (solventSaturation > 0 || !isSolubleInWater()) {
            priVars.setPrimaryVarsMeaningSolvent(PrimaryVariables::SolventMeaning::Ss);
            priVars[solventSaturationIdx] = solventSaturation;
        } else {
            priVars.setPrimaryVarsMeaningSolvent(PrimaryVariables::SolventMeaning::Rsolw);
            priVars[solventSaturationIdx] = solventRsw;
        } 
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if constexpr (enableSolvent)
            // do a plain unchopped Newton update
            newPv[solventSaturationIdx] = oldPv[solventSaturationIdx] - delta[solventSaturationIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables&,
                                     const EqVector&)
    {
        // do not consider consider the cange of solvent primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of solvents when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiSolventEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if constexpr (enableSolvent) {
            unsigned dofIdx = model.dofMapper().index(dof);

            const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
            outstream << priVars[solventSaturationIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if constexpr (enableSolvent) {
            unsigned dofIdx = model.dofMapper().index(dof);

            PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
            PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

            instream >> priVars0[solventSaturationIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1 = priVars0[solventSaturationIdx];
        }
    }

    static const SolventPvt& solventPvt()
    { 
        return params_.solventPvt_;    
    }


    static const Co2GasPvt& co2GasPvt()
    { 
            return params_.co2GasPvt_;
    }

    static const H2GasPvt& h2GasPvt()
    { 
            return params_.h2GasPvt_;
    }

    static const BrineCo2Pvt& brineCo2Pvt()
    { 
            return params_.brineCo2Pvt_;
    }

    static const BrineH2Pvt& brineH2Pvt()
    { 
            return params_.brineH2Pvt_;
    }

    static const TabulatedFunction& ssfnKrg(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.ssfnKrg_[satnumRegionIdx];
    }

    static const TabulatedFunction& ssfnKrs(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.ssfnKrs_[satnumRegionIdx];
    }

    static const TabulatedFunction& sof2Krn(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.sof2Krn_[satnumRegionIdx];
    }

    static const TabulatedFunction& misc(const ElementContext& elemCtx,
                                         unsigned scvIdx,
                                         unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.misc_[miscnumRegionIdx];
    }

    static const TabulatedFunction& pmisc(const ElementContext& elemCtx,
                                          unsigned scvIdx,
                                          unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.pmisc_[miscnumRegionIdx];
    }

    static const TabulatedFunction& msfnKrsg(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.msfnKrsg_[satnumRegionIdx];
    }

    static const TabulatedFunction& msfnKro(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.msfnKro_[satnumRegionIdx];
    }

    static const TabulatedFunction& sorwmis(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.sorwmis_[miscnumRegionIdx];
    }

    static const TabulatedFunction& sgcwmis(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.sgcwmis_[miscnumRegionIdx];
    }

    static const TabulatedFunction& tlPMixTable(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.tlPMixTable_[miscnumRegionIdx];
    }

    static const Scalar& tlMixParamViscosity(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.tlMixParamViscosity_[miscnumRegionIdx];
    }

    static const Scalar& tlMixParamDensity(const ElementContext& elemCtx,
                                           unsigned scvIdx,
                                           unsigned timeIdx)
    {
        unsigned miscnumRegionIdx = elemCtx.problem().miscnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.tlMixParamDensity_[miscnumRegionIdx];
    }

    static bool isMiscible()
    {
        return params_.isMiscible_;
    }

    template <class Value>
    static const Value solubilityLimit(unsigned pvtIdx, const Value& temperature, const Value& pressure, const Value& saltConcentration) 
    {
        if (!isSolubleInWater())
            return 0.0;

        assert(isCO2Sol() || isH2Sol());
        if (isCO2Sol())
            return brineCo2Pvt().saturatedGasDissolutionFactor(pvtIdx, temperature, pressure, saltConcentration);
        else
            return brineH2Pvt().saturatedGasDissolutionFactor(pvtIdx, temperature, pressure, saltConcentration);
    }

    static bool isSolubleInWater() 
    {
        return params_.rsSolw_active_;
    } 

    static bool isCO2Sol() 
    {
        return params_.co2sol_;
    }

    static bool isH2Sol() 
    {
        return params_.h2sol_;
    }  

private:
    static BlackOilSolventParams<Scalar> params_; // the krg(Fs) column of the SSFN table
};

template <class TypeTag, bool enableSolventV>
BlackOilSolventParams<typename BlackOilSolventModule<TypeTag, enableSolventV>::Scalar>
BlackOilSolventModule<TypeTag, enableSolventV>::params_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilSolventIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag, bool enableSolventV = getPropValue<TypeTag, Properties::EnableSolvent>()>
class BlackOilSolventIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using SolventModule = BlackOilSolventModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr double cutOff = 1e-12;


public:
    /*!
     * \brief Called before the saturation functions are doing their magic
     *
     * At this point, the saturations of the fluid state correspond to those if the phases
     * were pure hydrocarbons.
     */
    void solventPreSatFuncUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        auto& fs = asImp_().fluidState_;
        solventSaturation_ = 0.0;
        if (priVars.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
            solventSaturation_ = priVars.makeEvaluation(solventSaturationIdx, timeIdx, elemCtx.linearizationType());
        }

        hydrocarbonSaturation_ = fs.saturation(gasPhaseIdx);

        // apply a cut-off. Don't waste calculations if no solvent
        if (solventSaturation().value() < cutOff)
            return;

        // make the saturation of the gas phase which is used by the saturation functions
        // the sum of the solvent "saturation" and the saturation the hydrocarbon gas.
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_ + solventSaturation_);
    }

    /*!
     * \brief Called after the saturation functions have been doing their magic
     *
     * After this function, all saturations, pressures
     * and relative permeabilities must be final. (i.e., the "hydrocarbon
     * saturations".)
     */
    void solventPostSatFuncUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {
        // revert the gas "saturation" of the fluid state back to the saturation of the
        // hydrocarbon gas.
        auto& fs = asImp_().fluidState_;
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_);


        // update rsSolw. This needs to be done after the pressure is defined in the fluid state. 
        rsSolw_ = 0.0;
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        if (priVars.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
            rsSolw_ = SolventModule::solubilityLimit(asImp_().pvtRegionIndex(), fs.temperature(waterPhaseIdx), fs.pressure(waterPhaseIdx), fs.saltConcentration());
        } else if (priVars.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Rsolw) {
            rsSolw_ = priVars.makeEvaluation(solventSaturationIdx, timeIdx, elemCtx.linearizationType());
        }

        solventMobility_ = 0.0;

        // apply a cut-off. Don't waste calculations if no solvent
        if (solventSaturation().value() < cutOff)
            return;

        // Pressure effects on capillary pressure miscibility
        if (SolventModule::isMiscible()) {
            const Evaluation& p = fs.pressure(oilPhaseIdx); // or gas pressure?
            const Evaluation pmisc = SolventModule::pmisc(elemCtx, dofIdx, timeIdx).eval(p, /*extrapolate=*/true);
            const Evaluation& pgImisc = fs.pressure(gasPhaseIdx);

            // compute capillary pressure for miscible fluid
            const auto& problem = elemCtx.problem();
            Evaluation pgMisc = 0.0;
            std::array<Evaluation, numPhases> pC;
            const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, materialParams, fs);

            //oil is the reference phase for pressure
            const auto linearizationType = elemCtx.linearizationType();
            if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg)
                pgMisc = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx, linearizationType);
            else {
                const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx, linearizationType);
                pgMisc = po + (pC[gasPhaseIdx] - pC[oilPhaseIdx]);
            }

           fs.setPressure(gasPhaseIdx, pmisc * pgMisc + (1.0 - pmisc) * pgImisc);
        }


        Evaluation gasSolventSat = hydrocarbonSaturation_ + solventSaturation_;

        if (gasSolventSat.value() < cutOff) // avoid division by zero
            return;

        Evaluation Fhydgas = hydrocarbonSaturation_/gasSolventSat;
        Evaluation Fsolgas = solventSaturation_/gasSolventSat;

        // account for miscibility of oil and solvent
        if (SolventModule::isMiscible()) {
            const auto& misc = SolventModule::misc(elemCtx, dofIdx, timeIdx);
            const auto& pmisc = SolventModule::pmisc(elemCtx, dofIdx, timeIdx);
            const Evaluation& p = fs.pressure(oilPhaseIdx); // or gas pressure?
            const Evaluation miscibility = misc.eval(Fsolgas, /*extrapolate=*/true) * pmisc.eval(p, /*extrapolate=*/true);

            // TODO adjust endpoints of sn and ssg
            unsigned cellIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
            const auto& materialLawManager = elemCtx.problem().materialLawManager();
            const auto& scaledDrainageInfo =
                    materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);

            const Scalar& sgcr = scaledDrainageInfo.Sgcr;
            const Scalar& sogcr = scaledDrainageInfo.Sogcr;
            const Evaluation& sw = fs.saturation(waterPhaseIdx);
            const auto& sorwmis = SolventModule::sorwmis(elemCtx, dofIdx, timeIdx);
            const auto& sgcwmis = SolventModule::sgcwmis(elemCtx, dofIdx, timeIdx);

            Evaluation sor = miscibility * sorwmis.eval(sw,  /*extrapolate=*/true) + (1.0 - miscibility) * sogcr;
            Evaluation sgc = miscibility * sgcwmis.eval(sw,  /*extrapolate=*/true) + (1.0 - miscibility) * sgcr;

            const Evaluation oilGasSolventSat = gasSolventSat + fs.saturation(oilPhaseIdx);
            const Evaluation zero = 0.0;
            const Evaluation oilGasSolventEffSat = std::max(oilGasSolventSat - sor - sgc, zero);

            Evaluation F_totalGas = 0.0;
            if (oilGasSolventEffSat.value() > cutOff) {
                const Evaluation gasSolventEffSat = std::max(gasSolventSat - sgc, zero);
                F_totalGas = gasSolventEffSat / oilGasSolventEffSat;
            }
            const auto& msfnKro = SolventModule::msfnKro(elemCtx, dofIdx, timeIdx);
            const auto& msfnKrsg = SolventModule::msfnKrsg(elemCtx, dofIdx, timeIdx);
            const auto& sof2Krn = SolventModule::sof2Krn(elemCtx, dofIdx, timeIdx);

            const Evaluation mkrgt = msfnKrsg.eval(F_totalGas, /*extrapolate=*/true) * sof2Krn.eval(oilGasSolventSat, /*extrapolate=*/true);
            const Evaluation mkro = msfnKro.eval(F_totalGas, /*extrapolate=*/true) * sof2Krn.eval(oilGasSolventSat, /*extrapolate=*/true);

            Evaluation& kro = asImp_().mobility_[oilPhaseIdx];
            Evaluation& krg = asImp_().mobility_[gasPhaseIdx];

            // combine immiscible and miscible part of the relperm
            krg *= (1.0 - miscibility);
            krg += miscibility * mkrgt;
            kro *= (1.0 - miscibility);
            kro += miscibility * mkro;
        }



        // compute the mobility of the solvent "phase" and modify the gas phase
        const auto& ssfnKrg = SolventModule::ssfnKrg(elemCtx, dofIdx, timeIdx);
        const auto& ssfnKrs = SolventModule::ssfnKrs(elemCtx, dofIdx, timeIdx);

        Evaluation& krg = asImp_().mobility_[gasPhaseIdx];
        solventMobility_ = krg * ssfnKrs.eval(Fsolgas, /*extrapolate=*/true);
        krg *= ssfnKrg.eval(Fhydgas, /*extrapolate=*/true);

    }

    /*!
     * \brief Update the intensive PVT properties needed to handle solvents from the
     *        primary variables.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void solventPvtUpdate_(const ElementContext& elemCtx,
                           unsigned scvIdx,
                           unsigned timeIdx)
    {
        const auto& iq = asImp_();
        unsigned pvtRegionIdx = iq.pvtRegionIndex();
        const Evaluation& T = iq.fluidState().temperature(gasPhaseIdx);
        const Evaluation& p = iq.fluidState().pressure(gasPhaseIdx);

        const Evaluation rv = 0.0;
        const Evaluation rvw = 0.0;
        if (SolventModule::isCO2Sol() || SolventModule::isH2Sol() ){
            if (SolventModule::isCO2Sol()) {
                const auto& co2gasPvt = SolventModule::co2GasPvt();
                solventInvFormationVolumeFactor_ = co2gasPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, rv, rvw);
                solventRefDensity_ = co2gasPvt.gasReferenceDensity(pvtRegionIdx);
                solventViscosity_ = co2gasPvt.viscosity(pvtRegionIdx, T, p, rv, rvw);

                const auto& brineCo2Pvt = SolventModule::brineCo2Pvt();
                auto& fs = asImp_().fluidState_;
                const auto& bw = brineCo2Pvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, rsSolw());

                const auto denw = bw*brineCo2Pvt.waterReferenceDensity(pvtRegionIdx)
                        + rsSolw()*bw*brineCo2Pvt.gasReferenceDensity(pvtRegionIdx);
                fs.setDensity(waterPhaseIdx, denw);
                fs.setInvB(waterPhaseIdx, bw);
                Evaluation& mobw = asImp_().mobility_[waterPhaseIdx];
                const auto& muWat = fs.viscosity(waterPhaseIdx);
                const auto& muWatEff = brineCo2Pvt.viscosity(pvtRegionIdx, T, p, rsSolw());
                mobw *= muWat / muWatEff;
            } else { 
                const auto& h2gasPvt = SolventModule::h2GasPvt();
                solventInvFormationVolumeFactor_ = h2gasPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, rv, rvw);
                solventRefDensity_ = h2gasPvt.gasReferenceDensity(pvtRegionIdx);
                solventViscosity_ = h2gasPvt.viscosity(pvtRegionIdx, T, p, rv, rvw); 
            
                const auto& brineH2Pvt = SolventModule::brineH2Pvt();
                auto& fs = asImp_().fluidState_;
                const auto& bw = brineH2Pvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, rsSolw());

                const auto denw = bw*brineH2Pvt.waterReferenceDensity(pvtRegionIdx)
                        + rsSolw()*bw*brineH2Pvt.gasReferenceDensity(pvtRegionIdx);
                fs.setDensity(waterPhaseIdx, denw);
                fs.setInvB(waterPhaseIdx, bw);
                Evaluation& mobw = asImp_().mobility_[waterPhaseIdx];
                const auto& muWat = fs.viscosity(waterPhaseIdx);
                const auto& muWatEff = brineH2Pvt.viscosity(pvtRegionIdx, T, p, rsSolw());
                mobw *= muWat / muWatEff;
            }
        } else {
            const auto& solventPvt = SolventModule::solventPvt();
            solventInvFormationVolumeFactor_ = solventPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p);
            solventRefDensity_ = solventPvt.referenceDensity(pvtRegionIdx);
            solventViscosity_ = solventPvt.viscosity(pvtRegionIdx, T, p);
        }   

        solventDensity_ = solventInvFormationVolumeFactor_*solventRefDensity_;
        effectiveProperties(elemCtx, scvIdx, timeIdx);

        solventMobility_ /= solventViscosity_;


    }

    const Evaluation& solventSaturation() const
    { return solventSaturation_; }

    const Evaluation& rsSolw() const
    { return rsSolw_; }

    const Evaluation& solventDensity() const
    { return solventDensity_; }

    const Evaluation& solventViscosity() const
    { return solventViscosity_; }

    const Evaluation& solventMobility() const
    { return solventMobility_; }

    const Evaluation& solventInverseFormationVolumeFactor() const
    { return solventInvFormationVolumeFactor_; }

    // This could be stored pr pvtRegion instead
    const Scalar& solventRefDensity() const
    { return solventRefDensity_; }

private:
    // Computes the effective properties based on
    // Todd-Longstaff mixing model.
    void effectiveProperties(const ElementContext& elemCtx,
                             unsigned scvIdx,
                             unsigned timeIdx)
    {
        if (!SolventModule::isMiscible())
            return;

        // Don't waste calculations if no solvent
        // Apply a cut-off for small and negative solvent saturations
        if (solventSaturation() < cutOff)
            return;

        auto& fs = asImp_().fluidState_;

        // Compute effective saturations
        const auto& sorwmis = SolventModule::sorwmis(elemCtx, scvIdx, timeIdx);
        const auto& sgcwmis = SolventModule::sgcwmis(elemCtx, scvIdx, timeIdx);
        const Evaluation& sw = fs.saturation(waterPhaseIdx);

        const Evaluation zero = 0.0;
        const Evaluation oilEffSat = std::max(fs.saturation(oilPhaseIdx) - sorwmis.eval(sw,  /*extrapolate=*/true),zero);
        const Evaluation gasEffSat = std::max(fs.saturation(gasPhaseIdx) - sgcwmis.eval(sw,  /*extrapolate=*/true),zero);
        const Evaluation solventEffSat = std::max(solventSaturation() - sgcwmis.eval(sw,  /*extrapolate=*/true),zero);

        const Evaluation oilGasSolventEffSat =  oilEffSat + gasEffSat + solventEffSat;
        const Evaluation oilSolventEffSat = oilEffSat + solventEffSat;
        const Evaluation solventGasEffSat = solventEffSat + gasEffSat;

        // Compute effective viscosities
        const Evaluation& muGas = fs.viscosity(gasPhaseIdx);
        const Evaluation& muOil = fs.viscosity(oilPhaseIdx);
        const Evaluation& muSolvent = solventViscosity_;

        assert(muOil.value() > 0);
        assert(muGas.value() > 0);
        assert(muSolvent.value() > 0);
        const Evaluation muOilPow = pow(muOil, 0.25);
        const Evaluation muGasPow = pow(muGas, 0.25);
        const Evaluation muSolventPow = pow(muSolvent, 0.25);

        Evaluation muMixOilSolvent = muOil;
        if (oilSolventEffSat > cutOff)
            muMixOilSolvent *= muSolvent / pow(((oilEffSat / oilSolventEffSat) * muSolventPow) + ((solventEffSat / oilSolventEffSat) * muOilPow) , 4.0);

        Evaluation muMixSolventGas = muGas;
        if (solventGasEffSat > cutOff)
            muMixSolventGas *= muSolvent / pow(((gasEffSat / solventGasEffSat) * muSolventPow) + ((solventEffSat / solventGasEffSat) * muGasPow) , 4.0);

        Evaluation muMixSolventGasOil = muOil;
        if (oilGasSolventEffSat > cutOff)
            muMixSolventGasOil *= muSolvent * muGas / pow(((oilEffSat / oilGasSolventEffSat) * muSolventPow *  muGasPow)
                  + ((solventEffSat / oilGasSolventEffSat) * muOilPow *  muGasPow) + ((gasEffSat / oilGasSolventEffSat) * muSolventPow * muOilPow), 4.0);

        // Mixing parameter for viscosity
        // The pressureMixingParameter represent the miscibility of the solvent while the mixingParameterViscosity the effect of the porous media.
        // The pressureMixingParameter is not implemented in ecl100.
        const Evaluation& po = fs.pressure(oilPhaseIdx);
        const auto& tlPMixTable = SolventModule::tlPMixTable(elemCtx, scvIdx, timeIdx);
        const Evaluation tlMixParamMu = SolventModule::tlMixParamViscosity(elemCtx, scvIdx, timeIdx) * tlPMixTable.eval(po,  /*extrapolate=*/true);

        Evaluation muOilEff = pow(muOil,1.0 - tlMixParamMu) * pow(muMixOilSolvent, tlMixParamMu);
        Evaluation muGasEff = pow(muGas,1.0 - tlMixParamMu) * pow(muMixSolventGas, tlMixParamMu);
        Evaluation muSolventEff = pow(muSolvent,1.0 - tlMixParamMu) * pow(muMixSolventGasOil, tlMixParamMu);

        // Compute effective densities
        const Evaluation& rhoGas = fs.density(gasPhaseIdx);
        const Evaluation& rhoOil = fs.density(oilPhaseIdx);
        const Evaluation& rhoSolvent = solventDensity_;

        // Mixing parameter for density
        // The pressureMixingParameter represent the miscibility of the solvent while the mixingParameterDenisty the effect of the porous media.
        // The pressureMixingParameter is not implemented in ecl100.
        const Evaluation tlMixParamRho = SolventModule::tlMixParamDensity(elemCtx, scvIdx, timeIdx) * tlPMixTable.eval(po,  /*extrapolate=*/true);

        // compute effective viscosities for density calculations. These have to
        // be recomputed as a different mixing parameter may be used.
        const Evaluation muOilEffPow = pow(pow(muOil, 1.0 - tlMixParamRho) * pow(muMixOilSolvent, tlMixParamRho), 0.25);
        const Evaluation muGasEffPow = pow(pow(muGas, 1.0 - tlMixParamRho) * pow(muMixSolventGas, tlMixParamRho), 0.25);
        const Evaluation muSolventEffPow = pow(pow(muSolvent, 1.0 - tlMixParamRho) * pow(muMixSolventGasOil, tlMixParamRho), 0.25);

        const Evaluation oilGasEffSaturation = oilEffSat + gasEffSat;
        Evaluation sof = 0.0;
        Evaluation sgf = 0.0;
        if (oilGasEffSaturation.value() > cutOff) {
            sof = oilEffSat / oilGasEffSaturation;
            sgf = gasEffSat / oilGasEffSaturation;
        }

        const Evaluation muSolventOilGasPow = muSolventPow * ((sgf * muOilPow) + (sof * muGasPow));

        Evaluation rhoMixSolventGasOil = 0.0;
        if (oilGasSolventEffSat.value() > cutOff)
            rhoMixSolventGasOil = (rhoOil * oilEffSat / oilGasSolventEffSat) + (rhoGas * gasEffSat / oilGasSolventEffSat) + (rhoSolvent * solventEffSat / oilGasSolventEffSat);

        Evaluation rhoGasEff = 0.0;
        if (std::abs(muSolventPow.value() - muGasPow.value()) < cutOff)
            rhoGasEff = ((1.0 - tlMixParamRho) * rhoGas) + (tlMixParamRho * rhoMixSolventGasOil);
        else {
            const Evaluation solventGasEffFraction = (muGasPow * (muSolventPow - muGasEffPow)) / (muGasEffPow * (muSolventPow - muGasPow));
            rhoGasEff = (rhoGas * solventGasEffFraction) + (rhoSolvent * (1.0 - solventGasEffFraction));
        }

        Evaluation rhoOilEff = 0.0;
        if (std::abs(muOilPow.value() - muSolventPow.value()) < cutOff) {
            rhoOilEff = ((1.0 - tlMixParamRho) * rhoOil) + (tlMixParamRho * rhoMixSolventGasOil);
        }
        else {
            const Evaluation solventOilEffFraction = (muOilPow * (muOilEffPow - muSolventPow)) / (muOilEffPow * (muOilPow - muSolventPow));
            rhoOilEff = (rhoOil * solventOilEffFraction) + (rhoSolvent * (1.0 - solventOilEffFraction));
        }

        Evaluation rhoSolventEff = 0.0;
        if (std::abs((muSolventOilGasPow.value() - (muOilPow.value() * muGasPow.value()))) < cutOff)
            rhoSolventEff = ((1.0 - tlMixParamRho) * rhoSolvent) + (tlMixParamRho * rhoMixSolventGasOil);
        else {
            const Evaluation sfraction_se = (muSolventOilGasPow - (muOilPow * muGasPow * muSolventPow / muSolventEffPow)) / (muSolventOilGasPow - (muOilPow * muGasPow));
            rhoSolventEff = (rhoSolvent * sfraction_se) + (rhoGas * sgf * (1.0 - sfraction_se)) + (rhoOil * sof * (1.0 - sfraction_se));
        }

        unsigned pvtRegionIdx = asImp_().pvtRegionIndex();
        // compute invB from densities.
        const Evaluation bOilEff = rhoOilEff / (FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx) + FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx) * fs.Rs());
        const Evaluation bGasEff = rhoGasEff / (FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx) + FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx) * fs.Rv());
        const Evaluation bSolventEff = rhoSolventEff / solventRefDensity();

        // account for pressure effects
        const auto& pmiscTable = SolventModule::pmisc(elemCtx, scvIdx, timeIdx);
        const Evaluation pmisc = pmiscTable.eval(po, /*extrapolate=*/true);

        // copy the unmodified invB factors
        const Evaluation bo = fs.invB(oilPhaseIdx);
        const Evaluation bg = fs.invB(gasPhaseIdx);
        const Evaluation bs = solventInverseFormationVolumeFactor();

        const Evaluation bo_eff = pmisc * bOilEff + (1.0 - pmisc) * bo;
        const Evaluation bg_eff = pmisc * bGasEff + (1.0 - pmisc) * bg;
        const Evaluation bs_eff = pmisc * bSolventEff + (1.0 - pmisc) * bs;

        // set the densities
        fs.setDensity(oilPhaseIdx,
                      bo_eff
                      *(FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)
                        + FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)*fs.Rs()));
        fs.setDensity(gasPhaseIdx,
                      bg_eff
                      *(FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)
                        + FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)*fs.Rv()));
        solventDensity_ = bs_eff*solventRefDensity();

        // set the viscosity / mobility
        // TODO make it possible to store and modify the viscosity in fs directly

        // keep the mu*b interpolation
        Evaluation& mobo = asImp_().mobility_[oilPhaseIdx];
        muOilEff = bo_eff / (pmisc * bOilEff / muOilEff + (1.0 - pmisc) * bo / muOil);
        mobo *= muOil / muOilEff;

        Evaluation& mobg = asImp_().mobility_[gasPhaseIdx];
        muGasEff = bg_eff / (pmisc * bGasEff / muGasEff + (1.0 - pmisc) * bg / muGas);
        mobg *= muGas / muGasEff;

        // Update viscosity of solvent
        solventViscosity_ = bs_eff / (pmisc * bSolventEff / muSolventEff + (1.0 - pmisc) * bs / muSolvent);
    }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation hydrocarbonSaturation_;
    Evaluation solventSaturation_;
    Evaluation rsSolw_;
    Evaluation solventDensity_;
    Evaluation solventViscosity_;
    Evaluation solventMobility_;
    Evaluation solventInvFormationVolumeFactor_;

    Scalar solventRefDensity_;
};

template <class TypeTag>
class BlackOilSolventIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;


public:
    void solventPreSatFuncUpdate_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    { }

    void solventPostSatFuncUpdate_(const ElementContext&,
                                   unsigned,
                                   unsigned)
    { }

    void solventPvtUpdate_(const ElementContext&,
                           unsigned,
                           unsigned)
    { }

    const Evaluation& solventSaturation() const
    { throw std::runtime_error("solventSaturation() called but solvents are disabled"); }

    const Evaluation& rsSolw() const
    { throw std::runtime_error("rsSolw() called but solvents are disabled"); }

    const Evaluation& solventDensity() const
    { throw std::runtime_error("solventDensity() called but solvents are disabled"); }

    const Evaluation& solventViscosity() const
    { throw std::runtime_error("solventViscosity() called but solvents are disabled"); }

    const Evaluation& solventMobility() const
    { throw std::runtime_error("solventMobility() called but solvents are disabled"); }

    const Evaluation& solventInverseFormationVolumeFactor() const
     { throw std::runtime_error("solventInverseFormationVolumeFactor() called but solvents are disabled"); }

    const Scalar& solventRefDensity() const
     { throw std::runtime_error("solventRefDensity() called but solvents are disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilSolventExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableSolventV = getPropValue<TypeTag, Properties::EnableSolvent>()>
class BlackOilSolventExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using Toolbox = MathToolbox<Evaluation>;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;

public:
    /*!
     * \brief Method which calculates the volume flux of the polymer "phase" using the
     *        pressure potential gradient of the gas phase and the intrinsic permeability
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateVolumeFluxPerm(const ElementContext& elemCtx,
                              unsigned scvfIdx,
                              unsigned timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        PressureCallback<TypeTag> pressureCallback(elemCtx);

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const auto& faceNormal = scvf.normal();

        unsigned i = scvf.interiorIndex();
        unsigned j = scvf.exteriorIndex();

        // calculate the "raw" pressure gradient
        DimEvalVector solventPGrad;
        pressureCallback.setPhaseIndex(gasPhaseIdx);
        gradCalc.calculateGradient(solventPGrad,
                                   elemCtx,
                                   scvfIdx,
                                   pressureCallback);
        Valgrind::CheckDefined(solventPGrad);

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, i, timeIdx);
            const auto& gEx = elemCtx.problem().gravity(elemCtx, j, timeIdx);

            const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);
            const auto& intQuantsEx = elemCtx.intensiveQuantities(j, timeIdx);

            const auto& posIn = elemCtx.pos(i, timeIdx);
            const auto& posEx = elemCtx.pos(j, timeIdx);
            const auto& posFace = scvf.integrationPos();

            // the distance between the centers of the control volumes
            DimVector distVecIn(posIn);
            DimVector distVecEx(posEx);
            DimVector distVecTotal(posEx);

            distVecIn -= posFace;
            distVecEx -= posFace;
            distVecTotal -= posIn;
            Scalar absDistTotalSquared = distVecTotal.two_norm2();

            // calculate the hydrostatic pressure at the integration point of the face
            auto rhoIn = intQuantsIn.solventDensity();
            auto pStatIn = - rhoIn*(gIn*distVecIn);

            // the quantities on the exterior side of the face do not influence the
            // result for the TPFA scheme, so they can be treated as scalar values.
            Scalar rhoEx = Toolbox::value(intQuantsEx.solventDensity());
            Scalar pStatEx = - rhoEx*(gEx*distVecEx);

            // compute the hydrostatic gradient between the two control volumes (this
            // gradient exhibitis the same direction as the vector between the two
            // control volume centers and the length (pStaticExterior -
            // pStaticInterior)/distanceInteriorToExterior
            DimEvalVector f(distVecTotal);
            f *= (pStatEx - pStatIn)/absDistTotalSquared;

            // calculate the final potential gradient
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                solventPGrad[dimIdx] += f[dimIdx];

                if (!isfinite(solventPGrad[dimIdx]))
                    throw NumericalProblem("Non-finite potential gradient for solvent 'phase'");
            }
        }

        // determine the upstream and downstream DOFs
        Evaluation solventPGradNormal = 0.0;
        for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
            solventPGradNormal += solventPGrad[dimIdx]*faceNormal[dimIdx];

        if (solventPGradNormal > 0) {
            solventUpstreamDofIdx_ = j;
            solventDownstreamDofIdx_ = i;
        }
        else {
            solventUpstreamDofIdx_ = i;
            solventDownstreamDofIdx_ = j;
        }

        const auto& up = elemCtx.intensiveQuantities(solventUpstreamDofIdx_, timeIdx);

        // this is also slightly hacky because it assumes that the derivative of the
        // flux between two DOFs only depends on the primary variables in the
        // upstream direction. For non-TPFA flux approximation schemes, this is not
        // true...
        if (solventUpstreamDofIdx_ == i)
            solventVolumeFlux_ = solventPGradNormal*up.solventMobility();
        else
            solventVolumeFlux_ = solventPGradNormal*scalarValue(up.solventMobility());
    }

    /*!
     * \brief Method which calculates the volume flux of the polymer "phase" using the
     *        gas pressure potential difference between cells and transmissibilities
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateVolumeFluxTrans(const ElementContext& elemCtx,
                               unsigned scvfIdx,
                               unsigned timeIdx)
    {
        const ExtensiveQuantities& extQuants = asImp_();

        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        assert(interiorDofIdx != exteriorDofIdx);

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        unsigned I = elemCtx.globalSpaceIndex(interiorDofIdx, timeIdx);
        unsigned J = elemCtx.globalSpaceIndex(exteriorDofIdx, timeIdx);

        Scalar thpres = elemCtx.problem().thresholdPressure(I, J);
        Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        Scalar zIn = elemCtx.problem().dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        Scalar zEx = elemCtx.problem().dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);
        Scalar distZ = zIn - zEx;

        const Evaluation& rhoIn = intQuantsIn.solventDensity();
        Scalar rhoEx = Toolbox::value(intQuantsEx.solventDensity());
        const Evaluation& rhoAvg = rhoIn*0.5 + rhoEx*0.5;

        const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(gasPhaseIdx);
        Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(gasPhaseIdx));
        pressureExterior += distZ*g*rhoAvg;

        Evaluation pressureDiffSolvent = pressureExterior - pressureInterior;
        if (std::abs(scalarValue(pressureDiffSolvent)) > thpres) {
            if (pressureDiffSolvent < 0.0)
                pressureDiffSolvent += thpres;
            else
                pressureDiffSolvent -= thpres;
        }
        else
            pressureDiffSolvent = 0.0;

        if (pressureDiffSolvent > 0.0) {
            solventUpstreamDofIdx_ = exteriorDofIdx;
            solventDownstreamDofIdx_ = interiorDofIdx;
        }
        else if (pressureDiffSolvent < 0.0) {
            solventUpstreamDofIdx_ = interiorDofIdx;
            solventDownstreamDofIdx_ = exteriorDofIdx;
        }
        else {
            // pressure potential gradient is zero; force consistent upstream and
            // downstream indices over the intersection regardless of the side which it
            // is looked at.
            solventUpstreamDofIdx_ = std::min(interiorDofIdx, exteriorDofIdx);
            solventDownstreamDofIdx_ = std::max(interiorDofIdx, exteriorDofIdx);
            solventVolumeFlux_ = 0.0;
            return;
        }

        Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(solventUpstreamDofIdx_, timeIdx);
        if (solventUpstreamDofIdx_ == interiorDofIdx)
            solventVolumeFlux_ =
                up.solventMobility()
                *(-trans/faceArea)
                *pressureDiffSolvent;
        else
            solventVolumeFlux_ =
                scalarValue(up.solventMobility())
                *(-trans/faceArea)
                *pressureDiffSolvent;
    }

    unsigned solventUpstreamIndex() const
    { return solventUpstreamDofIdx_; }

    unsigned solventDownstreamIndex() const
    { return solventDownstreamDofIdx_; }

    const Evaluation& solventVolumeFlux() const
    { return solventVolumeFlux_; }

    void setSolventVolumeFlux(const Evaluation& solventVolumeFlux)
    { solventVolumeFlux_ = solventVolumeFlux; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation solventVolumeFlux_;
    unsigned solventUpstreamDofIdx_;
    unsigned solventDownstreamDofIdx_;
};

template <class TypeTag>
class BlackOilSolventExtensiveQuantities<TypeTag, false>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    void updateVolumeFluxPerm(const ElementContext&,
                              unsigned,
                              unsigned)
    { }

    void updateVolumeFluxTrans(const ElementContext&,
                              unsigned,
                              unsigned)
    { }

    unsigned solventUpstreamIndex() const
    { throw std::runtime_error("solventUpstreamIndex() called but solvents are disabled"); }

    unsigned solventDownstreamIndex() const
    { throw std::runtime_error("solventDownstreamIndex() called but solvents are disabled"); }

    const Evaluation& solventVolumeFlux() const
    { throw std::runtime_error("solventVolumeFlux() called but solvents are disabled"); }

    void setSolventVolumeFlux(const Evaluation& /* solventVolumeFlux */)
    { throw std::runtime_error("setSolventVolumeFlux() called but solvents are disabled"); }
};

} // namespace Opm

#endif
