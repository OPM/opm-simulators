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
#ifndef EWOMS_BLACK_OIL_EXTBO_MODULE_HH
#define EWOMS_BLACK_OIL_EXTBO_MODULE_HH

#include "blackoilproperties.hh"

//#include <opm/models/io/vtkBlackOilExtboModule.hh> //TODO: Missing ...
#include <opm/models/io/vtkblackoilsolventmodule.hh>

#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/MsfnTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PmiscTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/MiscTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SorwmisTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgcwmisTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TlpmixpaTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model.
 */
template <class TypeTag, bool enableExtboV = getPropValue<TypeTag, Properties::EnableExtbo>()>
class BlackOilExtboModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Opm::SolventPvt<Scalar> SolventPvt;

    typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;
    typedef typename Opm::UniformXTabulated2DFunction<Scalar> Tabulated2DFunction;

    static constexpr unsigned zFractionIdx = Indices::zFractionIdx;
    static constexpr unsigned contiZfracEqIdx = Indices::contiZfracEqIdx;
    static constexpr unsigned enableExtbo = enableExtboV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool blackoilConserveSurfaceVolume = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume);


public:
#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the solvent module
     */
    static void initFromState(const Opm::EclipseState& eclState)
    {
        // some sanity checks: if extended BO is enabled, the PVTSOL keyword must be
        // present, if extended BO is disabled the keyword must not be present.
        if (enableExtbo && !eclState.runspec().phases().active(Phase::ZFRACTION))
            throw std::runtime_error("Extended black oil treatment requested at compile "
                                     "time, but the deck does not contain the PVTSOL keyword");
        else if (!enableExtbo && eclState.runspec().phases().active(Phase::ZFRACTION))
            throw std::runtime_error("Extended black oil treatment disabled at compile time, but the deck "
                                     "contains the PVTSOL keyword");

        if (!eclState.runspec().phases().active(Phase::ZFRACTION))
            return; // solvent treatment is supposed to be disabled

        // pvt properties from kw PVTSOL:

        const auto& tableManager = eclState.getTableManager();
        const auto& pvtsolTables = tableManager.getPvtsolTables();

        size_t numPvtRegions = pvtsolTables.size();

        B_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        BG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        X_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        Y_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        VISCO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        VISCG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

        PBUB_RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
        PBUB_RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

        for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
          const auto& pvtsolTable = pvtsolTables[regionIdx];

          const auto& saturatedTable = pvtsolTable.getSaturatedTable();
          assert(saturatedTable.numRows() > 1);

          for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
            Scalar ZCO2 = saturatedTable.get("ZCO2", outerIdx);

            B_[regionIdx].appendXPos(ZCO2);
            BG_[regionIdx].appendXPos(ZCO2);

            RS_[regionIdx].appendXPos(ZCO2);
            RV_[regionIdx].appendXPos(ZCO2);

            X_[regionIdx].appendXPos(ZCO2);
            Y_[regionIdx].appendXPos(ZCO2);

            VISCO_[regionIdx].appendXPos(ZCO2);
            VISCG_[regionIdx].appendXPos(ZCO2);

            PBUB_RS_[regionIdx].appendXPos(ZCO2);
            PBUB_RV_[regionIdx].appendXPos(ZCO2);

            const auto& underSaturatedTable = pvtsolTable.getUnderSaturatedTable(outerIdx);
            size_t numRows = underSaturatedTable.numRows();

            for (unsigned innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
              Scalar po = underSaturatedTable.get("P", innerIdx);
              Scalar bo = underSaturatedTable.get("B_O", innerIdx);
              Scalar bg = underSaturatedTable.get("B_G", innerIdx);
              Scalar rs = underSaturatedTable.get("RS", innerIdx);
              Scalar rv = underSaturatedTable.get("RV", innerIdx);
              Scalar xv = underSaturatedTable.get("XVOL", innerIdx);
              Scalar yv = underSaturatedTable.get("YVOL", innerIdx);
              Scalar mo = underSaturatedTable.get("MU_O", innerIdx);
              Scalar mg = underSaturatedTable.get("MU_G", innerIdx);

              B_[regionIdx].appendSamplePoint(outerIdx,po,bo);
              BG_[regionIdx].appendSamplePoint(outerIdx,po,bg);

              RS_[regionIdx].appendSamplePoint(outerIdx,po,rs+1.0e-10*innerIdx);
              RV_[regionIdx].appendSamplePoint(outerIdx,po,rv);

              X_[regionIdx].appendSamplePoint(outerIdx,po,xv);
              Y_[regionIdx].appendSamplePoint(outerIdx,po,yv);

              VISCO_[regionIdx].appendSamplePoint(outerIdx,po,mo);
              VISCG_[regionIdx].appendSamplePoint(outerIdx,po,mg);

              // rs,rv -> pressure
              PBUB_RS_[regionIdx].appendSamplePoint(outerIdx, rs, po);
              PBUB_RV_[regionIdx].appendSamplePoint(outerIdx, rv, po);

            }
          }
        }

        // Reference density for pure z-component taken from kw SDENSITY
        const auto& sdensityTables = eclState.getTableManager().getSolventDensityTables();
        if (sdensityTables.size() == numPvtRegions) {
           zReferenceDensity_.resize(numPvtRegions);
           for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
             Scalar rhoRefS = sdensityTables[regionIdx].getSolventDensityColumn().front();
             zReferenceDensity_[regionIdx]=rhoRefS;
           }
        }
        else
           throw std::runtime_error("Extbo:  kw SDENSITY is missing or not aligned with NTPVT\n");
    }
#endif

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        if (!enableExtbo)
            // extBO have disabled at compile time
            return;

        //Opm::VtkBlackOilExtboModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableExtbo)
            // extBO have disabled at compile time
            return;

        //model.addOutputModule(new Opm::VtkBlackOilExtboModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableExtbo)
            // extBO have disabled at compile time
            return false;

        return pvIdx == zFractionIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "z_fraction";
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableExtbo)
            return false;

        return eqIdx == contiZfracEqIdx;
    }

    static std::string eqName(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableExtbo)
            return;

        if (blackoilConserveSurfaceVolume) {
            storage[contiZfracEqIdx] =
                    Toolbox::template decay<LhsEval>(intQuants.porosity())
                    * Toolbox::template decay<LhsEval>(intQuants.yvalue())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx))
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
            if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                storage[contiZfracEqIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.porosity())
                    * Toolbox::template decay<LhsEval>(intQuants.xvalue())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs())
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(oilPhaseIdx))
                    * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(oilPhaseIdx));
            }
            // Reg. terms ...
            storage[contiZfracEqIdx] += 0.1*0.00001*(1.0-Toolbox::template decay<LhsEval>(intQuants.zFraction()))
                                      + 0.1*0.00001*Toolbox::template decay<LhsEval>(intQuants.porosity())
                                                   * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx))
                                                   * Toolbox::template decay<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx));
            storage[contiZfracEqIdx-1] += 0.1*0.00001*Toolbox::template decay<LhsEval>(intQuants.zFraction());
        }
        else {
            std::abort();
        }

    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableExtbo)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        if (blackoilConserveSurfaceVolume) {
            unsigned inIdx = extQuants.interiorIndex();

            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(gasPhaseIdx));
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            const auto& fs = up.fluidState();
            if (upIdx == inIdx) {
                flux[contiZfracEqIdx] =
                    extQuants.volumeFlux(gasPhaseIdx)
                    * (up.yvalue())
                    * fs.invB(gasPhaseIdx);
            }
            else {
                flux[contiZfracEqIdx] =
                        extQuants.volumeFlux(gasPhaseIdx)
                        * (Opm::decay<Scalar>(up.yvalue()))
                        * Opm::decay<Scalar>(fs.invB(gasPhaseIdx));
            }
            if (FluidSystem::enableDissolvedGas()) { // account for dissolved z in oil phase
                unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(oilPhaseIdx));
                const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
                const auto& fs = up.fluidState();
                if (upIdx == inIdx) {
                    flux[contiZfracEqIdx] +=
                        extQuants.volumeFlux(oilPhaseIdx)
                        * up.xvalue()
                        * fs.Rs()
                        * fs.invB(oilPhaseIdx);
                }
                else {
                    flux[contiZfracEqIdx] +=
                        extQuants.volumeFlux(oilPhaseIdx)
                        * Opm::decay<Scalar>(up.xvalue())
                        * Opm::decay<Scalar>(fs.Rs())
                        * Opm::decay<Scalar>(fs.invB(oilPhaseIdx));
                }
            }
        }
        else {
            std::abort();
        }
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar zFraction)
    {
        if (!enableExtbo)
            return;

        priVars[zFractionIdx] = zFraction;
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enableExtbo)
            return;

        // do a plain unchopped Newton update
        newPv[zFractionIdx] = oldPv[zFractionIdx] - delta[zFractionIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
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
        return std::abs(Toolbox::scalarValue(resid[contiZfracEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableExtbo)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);

        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[zFractionIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableExtbo)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);

        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[zFractionIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[zFractionIdx];
    }

    static const Tabulated2DFunction& TableX(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return X_[pvtRegionIdx];
    }

    static const Scalar xvalue(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& xvalueTable = X_[pvtRegionIdx];
        return xvalueTable.eval(z, pressure);
    }

    static const Evaluation xvalue(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& xvalueTable = X_[pvtRegionIdx];
        return xvalueTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableY(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return Y_[pvtRegionIdx];
    }

    static const Scalar yvalue(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& yvalueTable = Y_[pvtRegionIdx];
        return yvalueTable.eval(z, pressure);
    }

    static const Evaluation yvalue(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& yvalueTable = Y_[pvtRegionIdx];
        return yvalueTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TablePBUB_RS(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return PBUB_RS_[pvtRegionIdx];
    }

    static const Scalar pbubRs(unsigned pvtRegionIdx, Scalar z, Scalar rs) {
        const auto& pbubRsTable = PBUB_RS_[pvtRegionIdx];
        return pbubRsTable.eval(z, rs);
    }

    static const Evaluation pbubRs(unsigned pvtRegionIdx, Evaluation z, Evaluation rs) {
        const auto& pbubRsTable = PBUB_RS_[pvtRegionIdx];
        return pbubRsTable.eval(z, rs);
    }

    static const Tabulated2DFunction& TablePBUB_RV(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return PBUB_RV_[pvtRegionIdx];
    }

    static const Scalar pbubRv(unsigned pvtRegionIdx, Scalar z, Scalar rv) {
        const auto& pbubRvTable = PBUB_RV_[pvtRegionIdx];
        return pbubRvTable.eval(z, rv);
    }

    static const Evaluation pbubRv(unsigned pvtRegionIdx, Evaluation z, Evaluation rv) {
        const auto& pbubRvTable = PBUB_RV_[pvtRegionIdx];
        return pbubRvTable.eval(z, rv);
    }

    static const Tabulated2DFunction& TableVISCO(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return VISCO_[pvtRegionIdx];
    }

    static const Scalar visco(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& viscoTable = VISCO_[pvtRegionIdx];
        return viscoTable.eval(z, pressure);
    }

    static const Evaluation visco(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& viscoTable = VISCO_[pvtRegionIdx];
        return viscoTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableVISCG(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return VISCG_[pvtRegionIdx];
    }

    static const Scalar viscg(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& viscgTable = VISCG_[pvtRegionIdx];
        return viscgTable.eval(z, pressure);
    }

    static const Evaluation viscg(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& viscgTable = VISCG_[pvtRegionIdx];
        return viscgTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableB(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return B_[pvtRegionIdx];
    }

    static const Scalar bo(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& boTable = B_[pvtRegionIdx];
        return boTable.eval(z, pressure);
    }

    static const Evaluation bo(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& boTable = B_[pvtRegionIdx];
        return boTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableBG(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return BG_[pvtRegionIdx];
    }

    static const Scalar bg(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& bgTable = BG_[pvtRegionIdx];
        return bgTable.eval(z, pressure);
    }

    static const Evaluation bg(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& bgTable = BG_[pvtRegionIdx];
        return bgTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableRs(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return RS_[pvtRegionIdx];
    }

    static const Scalar rs(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& rsTable = RS_[pvtRegionIdx];
        return rsTable.eval(z, pressure);
    }

    static const Evaluation rs(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& rsTable = RS_[pvtRegionIdx];
        return rsTable.eval(z, pressure);
    }

    static const Tabulated2DFunction& TableRv(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return RV_[pvtRegionIdx];
    }

    static const Scalar rv(unsigned pvtRegionIdx, Scalar pressure, Scalar z) {
        const auto& rvTable = RV_[pvtRegionIdx];
        return rvTable.eval(z, pressure);
    }

    static const Evaluation rv(unsigned pvtRegionIdx, Evaluation pressure, Evaluation z) {
        const auto& rvTable = RV_[pvtRegionIdx];
        return rvTable.eval(z, pressure);
    }

    static const TabulatedFunction& TablePsat(const ElementContext& elemCtx,
                                         unsigned scvIdx,
                                         unsigned timeIdx)
    {
        unsigned pvtRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return Psat_[pvtRegionIdx];
    }

    static const Scalar PsatMult(unsigned pvtRegionIdx, Scalar z) {
        const auto& PsatTable = Psat_[pvtRegionIdx];
        return PsatTable.eval(z);
    }

    static const Evaluation PsatMult(unsigned pvtRegionIdx, Evaluation z) {
        const auto& PsatTable = Psat_[pvtRegionIdx];
        return PsatTable.eval(z);
    }

    static const Scalar referenceDensity(unsigned regionIdx) {
        return zReferenceDensity_[regionIdx];
    }

private:

    static std::vector<Tabulated2DFunction> X_;
    static std::vector<Tabulated2DFunction> Y_;
    static std::vector<Tabulated2DFunction> PBUB_RS_;
    static std::vector<Tabulated2DFunction> PBUB_RV_;
    static std::vector<Tabulated2DFunction> VISCO_;
    static std::vector<Tabulated2DFunction> VISCG_;
    static std::vector<Tabulated2DFunction> B_;
    static std::vector<Tabulated2DFunction> BG_;
    static std::vector<Tabulated2DFunction> RS_;
    static std::vector<Tabulated2DFunction> RV_;
    static std::vector<TabulatedFunction> Psat_;

    static std::vector<Scalar> zReferenceDensity_;
};

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::X_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::Y_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::PBUB_RS_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::PBUB_RV_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::VISCO_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::VISCG_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::B_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::BG_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::RS_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Tabulated2DFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::RV_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::TabulatedFunction>
BlackOilExtboModule<TypeTag, enableExtboV>::Psat_;

template <class TypeTag, bool enableExtboV>
std::vector<typename BlackOilExtboModule<TypeTag, enableExtboV>::Scalar>
BlackOilExtboModule<TypeTag, enableExtboV>::zReferenceDensity_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilExtboIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag, bool enableExtboV = GET_PROP_VALUE(TypeTag, EnableExtbo)>
class BlackOilExtboIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilExtboModule<TypeTag> SolventModule;
    typedef BlackOilExtboModule<TypeTag> ExtboModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int zFractionIdx = Indices::zFractionIdx;
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
    void zFractionUpdate_(const ElementContext& elemCtx,
                          unsigned dofIdx,
                          unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        auto& fs = asImp_().fluidState_;

        zFraction_ = priVars.makeEvaluation(zFractionIdx, timeIdx);

        visco_ = ExtboModule::visco(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        viscg_ = ExtboModule::viscg(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bo_ = ExtboModule::bo(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        bg_ = ExtboModule::bg(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);

        bz_ = ExtboModule::bg(pvtRegionIdx, fs.pressure(oilPhaseIdx), 0.99);

        rs_ = ExtboModule::rs(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        if (FluidSystem::enableVaporizedOil())
            rv_ = ExtboModule::rv(pvtRegionIdx, fs.pressure(gasPhaseIdx), zFraction_);
        else
            rv_ = 0.0;

        xvalue_ = ExtboModule::xvalue(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);
        yvalue_ = ExtboModule::yvalue(pvtRegionIdx, fs.pressure(oilPhaseIdx), zFraction_);

        pbub_ = fs.pressure(oilPhaseIdx);

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
           static const Scalar thresholdWaterFilledCell = 1.0 - 1e-6;
           Scalar Sw = 0.0;
           if (Indices::waterEnabled)
              Sw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx).value();

           if (Sw >= thresholdWaterFilledCell)
              rs_ = 0.0;  // water only, zero rs_ ...
        }

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
           rs_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           const Evaluation zLim = 0.7; //TODO: Approx. one-phase liq region ...
           if (zFraction_ > zLim) {
             pbub_ = ExtboModule::pbubRs(pvtRegionIdx, zLim, rs_);
           } else {
             pbub_ = ExtboModule::pbubRs(pvtRegionIdx, zFraction_, rs_);
           }
           //TODO: Undersaturated compressibility ...
           bo_ = ExtboModule::bo(pvtRegionIdx, pbub_, zFraction_) - 4.0e-9*(fs.pressure(oilPhaseIdx)-pbub_);

           xvalue_ = ExtboModule::xvalue(pvtRegionIdx, pbub_, zFraction_);
        }

        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
           rv_ = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
           Evaluation rvsat = ExtboModule::rv(pvtRegionIdx, pbub_, zFraction_);
           //TODO: Undersaturated compressibility ...
           bg_ = ExtboModule::bg(pvtRegionIdx, pbub_, zFraction_) + 0.08*(rvsat-rv_);

           yvalue_ = ExtboModule::yvalue(pvtRegionIdx, pbub_, zFraction_);
        }
    }

    /*!
     * \brief Update the intensive PVT properties needed to handle solvents from the
     *        primary variables.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void zPvtUpdate_(const ElementContext& elemCtx,
                           unsigned scvIdx,
                           unsigned timeIdx)
    {
        const auto& iq = asImp_();
        auto& fs = asImp_().fluidState_;

        unsigned pvtRegionIdx = iq.pvtRegionIndex();
        zRefDensity_ = ExtboModule::referenceDensity(pvtRegionIdx);
        zPureInvFormationVolumeFactor_ = 1.0/bz_;

        fs.setInvB(oilPhaseIdx, 1.0/bo_);
        fs.setInvB(gasPhaseIdx, 1.0/bg_);

        const PrimaryVariables& priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        if (priVars.primaryVarsMeaning() != PrimaryVariables::Sw_po_Rs) {
           fs.setRs(rs_);
        }
        if (priVars.primaryVarsMeaning() != PrimaryVariables::Sw_pg_Rv) {
           fs.setRv(rv_);
        }

        fs.setDensity(oilPhaseIdx,
                      fs.invB(oilPhaseIdx)
                      *(FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)
                        + (1.0-xvalue_)*fs.Rs()*FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)
                        + xvalue_*fs.Rs()*zRefDensity_ ));
        fs.setDensity(gasPhaseIdx,
                      fs.invB(gasPhaseIdx)
                      *(FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx)*(1.0-yvalue_)+yvalue_*zRefDensity_
                        + FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx)*fs.Rv()));
    }

    const Evaluation& zFraction() const
    { return zFraction_; }

    const Evaluation& xvalue() const
    { return xvalue_; }

    const Evaluation& yvalue() const
    { return yvalue_; }

    const Evaluation& pbub() const
    { return pbub_; }

    const Evaluation& visco() const
    { return visco_; }

    const Evaluation& viscg() const
    { return viscg_; }

    const Evaluation& bo() const
    { return bo_; }

    const Evaluation& bg() const
    { return bg_; }

    const Evaluation& rs() const
    { return rs_; }

    const Evaluation& rv() const
    { return rv_; }

    const Evaluation& zPureInvFormationVolumeFactor() const
    { return zPureInvFormationVolumeFactor_; }

    // This could be stored pr pvtRegion instead
    const Scalar& zRefDensity() const
    { return zRefDensity_; }

private:


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation hydrocarbonSaturation_;
    Evaluation zFraction_;
    Evaluation xvalue_;
    Evaluation yvalue_;
    Evaluation pbub_;
    Evaluation visco_;
    Evaluation viscg_;
    Evaluation bo_;
    Evaluation bg_;
    Evaluation bz_;
    Evaluation rs_;
    Evaluation rv_;
    Evaluation zPureInvFormationVolumeFactor_;

    Scalar zRefDensity_;
};

template <class TypeTag>
class BlackOilExtboIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


public:

    void zPvtUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                     unsigned scvIdx OPM_UNUSED,
                     unsigned timeIdx OPM_UNUSED)
    { }

    void zFractionUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                          unsigned scvIdx OPM_UNUSED,
                          unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& xvalue() const
    { throw std::runtime_error("xvalue() called but extbo is disabled"); }

    const Evaluation& yvalue() const
    { throw std::runtime_error("yvalue() called but extbo is disabled"); }

    const Evaluation& visco() const
    { throw std::runtime_error("visco() called but extbo is disabled"); }

    const Evaluation& viscg() const
    { throw std::runtime_error("viscg() called but extbo is disabled"); }

    const Evaluation& rs() const
    { throw std::runtime_error("rs() called but extbo is disabled"); }

    const Evaluation& rv() const
    { throw std::runtime_error("rv() called but extbo is disabled"); }

    const Evaluation& zPureInvFormationVolumeFactor() const
    { throw std::runtime_error("zPureInvFormationVolumeFactor() called but extbo is disabled"); }

    const Evaluation& zFraction() const
    { throw std::runtime_error("zFraction() called but extbo is disabled"); }

    const Evaluation& zInverseFormationVolumeFactor() const
     { throw std::runtime_error("zInverseFormationVolumeFactor() called but extbo is disabled"); }

    const Scalar& zRefDensity() const
     { throw std::runtime_error("zRefDensity() called but extbo is disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilExtboExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableExtboV = GET_PROP_VALUE(TypeTag, EnableExtbo)>
class BlackOilExtboExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilExtboExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:

};

} // namespace Opm

#endif
