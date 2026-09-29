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

#include <config.h>
#include <opm/models/blackoil/blackoilextboparams.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PmiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SorwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgcwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TlpmixpaTable.hpp>

#include <cstddef>
#include <stdexcept>

namespace {

template<bool enableExtbo>
void verifyState(const Opm::EclipseState& eclState)
{
    // some sanity checks: if extended BO is enabled, the PVTSOL keyword must be
    // present, if extended BO is disabled the keyword must not be present.
    if constexpr (enableExtbo) {
        if (!eclState.runspec().phases().active(Opm::Phase::ZFRACTION)) {
            throw std::runtime_error("Extended black oil treatment requested at compile "
                                     "time, but the deck does not contain the PVTSOL keyword");
        }
    }
    else {
        if (eclState.runspec().phases().active(Opm::Phase::ZFRACTION)) {
            throw std::runtime_error("Extended black oil treatment disabled at compile time, but the deck "
                                     "contains the PVTSOL keyword");
        }
    }
}

template<class Scalar>
std::vector<Scalar>
parseSdensity(const Opm::EclipseState& eclState,
              const std::size_t numPvtRegions)
{
    const auto& sdensityTables = eclState.getTableManager().getSolventDensityTables();
    if (sdensityTables.size() == numPvtRegions) {
        std::vector<Scalar> zReferenceDensity(numPvtRegions);
        for (std::size_t regionIdx = 0; regionIdx < numPvtRegions; ++regionIdx) {
            const Scalar rhoRefS = sdensityTables[regionIdx].getSolventDensityColumn().front();
            zReferenceDensity[regionIdx] = rhoRefS;
        }
        return zReferenceDensity;
    }
    else {
        throw std::runtime_error("Extbo: kw SDENSITY is missing or not aligned with NTPVT\n");
    }
}

}

namespace Opm {

template<class Scalar>
template<bool enableExtbo>
void BlackOilExtboParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    verifyState<enableExtbo>(eclState);
    if (!eclState.runspec().phases().active(Phase::ZFRACTION)) {
        return; // solvent treatment is supposed to be disabled
    }

    // pvt properties from kw PVTSOL:
    const auto& tableManager = eclState.getTableManager();
    const auto& pvtsolTables = tableManager.getPvtsolTables();

    const std::size_t numPvtRegions = pvtsolTables.size();

    BO_.resize(numPvtRegions);
    BG_.resize(numPvtRegions);
    RS_.resize(numPvtRegions);
    RV_.resize(numPvtRegions);
    X_.resize(numPvtRegions);
    Y_.resize(numPvtRegions);
    VISCO_.resize(numPvtRegions);
    VISCG_.resize(numPvtRegions);

    PBUB_RS_.resize(numPvtRegions);
    PBUB_RV_.resize(numPvtRegions);

    zLim_.resize(numPvtRegions);

    const bool extractCmpFromPvt = true; //<false>: Default values used in [*]
    oilCmp_.resize(numPvtRegions);
    gasCmp_.resize(numPvtRegions);

    for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++regionIdx) {
        const auto& pvtsolTable = pvtsolTables[regionIdx];

        const auto& saturatedTable = pvtsolTable.getSaturatedTable();
        assert(saturatedTable.numRows() > 1);

        std::vector<Scalar> oilCmp(saturatedTable.numRows(), -4.0e-9); //Default values used in [*]
        std::vector<Scalar> gasCmp(saturatedTable.numRows(), -0.08);   //-------------"-------------
        zLim_[regionIdx] = 0.7;                                //-------------"-------------
        std::vector<Scalar> zArg(saturatedTable.numRows(), 0.0);

        Tabulated2DFunctionBuilder BOBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder BGBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder RSBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder RVBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder XBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder YBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder VISCOBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder VISCGBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder PBUB_RSBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};
        Tabulated2DFunctionBuilder PBUB_RVBuilder{Tabulated2DFunction::InterpolationPolicy::LeftExtreme};

        for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++outerIdx) {
            Scalar ZCO2 = saturatedTable.get("ZCO2", outerIdx);

            zArg[outerIdx] = ZCO2;

            BOBuilder.appendXPos(ZCO2);
            BGBuilder.appendXPos(ZCO2);

            RSBuilder.appendXPos(ZCO2);
            RVBuilder.appendXPos(ZCO2);

            XBuilder.appendXPos(ZCO2);
            YBuilder.appendXPos(ZCO2);

            VISCOBuilder.appendXPos(ZCO2);
            VISCGBuilder.appendXPos(ZCO2);

            PBUB_RSBuilder.appendXPos(ZCO2);
            PBUB_RVBuilder.appendXPos(ZCO2);

            const auto& underSaturatedTable = pvtsolTable.getUnderSaturatedTable(outerIdx);
            const std::size_t numRows = underSaturatedTable.numRows();

            Scalar bo0 = 0.0;
            Scalar po0 = 0.0;
            for (unsigned innerIdx = 0; innerIdx < numRows; ++innerIdx) {
                Scalar po = underSaturatedTable.get("P", innerIdx);
                Scalar bo = underSaturatedTable.get("B_O", innerIdx);
                Scalar bg = underSaturatedTable.get("B_G", innerIdx);
                Scalar rs = underSaturatedTable.get("RS", innerIdx) + innerIdx * 1.0e-10;
                Scalar rv = underSaturatedTable.get("RV", innerIdx) + innerIdx * 1.0e-10;
                Scalar xv = underSaturatedTable.get("XVOL", innerIdx);
                Scalar yv = underSaturatedTable.get("YVOL", innerIdx);
                Scalar mo = underSaturatedTable.get("MU_O", innerIdx);
                Scalar mg = underSaturatedTable.get("MU_G", innerIdx);

                if (bo0 > bo) { // This is undersaturated oil-phase for ZCO2 <= zLim ...
                    // Here we assume tabulated bo to decay beyond boiling point
                    if (extractCmpFromPvt) {
                        const Scalar cmpFactor = (bo - bo0) / (po - po0);
                        oilCmp[outerIdx] = cmpFactor;
                        zLim_[regionIdx] = ZCO2;
                        //std::cout << "### cmpFactorOil: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                    }
                    break;
                } else if (bo0 == bo) { // This is undersaturated gas-phase for ZCO2 > zLim ...
                    // Here we assume tabulated bo to be constant extrapolated beyond dew point
                    if (innerIdx+1 < numRows && ZCO2<1.0 && extractCmpFromPvt) {
                        const Scalar rvNxt = underSaturatedTable.get("RV", innerIdx + 1) + innerIdx * 1.0e-10;
                        const Scalar bgNxt = underSaturatedTable.get("B_G", innerIdx + 1);
                        const Scalar cmpFactor = (bgNxt - bg) / (rvNxt - rv);
                        gasCmp[outerIdx] = cmpFactor;
                        //std::cout << "### cmpFactorGas: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                    }

                    BOBuilder.appendSamplePoint(outerIdx, po, bo);
                    BGBuilder.appendSamplePoint(outerIdx, po, bg);
                    RSBuilder.appendSamplePoint(outerIdx, po, rs);
                    RVBuilder.appendSamplePoint(outerIdx, po, rv);
                    XBuilder.appendSamplePoint(outerIdx, po, xv);
                    YBuilder.appendSamplePoint(outerIdx, po, yv);
                    VISCOBuilder.appendSamplePoint(outerIdx, po, mo);
                    VISCGBuilder.appendSamplePoint(outerIdx, po, mg);
                    break;
                }

                bo0 = bo;
                po0 = po;

                BOBuilder.appendSamplePoint(outerIdx, po, bo);
                BGBuilder.appendSamplePoint(outerIdx, po, bg);

                RSBuilder.appendSamplePoint(outerIdx, po, rs);
                RVBuilder.appendSamplePoint(outerIdx, po, rv);

                XBuilder.appendSamplePoint(outerIdx, po, xv);
                YBuilder.appendSamplePoint(outerIdx, po, yv);

                VISCOBuilder.appendSamplePoint(outerIdx, po, mo);
                VISCGBuilder.appendSamplePoint(outerIdx, po, mg);

                       // rs,rv -> pressure
                PBUB_RSBuilder.appendSamplePoint(outerIdx, rs, po);
                PBUB_RVBuilder.appendSamplePoint(outerIdx, rv, po);
            }
        }
        oilCmp_[regionIdx].setXYContainers(zArg, oilCmp, /*sortInput=*/false);
        gasCmp_[regionIdx].setXYContainers(zArg, gasCmp, /*sortInput=*/false);

        BO_[regionIdx] = std::move(BOBuilder).build();
        BG_[regionIdx] = std::move(BGBuilder).build();
        RS_[regionIdx] = std::move(RSBuilder).build();
        RV_[regionIdx] = std::move(RVBuilder).build();
        X_[regionIdx] = std::move(XBuilder).build();
        Y_[regionIdx] = std::move(YBuilder).build();
        VISCO_[regionIdx] = std::move(VISCOBuilder).build();
        VISCG_[regionIdx] = std::move(VISCGBuilder).build();
        PBUB_RS_[regionIdx] = std::move(PBUB_RSBuilder).build();
        PBUB_RV_[regionIdx] = std::move(PBUB_RVBuilder).build();
    }

    // Reference density for pure z-component taken from kw SDENSITY
    zReferenceDensity_ = parseSdensity<Scalar>(eclState, numPvtRegions);
}

#define INSTANTIATE_TYPE(T)                                                          \
    template struct BlackOilExtboParams<T>;                                          \
    template void BlackOilExtboParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilExtboParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
