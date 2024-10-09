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

#include <cstddef>
#include <stdexcept>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableExtbo>
void BlackOilExtboParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if extended BO is enabled, the PVTSOL keyword must be
    // present, if extended BO is disabled the keyword must not be present.
    if constexpr (enableExtbo) {
        if (!eclState.runspec().phases().active(Phase::ZFRACTION)) {
            throw std::runtime_error("Extended black oil treatment requested at compile "
                                     "time, but the deck does not contain the PVTSOL keyword");
        }
    }
    else {
        if (!enableExtbo && eclState.runspec().phases().active(Phase::ZFRACTION)) {
            throw std::runtime_error("Extended black oil treatment disabled at compile time, but the deck "
                                     "contains the PVTSOL keyword");
        }
    }

    if (!eclState.runspec().phases().active(Phase::ZFRACTION)) {
        return; // solvent treatment is supposed to be disabled
    }

    // pvt properties from kw PVTSOL:

    const auto& tableManager = eclState.getTableManager();
    const auto& pvtsolTables = tableManager.getPvtsolTables();

    std::size_t numPvtRegions = pvtsolTables.size();

    BO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    BG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    X_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    Y_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    VISCO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    VISCG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

    PBUB_RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    PBUB_RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

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

        for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++outerIdx) {
            Scalar ZCO2 = saturatedTable.get("ZCO2", outerIdx);

            zArg[outerIdx] = ZCO2;

            BO_[regionIdx].appendXPos(ZCO2);
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
            std::size_t numRows = underSaturatedTable.numRows();

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
                        Scalar cmpFactor = (bo - bo0) / (po - po0);
                        oilCmp[outerIdx] = cmpFactor;
                        zLim_[regionIdx] = ZCO2;
                        //std::cout << "### cmpFactorOil: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                    }
                    break;
                } else if (bo0 == bo) { // This is undersaturated gas-phase for ZCO2 > zLim ...
                    // Here we assume tabulated bo to be constant extrapolated beyond dew point
                    if (innerIdx+1 < numRows && ZCO2<1.0 && extractCmpFromPvt) {
                        Scalar rvNxt = underSaturatedTable.get("RV", innerIdx + 1) + innerIdx * 1.0e-10;
                        Scalar bgNxt = underSaturatedTable.get("B_G", innerIdx + 1);
                        Scalar cmpFactor = (bgNxt - bg) / (rvNxt - rv);
                        gasCmp[outerIdx] = cmpFactor;
                        //std::cout << "### cmpFactorGas: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
                    }

                    BO_[regionIdx].appendSamplePoint(outerIdx, po, bo);
                    BG_[regionIdx].appendSamplePoint(outerIdx, po, bg);
                    RS_[regionIdx].appendSamplePoint(outerIdx, po, rs);
                    RV_[regionIdx].appendSamplePoint(outerIdx, po, rv);
                    X_[regionIdx].appendSamplePoint(outerIdx, po, xv);
                    Y_[regionIdx].appendSamplePoint(outerIdx, po, yv);
                    VISCO_[regionIdx].appendSamplePoint(outerIdx, po, mo);
                    VISCG_[regionIdx].appendSamplePoint(outerIdx, po, mg);
                    break;
                }

                bo0 = bo;
                po0 = po;

                BO_[regionIdx].appendSamplePoint(outerIdx, po, bo);
                BG_[regionIdx].appendSamplePoint(outerIdx, po, bg);

                RS_[regionIdx].appendSamplePoint(outerIdx, po, rs);
                RV_[regionIdx].appendSamplePoint(outerIdx, po, rv);

                X_[regionIdx].appendSamplePoint(outerIdx, po, xv);
                Y_[regionIdx].appendSamplePoint(outerIdx, po, yv);

                VISCO_[regionIdx].appendSamplePoint(outerIdx, po, mo);
                VISCG_[regionIdx].appendSamplePoint(outerIdx, po, mg);

                       // rs,rv -> pressure
                PBUB_RS_[regionIdx].appendSamplePoint(outerIdx, rs, po);
                PBUB_RV_[regionIdx].appendSamplePoint(outerIdx, rv, po);
            }
        }
        oilCmp_[regionIdx].setXYContainers(zArg, oilCmp, /*sortInput=*/false);
        gasCmp_[regionIdx].setXYContainers(zArg, gasCmp, /*sortInput=*/false);
    }

           // Reference density for pure z-component taken from kw SDENSITY
    const auto& sdensityTables = eclState.getTableManager().getSolventDensityTables();
    if (sdensityTables.size() == numPvtRegions) {
        zReferenceDensity_.resize(numPvtRegions);
        for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++regionIdx) {
            Scalar rhoRefS = sdensityTables[regionIdx].getSolventDensityColumn().front();
            zReferenceDensity_[regionIdx] = rhoRefS;
        }
    }
    else
        throw std::runtime_error("Extbo:  kw SDENSITY is missing or not aligned with NTPVT\n");
}
#endif

#define INSTANTIATE_TYPE(T)                                                          \
    template struct BlackOilExtboParams<T>;                                          \
    template void BlackOilExtboParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilExtboParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
