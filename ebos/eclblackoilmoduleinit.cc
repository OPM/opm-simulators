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
#include <ebos/eclblackoilmoduleinit.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FoamadsTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FoammobTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtwsaltTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PmiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltSolubilityTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgcwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SorwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TlpmixpaTable.hpp>

#include <opm/models/blackoil/blackoilbrineparams.hh>
#include <opm/models/blackoil/blackoilextboparams.hh>
#include <opm/models/blackoil/blackoilfoamparams.hh>

#include <cassert>
#include <stdexcept>

namespace Opm {

template<bool enableSaltPrecipitation, class Scalar>
BlackOilBrineParams<Scalar>
setupBrineParams(bool enableBrine,
                 const EclipseState& eclState)
{
    BlackOilBrineParams<Scalar> params;

    // some sanity checks: if brine are enabled, the BRINE keyword must be
    // present, if brine are disabled the keyword must not be present.
    if (enableBrine && !eclState.runspec().phases().active(Phase::BRINE)) {
        throw std::runtime_error("Non-trivial brine treatment requested at compile time, but "
                                 "the deck does not contain the BRINE keyword");
    }
    else if (!enableBrine && eclState.runspec().phases().active(Phase::BRINE)) {
        throw std::runtime_error("Brine treatment disabled at compile time, but the deck "
                                 "contains the BRINE keyword");
    }

    if (!eclState.runspec().phases().active(Phase::BRINE))
        return params; // brine treatment is supposed to be disabled

    const auto& tableManager = eclState.getTableManager();

    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    params.referencePressure_.resize(numPvtRegions);

    const auto& pvtwsaltTables = tableManager.getPvtwSaltTables();

    // initialize the objects which deal with the BDENSITY keyword
    const auto& bdensityTables = tableManager.getBrineDensityTables();
    if (!bdensityTables.empty()) {
        params.bdensityTable_.resize(numPvtRegions);
        assert(numPvtRegions == bdensityTables.size());
        for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
            const auto& bdensityTable = bdensityTables[pvtRegionIdx];
            const auto& pvtwsaltTable = pvtwsaltTables[pvtRegionIdx];
            const auto& c = pvtwsaltTable.getSaltConcentrationColumn();
            params.bdensityTable_[pvtRegionIdx].setXYContainers(c, bdensityTable);
        }
    }

    if constexpr (enableSaltPrecipitation) {
        const TableContainer& permfactTables = tableManager.getPermfactTables();
        params.permfactTable_.resize(numPvtRegions);
        for (size_t i = 0; i < permfactTables.size(); ++i) {
            const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
            params.permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
        }

        const TableContainer& saltsolTables = tableManager.getSaltsolTables();
        if (!saltsolTables.empty()) {
            params.saltsolTable_.resize(numPvtRegions);
            assert(numPvtRegions == saltsolTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const SaltsolTable& saltsolTable = saltsolTables.getTable<SaltsolTable>(pvtRegionIdx );
                params.saltsolTable_[pvtRegionIdx] = saltsolTable.getSaltsolColumn().front();
            }
        }
    }

    return params;
}

template<class Scalar>
BlackOilExtboParams<Scalar> setupExtboParams(bool enableExtbo,
                                             const EclipseState& eclState)
{
    BlackOilExtboParams<Scalar> params;
    // some sanity checks: if extended BO is enabled, the PVTSOL keyword must be
    // present, if extended BO is disabled the keyword must not be present.
    if (enableExtbo && !eclState.runspec().phases().active(Phase::ZFRACTION))
        throw std::runtime_error("Extended black oil treatment requested at compile "
                                 "time, but the deck does not contain the PVTSOL keyword");
    else if (!enableExtbo && eclState.runspec().phases().active(Phase::ZFRACTION))
        throw std::runtime_error("Extended black oil treatment disabled at compile time, but the deck "
                                 "contains the PVTSOL keyword");

    if (!eclState.runspec().phases().active(Phase::ZFRACTION))
        return params; // solvent treatment is supposed to be disabled

    // pvt properties from kw PVTSOL:

    const auto& tableManager = eclState.getTableManager();
    const auto& pvtsolTables = tableManager.getPvtsolTables();

    size_t numPvtRegions = pvtsolTables.size();

    using Tabulated2DFunction = typename BlackOilExtboParams<Scalar>::Tabulated2DFunction;

    params.BO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.BG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.X_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.Y_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.VISCO_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.VISCG_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

    params.PBUB_RS_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});
    params.PBUB_RV_.resize(numPvtRegions, Tabulated2DFunction{Tabulated2DFunction::InterpolationPolicy::LeftExtreme});

    params.zLim_.resize(numPvtRegions);

    const bool extractCmpFromPvt = true; //<false>: Default values used in [*]
    params.oilCmp_.resize(numPvtRegions);
    params.gasCmp_.resize(numPvtRegions);

    for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
      const auto& pvtsolTable = pvtsolTables[regionIdx];

      const auto& saturatedTable = pvtsolTable.getSaturatedTable();
      assert(saturatedTable.numRows() > 1);

      std::vector<Scalar> oilCmp(saturatedTable.numRows(), -4.0e-9); //Default values used in [*]
      std::vector<Scalar> gasCmp(saturatedTable.numRows(), -0.08);   //-------------"-------------
      params.zLim_[regionIdx] = 0.7;                                //-------------"-------------
      std::vector<Scalar> zArg(saturatedTable.numRows(), 0.0);

      for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
        Scalar ZCO2 = saturatedTable.get("ZCO2", outerIdx);

        zArg[outerIdx] = ZCO2;

        params.BO_[regionIdx].appendXPos(ZCO2);
        params.BG_[regionIdx].appendXPos(ZCO2);

        params.RS_[regionIdx].appendXPos(ZCO2);
        params.RV_[regionIdx].appendXPos(ZCO2);

        params.X_[regionIdx].appendXPos(ZCO2);
        params.Y_[regionIdx].appendXPos(ZCO2);

        params.VISCO_[regionIdx].appendXPos(ZCO2);
        params.VISCG_[regionIdx].appendXPos(ZCO2);

        params.PBUB_RS_[regionIdx].appendXPos(ZCO2);
        params.PBUB_RV_[regionIdx].appendXPos(ZCO2);

        const auto& underSaturatedTable = pvtsolTable.getUnderSaturatedTable(outerIdx);
        size_t numRows = underSaturatedTable.numRows();

        Scalar bo0 = 0.0;
        Scalar po0 = 0.0;
        for (unsigned innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
          Scalar po = underSaturatedTable.get("P", innerIdx);
          Scalar bo = underSaturatedTable.get("B_O", innerIdx);
          Scalar bg = underSaturatedTable.get("B_G", innerIdx);
          Scalar rs = underSaturatedTable.get("RS", innerIdx)+innerIdx*1.0e-10;
          Scalar rv = underSaturatedTable.get("RV", innerIdx)+innerIdx*1.0e-10;
          Scalar xv = underSaturatedTable.get("XVOL", innerIdx);
          Scalar yv = underSaturatedTable.get("YVOL", innerIdx);
          Scalar mo = underSaturatedTable.get("MU_O", innerIdx);
          Scalar mg = underSaturatedTable.get("MU_G", innerIdx);

          if (bo0 > bo) { // This is undersaturated oil-phase for ZCO2 <= zLim ...
                          // Here we assume tabulated bo to decay beyond boiling point
              if (extractCmpFromPvt) {
                  Scalar cmpFactor = (bo-bo0)/(po-po0);
                  oilCmp[outerIdx] = cmpFactor;
                  params.zLim_[regionIdx] = ZCO2;
                  //std::cout << "### cmpFactorOil: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
              }
              break;
          } else if (bo0 == bo) { // This is undersaturated gas-phase for ZCO2 > zLim ...
                                  // Here we assume tabulated bo to be constant extrapolated beyond dew point
              if (innerIdx+1 < numRows && ZCO2<1.0 && extractCmpFromPvt) {
                Scalar rvNxt = underSaturatedTable.get("RV", innerIdx+1)+innerIdx*1.0e-10;
                Scalar bgNxt = underSaturatedTable.get("B_G", innerIdx+1);
                Scalar cmpFactor = (bgNxt-bg)/(rvNxt-rv);
                gasCmp[outerIdx] = cmpFactor;
                //std::cout << "### cmpFactorGas: " << cmpFactor << "  zLim: " << zLim_[regionIdx] << std::endl;
              }

              params.BO_[regionIdx].appendSamplePoint(outerIdx,po,bo);
              params.BG_[regionIdx].appendSamplePoint(outerIdx,po,bg);
              params.RS_[regionIdx].appendSamplePoint(outerIdx,po,rs);
              params.RV_[regionIdx].appendSamplePoint(outerIdx,po,rv);
              params.X_[regionIdx].appendSamplePoint(outerIdx,po,xv);
              params.Y_[regionIdx].appendSamplePoint(outerIdx,po,yv);
              params.VISCO_[regionIdx].appendSamplePoint(outerIdx,po,mo);
              params.VISCG_[regionIdx].appendSamplePoint(outerIdx,po,mg);
              break;
          }

          bo0 = bo;
          po0 = po;

          params.BO_[regionIdx].appendSamplePoint(outerIdx,po,bo);
          params.BG_[regionIdx].appendSamplePoint(outerIdx,po,bg);

          params.RS_[regionIdx].appendSamplePoint(outerIdx,po,rs);
          params.RV_[regionIdx].appendSamplePoint(outerIdx,po,rv);

          params.X_[regionIdx].appendSamplePoint(outerIdx,po,xv);
          params.Y_[regionIdx].appendSamplePoint(outerIdx,po,yv);

          params.VISCO_[regionIdx].appendSamplePoint(outerIdx,po,mo);
          params.VISCG_[regionIdx].appendSamplePoint(outerIdx,po,mg);

          // rs,rv -> pressure
          params.PBUB_RS_[regionIdx].appendSamplePoint(outerIdx, rs, po);
          params.PBUB_RV_[regionIdx].appendSamplePoint(outerIdx, rv, po);
        }
      }
      params.oilCmp_[regionIdx].setXYContainers(zArg, oilCmp, /*sortInput=*/false);
      params.gasCmp_[regionIdx].setXYContainers(zArg, gasCmp, /*sortInput=*/false);
    }

    // Reference density for pure z-component taken from kw SDENSITY
    const auto& sdensityTables = eclState.getTableManager().getSolventDensityTables();
    if (sdensityTables.size() == numPvtRegions) {
       params.zReferenceDensity_.resize(numPvtRegions);
       for (unsigned regionIdx = 0; regionIdx < numPvtRegions; ++ regionIdx) {
         Scalar rhoRefS = sdensityTables[regionIdx].getSolventDensityColumn().front();
         params.zReferenceDensity_[regionIdx] = rhoRefS;
       }
    }
    else
       throw std::runtime_error("Extbo:  kw SDENSITY is missing or not aligned with NTPVT\n");

    return params;
}

template<class Scalar>
BlackOilFoamParams<Scalar> setupFoamParams(bool enableFoam,
                                           const EclipseState& eclState)
{
    BlackOilFoamParams<Scalar> params;
    // some sanity checks: if foam is enabled, the FOAM keyword must be
    // present, if foam is disabled the keyword must not be present.
    if (enableFoam && !eclState.runspec().phases().active(Phase::FOAM)) {
        throw std::runtime_error("Non-trivial foam treatment requested at compile time, but "
                                 "the deck does not contain the FOAM keyword");
    }
    else if (!enableFoam && eclState.runspec().phases().active(Phase::FOAM)) {
        throw std::runtime_error("Foam treatment disabled at compile time, but the deck "
                                 "contains the FOAM keyword");
    }

    if (!eclState.runspec().phases().active(Phase::FOAM)) {
        return params; // foam treatment is supposed to be disabled
    }

    // Check that only implemented options are used.
    // We only support the default values of FOAMOPTS (GAS, TAB).
    if (eclState.getInitConfig().getFoamConfig().getTransportPhase() != Phase::GAS) {
        throw std::runtime_error("In FOAMOPTS, only GAS is allowed for the transport phase.");
    }
    if (eclState.getInitConfig().getFoamConfig().getMobilityModel() != FoamConfig::MobilityModel::TAB) {
        throw std::runtime_error("In FOAGMOPTS, only TAB is allowed for the gas mobility factor reduction model.");
    }

    const auto& tableManager = eclState.getTableManager();
    const unsigned int numSatRegions = tableManager.getTabdims().getNumSatTables();
    params.setNumSatRegions(numSatRegions);
    const unsigned int numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    params.gasMobilityMultiplierTable_.resize(numPvtRegions);

    // Get and check FOAMROCK data.
    const FoamConfig& foamConf = eclState.getInitConfig().getFoamConfig();
    if (numSatRegions != foamConf.size()) {
        throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the number of elements "
                                 "in FoamConfig, which typically corresponds to the number of records in FOAMROCK.");
    }

    // Get and check FOAMADS data.
    const auto& foamadsTables = tableManager.getFoamadsTables();
    if (foamadsTables.empty()) {
        throw std::runtime_error("FOAMADS must be specified in FOAM runs");
    }
    if (numSatRegions != foamadsTables.size()) {
        throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the "
                                 "number of FOAMADS tables.");
    }

    // Set data that vary with saturation region.
    for (std::size_t satReg = 0; satReg < numSatRegions; ++satReg) {
        const auto& rec = foamConf.getRecord(satReg);
        params.foamCoefficients_[satReg] = typename BlackOilFoamParams<Scalar>::FoamCoefficients();
        params.foamCoefficients_[satReg].fm_min = rec.minimumSurfactantConcentration();
        params.foamCoefficients_[satReg].fm_surf = rec.referenceSurfactantConcentration();
        params.foamCoefficients_[satReg].ep_surf = rec.exponent();
        params.foamRockDensity_[satReg] = rec.rockDensity();
        params.foamAllowDesorption_[satReg] = rec.allowDesorption();
        const auto& foamadsTable = foamadsTables.template getTable<FoamadsTable>(satReg);
        const auto& conc = foamadsTable.getFoamConcentrationColumn();
        const auto& ads = foamadsTable.getAdsorbedFoamColumn();
        params.adsorbedFoamTable_[satReg].setXYContainers(conc, ads);
    }

    // Get and check FOAMMOB data.
    const auto& foammobTables = tableManager.getFoammobTables();
    if (foammobTables.empty()) {
        // When in the future adding support for the functional
        // model, FOAMMOB will not be required anymore (functional
        // family of keywords can be used instead, FOAMFSC etc.).
        throw std::runtime_error("FOAMMOB must be specified in FOAM runs");
    }
    if (numPvtRegions != foammobTables.size()) {
        throw std::runtime_error("Inconsistent sizes, number of PVT regions differ from the "
                                 "number of FOAMMOB tables.");
    }

    // Set data that vary with PVT region.
    for (std::size_t pvtReg = 0; pvtReg < numPvtRegions; ++pvtReg) {
        const auto& foammobTable = foammobTables.template getTable<FoammobTable>(pvtReg);
        const auto& conc = foammobTable.getFoamConcentrationColumn();
        const auto& mobMult = foammobTable.getMobilityMultiplierColumn();
        params.gasMobilityMultiplierTable_[pvtReg].setXYContainers(conc, mobMult);
    }

    return params;
}

template BlackOilBrineParams<double>
setupBrineParams<false,double>(bool, const EclipseState&);
template BlackOilBrineParams<double>
setupBrineParams<true,double>(bool, const EclipseState&);

template BlackOilExtboParams<double>
setupExtboParams<double>(bool, const EclipseState&);

template BlackOilFoamParams<double>
setupFoamParams<double>(bool, const EclipseState&);

}
