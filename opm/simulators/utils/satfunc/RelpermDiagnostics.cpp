/*
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/utils/satfunc/RelpermDiagnostics.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SlgofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof3Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/WsfTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/GsfTable.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/polyhedralgrid/levelcartesianindexmapper.hh>

#ifdef HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/3d/gridview.hh>
#include <opm/simulators/flow/AluGridLevelCartesianIndexMapper.hpp>
#endif // HAVE_DUNE_ALUGRID

#include <algorithm>
#include <functional>

#include <fmt/format.h>

namespace Opm {

    bool RelpermDiagnostics::phaseCheck_(const EclipseState& es)
    {
        const auto& phases = es.runspec().phases();
        bool hasWater   = phases.active( Phase::WATER );
        bool hasGas     = phases.active( Phase::GAS );
        bool hasOil     = phases.active( Phase::OIL );
        bool hasSolvent = phases.active( Phase::SOLVENT );

        if (hasWater && !hasGas && !hasOil && !hasSolvent) {
            OpmLog::info("System:  Single phase Water system. Nothing to check");
            return false;
        }

        if (!hasWater && hasGas && !hasOil && !hasSolvent) {
            OpmLog::info("System:  Single phase Gas system. Nothing to check");
            return false;
        }

        if (!hasWater && !hasGas && hasOil && !hasSolvent) {
            OpmLog::info("System: Single phase Oil system. Nothing to check");
            return false;
        }
        if (hasWater && hasGas && !hasOil && !hasSolvent) {
            OpmLog::info("System: Water-Gas system.");
            fluidSystem_ = FluidSystem::WaterGas;
        }
        if (hasWater && hasOil && !hasGas && !hasSolvent) {
            OpmLog::info("System: Oil-Water system.");
            fluidSystem_ = FluidSystem::OilWater;
        }
        if (hasOil && hasGas && !hasWater && !hasSolvent) {
            OpmLog::info("System: Oil-Gas system.");
            fluidSystem_ = FluidSystem::OilGas;
        }
        if (hasOil && hasWater && hasGas && !hasSolvent) {
            OpmLog::info("System: Black-oil system.");
            fluidSystem_ = FluidSystem::BlackOil;
        }
        if (hasSolvent) {
            OpmLog::info("System: Solvent model.");
            fluidSystem_ = FluidSystem::Solvent;
        }
        return true;
    }

    void RelpermDiagnostics::satFamilyCheck_(const EclipseState& eclState)
    {
        const auto& tableManager = eclState.getTableManager();
        const TableContainer& swofTables = tableManager.getSwofTables();
        const TableContainer& slgofTables= tableManager.getSlgofTables();
        const TableContainer& sgofTables = tableManager.getSgofTables();
        const TableContainer& swfnTables = tableManager.getSwfnTables();
        const TableContainer& sgfnTables = tableManager.getSgfnTables();
        const TableContainer& sof3Tables = tableManager.getSof3Tables();
        const TableContainer& sof2Tables = tableManager.getSof2Tables();
        const TableContainer& sgwfnTables= tableManager.getSgwfnTables();

        const SwofletTable& swofletTable = tableManager.getSwofletTable();
        const SgofletTable& sgofletTable = tableManager.getSgofletTable();

        const TableContainer& gsfTables = tableManager.getGsfTables();
        const TableContainer& wsfTables = tableManager.getWsfTables();

        const auto& phases = eclState.runspec().phases();
        const bool waterActive = phases.active(Phase::WATER);
        const bool gasActive   = phases.active(Phase::GAS);
        const bool oilActive   = phases.active(Phase::OIL);

        std::array<bool,3> family{
            oilActive,
            true,
            !gsfTables.empty() && !wsfTables.empty()
        };
        auto& [family1, family2, family3] = family;

        // Family I test.
        if (waterActive) {
            family1 = family1 && (!swofTables.empty() || !swofletTable.empty());
        }
        if (gasActive) {
            family1 = family1 && ((!sgofTables.empty() || !sgofletTable.empty()) || !slgofTables.empty());
        }

        // Family II test.
        if (waterActive) {
            family2 = family2 && (!swfnTables.empty() || !sgwfnTables.empty());
        }
        if (oilActive) {
            family2 = family2 && (!sof3Tables.empty() || !sof2Tables.empty());
        }
        if (gasActive) {
            family2 = family2 && (!sgfnTables.empty() || !sgwfnTables.empty());
        }

        analyzeFamily(eclState, family);
    }

    void RelpermDiagnostics::analyzeFamily(const EclipseState& eclState,
                                           const std::array<bool,3>& family)
    {
        const auto& [family1, family2, family3] = family;
        if (family3) {
            const auto& phases = eclState.runspec().phases();
            const bool co2store = eclState.runspec().co2Storage();
            const bool h2store = eclState.runspec().h2Storage();
            if (!((co2store || h2store) && phases.active(Phase::GAS) && phases.active(Phase::WATER))) {
                OpmLog::info(
                    "Relative permeability input format: Saturation Family III. \n"
                    "Only valid for CO2STORE or H2STORE cases with GAS and WATER."
                );
            }
            satFamily_ = SaturationFunctionFamily::FamilyIII;
            OpmLog::info("Relative permeability input format: Saturation Family III (GSF/WSF).");
            return;
        }

        if (family1 && family2) {
            OpmLog::error(
                "Saturation families should not be mixed.\n"
                "Use either SGOF and SWOF or SGFN, SWFN and SOF3"
            );
        }

        if (!family1 && !family2) {
            OpmLog::error(
                "Saturations function must be specified using either \n"
                "family 1, family 2 or family3 keywords \n"
                "Use either SGOF and SWOF or SGFN, SWFN and SOF3."
            );
        }

        if (family1 && !family2) {
            satFamily_ = SaturationFunctionFamily::FamilyI;
            OpmLog::info("Relative permeability input format: Saturation Family I.");
        }
        if (!family1 && family2) {
            satFamily_ = SaturationFunctionFamily::FamilyII;
            OpmLog::info("Relative permeability input format: Saturation Family II.");
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_(const SwofTable& swofTables,
                                         const std::size_t satnumIdx)
    {
        const auto& sw = swofTables.getSwColumn();
        const auto& krw = swofTables.getKrwColumn();
        const auto& krow = swofTables.getKrowColumn();

        // Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SWOF table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }
        // TODO check endpoint sw.back() == 1. - Sor.
        // Check krw column.
        if (krw.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SWOF table SATNUM = {}, first value of krw should be 0.",
                satnumIdx
            ));
        }
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SWOF table SATNUM = {}, krw should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krow column.
        if (krow.front() > 1.0 || krow.back() < 0.0) {
            OpmLog::error(fmt::format(
                "In SWOF table SATNUM = {}, krow should be in range [0, 1].",
                satnumIdx
            ));
        }
        /// TODO check if run with gas.
    }

    template<>
    void RelpermDiagnostics::checkTable_<SgofTable>(const SgofTable& sgofTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& sg = sgofTables.getSgColumn();
        const auto& krg = sgofTables.getKrgColumn();
        const auto& krog = sgofTables.getKrogColumn();

        // Check sw column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGOF table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }
        if (sg.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SGOF table SATNUM = {}, first value of sg should be 0.",
                satnumIdx
            ));
        }
        // TODO check endpoint sw.back() == 1. - Sor.
        // Check krw column.
        if (krg.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SGOF table SATNUM = {}, first value of krg should be 0.",
                satnumIdx
            ));
        }
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGOF table SATNUM = {}, krg should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krow column.
        if (krog.front() > 1.0 || krog.back() < 0.0) {
            OpmLog::error(fmt::format(
                "In SGOF table SATNUM = {}, krog should be in range [0, 1].",
                satnumIdx
            ));
        }
        // TODO check if run with water.
    }

    template<>
    void RelpermDiagnostics::checkTable_<SlgofTable>(const SlgofTable& slgofTables,
                                                     const std::size_t satnumIdx)
    {
        const auto& sl = slgofTables.getSlColumn();
        const auto& krg = slgofTables.getKrgColumn();
        const auto& krog = slgofTables.getKrogColumn();

        // Check sl column.
        // TODO first value means sl = swco + sor
        if (sl.front() < 0.0 || sl.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SLGOF table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }
        if (sl.back() != 1.0) {
            OpmLog::error(fmt::format(
                "In SLGOF table SATNUM = {}, last value of sl should be 1.",
                satnumIdx
            ));
        }

        if (krg.front() > 1.0 || krg.back() < 0) {
            OpmLog::error(fmt::format(
                "In SLGOF table SATNUM = {}, krg should be in range [0, 1].",
                satnumIdx
            ));
        }
        if (krg.back() != 0.0) {
            OpmLog::error(fmt::format(
                "In SLGOF table SATNUM = {}, last value of krg should be 0.",
                satnumIdx
            ));
        }

        if (krog.front() < 0.0 || krog.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SLGOF table SATNUM = {}, krog should be in range [0, 1].",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SwfnTable>(const SwfnTable& swfnTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& sw = swfnTables.getSwColumn();
        const auto& krw = swfnTables.getKrwColumn();

        // Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SWFN table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krw column.
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SWFN table SATNUM = {}, krw should be in range [0,1].",
                satnumIdx
            ));
        }

        if (krw.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SWFN table SATNUM = {}, first value of krw should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<WsfTable>(const WsfTable& wsfTables,
                                                   const std::size_t satnumIdx)
    {
        const auto& sw = wsfTables.getSwColumn();
        const auto& krw = wsfTables.getKrwColumn();

        // Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In WSF table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krw column.
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In WSF table SATNUM = {}, krw should be in range [0,1].",
                satnumIdx
            ));
        }

        if (krw.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In WSF table SATNUM = {}, first value of krw should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SgfnTable>(const SgfnTable& sgfnTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& sg = sgfnTables.getSgColumn();
        const auto& krg = sgfnTables.getKrgColumn();

        // Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGFN table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGFN table SATNUM = {}, krg should be in range [0,1].",
                satnumIdx
            ));
        }
        if (krg.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SGFN table SATNUM = {}, first value of krg should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<GsfTable>(const GsfTable& gsfTables,
                                                   const std::size_t satnumIdx)
    {
        const auto& sg = gsfTables.getSgColumn();
        const auto& krg = gsfTables.getKrgColumn();

        // Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In GSF table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In GSF table SATNUM = {}, krg should be in range [0,1].",
                satnumIdx
            ));
        }
        if (krg.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In GSF table SATNUM = {}, first value of krg should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<Sof3Table>(const Sof3Table& sof3Tables,
                                                    const std::size_t satnumIdx)
    {
        const auto& so = sof3Tables.getSoColumn();
        const auto& krow = sof3Tables.getKrowColumn();
        const auto& krog = sof3Tables.getKrogColumn();

        // Check so column.
        // TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krow column.
        if (krow.front() < 0.0 || krow.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, krow should be in range [0,1].",
                satnumIdx
            ));
        }
        if (krow.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, first value of krow should be 0.",
                satnumIdx
            ));
        }

        // Check krog column.
        if (krog.front() < 0.0 || krog.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, krog should be in range [0,1].",
                satnumIdx
            ));
        }

        if (krog.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, first value of krog should be 0.",
                satnumIdx
            ));
        }

        if (krog.back() != krow.back()) {
            OpmLog::error(fmt::format(
                "In SOF3 table SATNUM = {}, max value of krog and krow should be the same.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<Sof2Table>(const Sof2Table& sof2Tables,
                                                    const std::size_t satnumIdx)
    {
        const auto& so = sof2Tables.getSoColumn();
        const auto& kro = sof2Tables.getKroColumn();

        // Check so column.
        // TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SOF2 table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krow column.
        if (kro.front() < 0.0 || kro.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SOF2 table SATNUM = {}, krow should be in range [0,1].",
                satnumIdx
            ));
        }
        if (kro.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SOF2 table SATNUM = {}, first value of krow should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SgwfnTable>(const SgwfnTable& sgwfnTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& sg = sgwfnTables.getSgColumn();
        const auto& krg = sgwfnTables.getKrgColumn();
        const auto& krgw = sgwfnTables.getKrgwColumn();

        // Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGWFN table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGWFN table SATNUM = {}, krg should be in range [0,1].",
                satnumIdx
            ));
        }
        if (krg.front() != 0.0) {
            OpmLog::error(fmt::format(
                "In SGWFN table SATNUM = {}, first value of krg should be 0.",
                satnumIdx
            ));
        }

        // Check krgw column.
        // TODO check saturation sw = 1. - sg
        if (krgw.front() > 1.0 || krgw.back() < 0.0) {
            OpmLog::error(fmt::format(
                "In SGWFN table SATNUM = {}, krgw should be in range [0,1].",
                satnumIdx
            ));
        }
        if (krgw.back() != 0.0) {
            OpmLog::error(fmt::format(
                "In SGWFN table SATNUM = {}, last value of krgw should be 0.",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SgcwmisTable>(const SgcwmisTable& sgcwmisTables,
                                                       const std::size_t satnumIdx)
    {
        const auto& sw = sgcwmisTables.getWaterSaturationColumn();
        const auto& sgc = sgcwmisTables.getMiscibleResidualGasColumn();

        // Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGCWMIS table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check critical gas column.
        if (sgc.front() < 0.0 || sgc.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SGCWMIS table SATNUM = {}, critical gas saturation should be in range [0,1].",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SorwmisTable>(const SorwmisTable& sorwmisTables,
                                                       const std::size_t satnumIdx)
    {
        const auto& sw = sorwmisTables.getWaterSaturationColumn();
        const auto& sor = sorwmisTables.getMiscibleResidualOilColumn();

        // Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SORWMIS table SATNUM = {}, saturation should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check critical oil column.
        if (sor.front() < 0.0 || sor.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SORWMIS table SATNUM = {}, critical oil saturation should be in range [0,1].",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<SsfnTable>(const SsfnTable& ssfnTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& frac = ssfnTables.getSolventFractionColumn();
        const auto& krgm = ssfnTables.getGasRelPermMultiplierColumn();
        const auto& krsm = ssfnTables.getSolventRelPermMultiplierColumn();

        // Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SSFN table SATNUM = {}, phase fraction should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check gas relperm multiplier column.
        if (krgm.front() < 0.0 || krgm.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SSFN table SATNUM = {}, gas relative permeability multiplier should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check solvent relperm multiplier column.
        if (krsm.front() < 0.0 || krsm.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In SSFN table SATNUM = {}, solvent relative permeability multiplier should be in range [0,1].",
                satnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<MiscTable>(const MiscTable& miscTables,
                                                    const std::size_t miscnumIdx)
    {
        const auto& frac = miscTables.getSolventFractionColumn();
        const auto& misc = miscTables.getMiscibilityColumn();

        // Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In MISC table MISCNUM = {}, phase fraction should be in range [0,1].",
                miscnumIdx
            ));
        }

        // Check miscibility column.
        if (misc.front() < 0.0 || misc.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In MISC table MISCNUM = {}, miscibility should be in range [0,1].",
                miscnumIdx
            ));
        }
    }

    template<>
    void RelpermDiagnostics::checkTable_<MsfnTable>(const MsfnTable& msfnTables,
                                                    const std::size_t satnumIdx)
    {
        const auto& frac = msfnTables.getGasPhaseFractionColumn();
        const auto& krgsm = msfnTables.getGasSolventRelpermMultiplierColumn();
        const auto& krom = msfnTables.getOilRelpermMultiplierColumn();

        // Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In MSFN table SATNUM = {}, total gas fraction should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check gas_solvent relperm multiplier column.
        if (krgsm.front() < 0.0 || krgsm.back() > 1.0) {
            OpmLog::error(fmt::format(
                "In MSFN table SATNUM = {}, gas+solvent relative permeability "
                "multiplier should be in range [0,1].",
                satnumIdx
            ));
        }

        // Check oil relperm multiplier column.
        if (krom.front() > 1.0 || krom.back() < 0.0) {
            OpmLog::error(fmt::format(
                "In MSFN table SATNUM = {}, oil relative permeability "
                "multiplier should be in range [0,1].",
                satnumIdx
            ));
        }
    }

    void RelpermDiagnostics::tableCheck_(const EclipseState& eclState)
    {
        const auto numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        OpmLog::info(fmt::format("Number of saturation regions: {}\n", numSatRegions));

        struct TableEntry
        {
            std::string name;
            const TableContainer& table;
            std::function<void(const TableContainer&, const std::size_t)> f;
        };

        auto mkop = [this]<class Type>()
        {
            return [this](const TableContainer& table, const std::size_t idx) {
                this->checkTable_<Type>(table.getTable<Type>(idx), idx + 1);
            };
        };

        const auto& tableManager = eclState.getTableManager();

        const auto tableChecks = std::array{
            TableEntry{"GSF", tableManager.getGsfTables(), mkop.operator()<GsfTable>()},
            TableEntry{"MSFN", tableManager.getMsfnTables(), mkop.operator()<MsfnTable>()},
            TableEntry{"SGCWMIS", tableManager.getSgcwmisTables(), mkop.operator()<SgcwmisTable>()},
            TableEntry{"SGFN", tableManager.getSgfnTables(), mkop.operator()<SgfnTable>()},
            TableEntry{"SGOF", tableManager.getSgofTables(), mkop.operator()<SgofTable>()},
            TableEntry{"SGWFN", tableManager.getSgwfnTables(), mkop.operator()<SgwfnTable>()},
            TableEntry{"SLGOF", tableManager.getSlgofTables(), mkop.operator()<SlgofTable>()},
            TableEntry{"SOF2", tableManager.getSof2Tables(), mkop.operator()<Sof2Table>()},
            TableEntry{"SOF3", tableManager.getSof3Tables(), mkop.operator()<Sof3Table>()},
            TableEntry{"SORWMIS", tableManager.getSorwmisTables(), mkop.operator()<SorwmisTable>()},
            TableEntry{"SSFN", tableManager.getSsfnTables(), mkop.operator()<SsfnTable>()},
            TableEntry{"SWFN", tableManager.getSwfnTables(), mkop.operator()<SwfnTable>()},
            TableEntry{"SWOF", tableManager.getSwofTables(), mkop.operator()<SwofTable>()},
            TableEntry{"WSF", tableManager.getWsfTables(), mkop.operator()<WsfTable>()},
        };

        for (std::size_t satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            std::ranges::for_each(tableChecks,
                                  [&tableManager, satnumIdx](const auto& input)
                                  {
                                      if (tableManager.hasTables(input.name)) {
                                          input.f(input.table, satnumIdx);
                                      }
                                  });
        }

        if (tableManager.hasTables("MISC")) {
            const auto& miscTables = tableManager.getMiscTables();
            const auto numMiscNumIdx = miscTables.size();
            OpmLog::info(fmt::format("Number of misc regions: {}\n", numMiscNumIdx));
            for (std::size_t miscNumIdx = 0; miscNumIdx < numMiscNumIdx; ++miscNumIdx) {
                checkTable_<MiscTable>(miscTables.getTable<MiscTable>(miscNumIdx), miscNumIdx+1);
            }
        }
    }

    void RelpermDiagnostics::unscaledEndPointsCheck_(const EclipseState& eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        const auto& runspec       = eclState.runspec();
        const auto numSatRegions = runspec.tabdims().getNumSatTables();

        if (numSatRegions < 1) {
            return;
        }

        unscaledEpsInfo_.resize(numSatRegions);

        const auto& tables = eclState.getTableManager();
        const auto& phases = runspec.phases();
        const auto tolcrit = runspec.saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto rtep =
            satfunc::getRawTableEndpoints(tables, phases, tolcrit);

        const auto rfunc =
            satfunc::getRawFunctionValues(tables, phases, rtep);

        // std::cout << "***************\nEnd-Points In all the Tables\n";
        for (std::size_t satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            this->unscaledEpsInfo_[satnumIdx]
                .extractUnscaled(rtep, rfunc, satnumIdx);

            /// Consistency check.
            if (unscaledEpsInfo_[satnumIdx].Sgu > (1. - unscaledEpsInfo_[satnumIdx].Swl)) {
                OpmLog::warning(fmt::format(
                    "In saturation table SATNUM = {}, Sgmax should not exceed 1-Swco.",
                    satnumIdx + 1
                ));
            }
            if (unscaledEpsInfo_[satnumIdx].Sgl > (1. - unscaledEpsInfo_[satnumIdx].Swu)) {
                OpmLog::warning(fmt::format(
                    "In saturation table SATNUM = {}, Sgco should not exceed 1-Swmax.",
                    satnumIdx + 1
                ));
            }

            // Krow(Sou) == Krog(Sou) for three-phase
            // means Krow(Swco) == Krog(Sgco)
            if (fluidSystem_ == FluidSystem::BlackOil) {
                blackoilChecks(eclState, satnumIdx);
            }

            // Krw(Sw=0)=Krg(Sg=0)=Krow(So=0)=Krog(So=0)=0.
            // Mobile fluid requirements
            if (((unscaledEpsInfo_[satnumIdx].Sowcr + unscaledEpsInfo_[satnumIdx].Swcr)-1) >= 0) {
                OpmLog::warning(fmt::format(
                    "In saturation table SATNUM = {}, Sowcr + Swcr should be less than 1.",
                    satnumIdx + 1
                ));
            }
            if (((unscaledEpsInfo_[satnumIdx].Sogcr +
                  unscaledEpsInfo_[satnumIdx].Sgcr +
                  unscaledEpsInfo_[satnumIdx].Swl) - 1) > 0)
            {
                OpmLog::warning(fmt::format(
                    "In saturation table SATNUM = {}, Sogcr + Sgcr + Swco should be less than 1.",
                    satnumIdx + 1
                ));
            }
        }
    }

    void RelpermDiagnostics::blackoilChecks(const EclipseState& eclState,
                                            const std::size_t satnumIdx)
    {
        const auto& tables = eclState.getTableManager();
        const TableContainer&  sgofTables = tables.getSgofTables();
        const SgofletTable&  sgofletTables = tables.getSgofletTable();
        const TableContainer& slgofTables = tables.getSlgofTables();
        const TableContainer&  swofTables = tables.getSwofTables();
        const SwofletTable&  swofletTables = tables.getSwofletTable();
        const TableContainer&  sof3Tables = tables.getSof3Tables();

        // Krow(Sou) == Krog(Sou) for three-phase
        // means Krow(Swco) == Krog(Sgco)
        double krow_value = 1e20;
        double krog_value = 1e-20;
        if (satFamily_ == SaturationFunctionFamily::FamilyI) {
          if (!sgofTables.empty()) {
              const auto& table = sgofTables.getTable<SgofTable>(satnumIdx);
              krog_value = table.evaluate( "KROG" , unscaledEpsInfo_[satnumIdx].Sgl );
          }
          else if (!sgofletTables.empty()) {
              krog_value = sgofletTables[satnumIdx].krt2_relperm;
          }
          else {
              assert(!slgofTables.empty());
              const auto& table = slgofTables.getTable<SlgofTable>(satnumIdx);
              krog_value = table.evaluate( "KROG" , unscaledEpsInfo_[satnumIdx].Sgl );
          }
          if (!swofTables.empty()) {
              const auto& table = swofTables.getTable<SwofTable>(satnumIdx);
              krow_value = table.evaluate("KROW" , unscaledEpsInfo_[satnumIdx].Swl);
          }
          else {
              assert(!swofletTables.empty());
              krow_value = swofletTables[satnumIdx].krt2_relperm;
          }
        }
        if (satFamily_ == SaturationFunctionFamily::FamilyII) {
            assert(!sof3Tables.empty());
            const auto& table = sof3Tables.getTable<Sof3Table>(satnumIdx);
            const double Sou = 1.- unscaledEpsInfo_[satnumIdx].Swl - unscaledEpsInfo_[satnumIdx].Sgl;

            krow_value = table.evaluate("KROW" , Sou);
            krog_value = table.evaluate("KROG" , Sou);
        }
        if (krow_value != krog_value) {
            OpmLog::warning(fmt::format(
                "In saturation table SATNUM = {}, Krow(Somax) should be equal to Krog(Somax).",
                satnumIdx + 1
            ));
        }
    }

    template <class LevelCartesianIndexMapper>
    void RelpermDiagnostics::diagnosis(const EclipseState& eclState,
                                       const LevelCartesianIndexMapper& levelCartesianIndexMapper)
    {
        OpmLog::info("\n===============Saturation Functions Diagnostics===============\n");
        bool doDiagnostics = phaseCheck_(eclState);
        if (!doDiagnostics) { // no diagnostics needed for single phase problems
            return;
        }
        satFamilyCheck_(eclState);
        tableCheck_(eclState);
        unscaledEndPointsCheck_(eclState);
        scaledEndPointsCheck_(eclState, levelCartesianIndexMapper);
    }

    template <class LevelCartesianIndexMapper>
    void RelpermDiagnostics::scaledEndPointsCheck_(const EclipseState& eclState,
                                                   const LevelCartesianIndexMapper& levelCartesianIndexMapper)
    {
        // All end points are subject to round-off errors, checks should account for it
        const float tolerance = 1e-6;
        const int nc = levelCartesianIndexMapper.compressedSize(0);
        const bool threepoint = eclState.runspec().endpointScaling().threepoint();
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties(eclState, false);
        const std::string tag = "Scaled endpoints";
        for (int c = 0; c < nc; ++c) {
            const std::string satnumIdx = std::to_string(epsGridProperties.satRegion(c));
            std::string cellIdx;
            {
                std::array<int, 3> ijk;
                levelCartesianIndexMapper.cartesianCoordinate(c, ijk, 0);
                cellIdx = "(" + std::to_string(ijk[0]) + ", " +
                    std::to_string(ijk[1]) + ", " +
                    std::to_string(ijk[2]) + ")";
            }
            scaledEpsInfo_[c].extractScaled(eclState, epsGridProperties, c);

            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl + tolerance)) {
                OpmLog::warning(tag,
                                fmt::format(
                                    "For scaled endpoints input, cell {} SATNUM = {}, "
                                    "SGU exceed 1.0 - SWL",
                                    cellIdx,
                                    satnumIdx
                                ));
            }

            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu + tolerance)) {
                OpmLog::warning(tag,
                                fmt::format(
                                    "For scaled endpoints input, cell {} SATNUM = {}, "
                                    "SGL exceed 1.0 - SWU",
                                    cellIdx,
                                    satnumIdx
                                ));
            }

            if (threepoint && fluidSystem_ == FluidSystem::BlackOil) {
                // Mobilility check.
                if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= (1.0 + tolerance)) {
                    OpmLog::warning(tag,
                                    fmt::format(
                                        "For scaled endpoints input, cell {} SATNUM = {}, "
                                        "SOWCR + SWCR exceed 1.0",
                                        cellIdx,
                                        satnumIdx
                                    ));
                }

            if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= (1.0 + tolerance)) {
                    OpmLog::warning(tag,
                                    fmt::format(
                                        "For scaled endpoints input, cell {} SATNUM = {}, "
                                        "SOGCR + SGCR + SWL exceed 1.0",
                                        cellIdx,
                                        satnumIdx
                                    ));
                }
            }
        }
    }

#define INSTANCE_DIAGNOSIS(...) \
    template void RelpermDiagnostics::diagnosis<LevelCartesianIndexMapper<__VA_ARGS__>>(const EclipseState&, const LevelCartesianIndexMapper<__VA_ARGS__>&); \
    template void RelpermDiagnostics::scaledEndPointsCheck_<LevelCartesianIndexMapper<__VA_ARGS__>>(const EclipseState&, const LevelCartesianIndexMapper<__VA_ARGS__>&);

    INSTANCE_DIAGNOSIS(Dune::CpGrid)
    INSTANCE_DIAGNOSIS(Dune::PolyhedralGrid<3,3>)

#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    INSTANCE_DIAGNOSIS(Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>)
#else
    INSTANCE_DIAGNOSIS(Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>)
#endif // HAVE_MPI
#endif // HAVE_DUNE_ALUGRID

} //namespace Opm
