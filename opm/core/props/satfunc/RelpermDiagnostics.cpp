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

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SgwfnTable.hpp>

namespace Opm{

    RelpermDiagnostics::RelpermDiagnostics(std::string& logFile)
    {
        streamLog_ = std::make_shared<Opm::StreamLog>(logFile, Opm::Log::DefaultMessageTypes);
    }


    RelpermDiagnostics::Counter::Counter()
        :error(0)
        ,warning(0)
        ,problem(0)
        ,bug(0)
    {
    }






    std::shared_ptr<Opm::StreamLog>
    RelpermDiagnostics::getOpmLog() const
    {
        return streamLog_;
    }




    std::vector<std::string>
    RelpermDiagnostics::getMessages() const
    {
        return messages_;
    }




    void RelpermDiagnostics::phaseCheck_(DeckConstPtr deck)
    {
        bool hasWater = deck->hasKeyword("WATER");
        bool hasGas = deck->hasKeyword("GAS");
        bool hasOil = deck->hasKeyword("OIL");
        bool hasSolvent = deck->hasKeyword("SOLVENT");
            
        if (hasWater && hasGas && !hasOil && !hasSolvent) {
            const std::string msg = "System:  Water-Gas system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::WaterGas;
        }
        if (hasWater && hasOil && !hasGas && !hasSolvent) { 
            const std::string msg = "System:  Oil-Water system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::OilWater; 
        }
        if (hasOil && hasGas && !hasWater && !hasSolvent) { 
            const std::string msg = "System:  Oil-Gas system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::OilGas; 
        }
        if (hasOil && hasWater && hasGas && !hasSolvent) {
            const std::string msg = "System:  Black-oil system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::BlackOil;
        }
        if (hasSolvent) {
            const std::string msg = "System:  Solvent model.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::Solvent;
        }
    }





    void RelpermDiagnostics::satFamilyCheck_(Opm::EclipseStateConstPtr eclState)
    {
        const auto& tableManager = eclState->getTableManager();
        const TableContainer& swofTables = tableManager.getSwofTables();
        const TableContainer& slgofTables= tableManager.getSlgofTables();
        const TableContainer& sgofTables = tableManager.getSgofTables();
        const TableContainer& swfnTables = tableManager.getSwfnTables();
        const TableContainer& sgfnTables = tableManager.getSgfnTables();
        const TableContainer& sof3Tables = tableManager.getSof3Tables();
        const TableContainer& sof2Tables = tableManager.getSof2Tables();
        const TableContainer& sgwfnTables= tableManager.getSgwfnTables();

        
        bool family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
        bool family2 = ((!swfnTables.empty() && !sgfnTables.empty()) || !sgwfnTables.empty()) && (!sof3Tables.empty() || !sof2Tables.empty());

        if (family1 && family2) {
            const std::string msg = "-- Error:   Saturation families should not be mixed.\n Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (!family1 && !family2) {
            const std::string msg = "-- Error,   Saturations function must be specified using either \n \
                             family 1 or family 2 keywords \n \
                             Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (family1 && !family2) {
            satFamily_ = SaturationFunctionFamily::FamilyI;
            const std::string msg = "Relative permeability input format: Saturation Family I.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
        } 
        if (!family1 && family2) {
            satFamily_ = SaturationFunctionFamily::FamilyII;
            const std::string msg = "Relative permeability input format: Saturation Family II.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
        }
    }


 

    void RelpermDiagnostics::tableCheck_(EclipseStateConstPtr eclState, 
                                         std::shared_ptr< const Deck > deck)
    {
        const int numSatRegions = deck->getKeyword("TABDIMS").getRecord(0).getItem("NTSFUN").get< int >(0);
        const std::string msg = "Number of saturation regions: " + std::to_string(numSatRegions) + "\n";
        std::cout << msg << std::endl;
        streamLog_->addMessage(Log::MessageType::Info, msg);
        const auto& tableManager = eclState->getTableManager();
        const TableContainer& swofTables    = tableManager.getSwofTables();
        const TableContainer& slgofTables   = tableManager.getSlgofTables();
        const TableContainer& sgofTables    = tableManager.getSgofTables();
        const TableContainer& swfnTables    = tableManager.getSwfnTables();
        const TableContainer& sgfnTables    = tableManager.getSgfnTables();
        const TableContainer& sof3Tables    = tableManager.getSof3Tables();
        const TableContainer& sof2Tables    = tableManager.getSof2Tables();
        const TableContainer& sgwfnTables   = tableManager.getSgwfnTables();
        const TableContainer& sgcwmisTables = tableManager.getSgcwmisTables();
        const TableContainer& sorwmisTables = tableManager.getSorwmisTables();
        const TableContainer& ssfnTables    = tableManager.getSsfnTables();
        const TableContainer& miscTables    = tableManager.getMiscTables();
        const TableContainer& msfnTables    = tableManager.getMsfnTables();
        
        for (int satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            if (deck->hasKeyword("SWOF")) {
                swofTableCheck_(swofTables.getTable<SwofTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SGOF")) {
                sgofTableCheck_(sgofTables.getTable<SgofTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SLGOF")) {
                slgofTableCheck_(slgofTables.getTable<SlgofTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SWFN")) {
                swfnTableCheck_(swfnTables.getTable<SwfnTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SGFN")) {
                sgfnTableCheck_(sgfnTables.getTable<SgfnTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SOF3")) {
                sof3TableCheck_(sof3Tables.getTable<Sof3Table>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SOF2")) {
                sof2TableCheck_(sof2Tables.getTable<Sof2Table>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SGWFN")) {
                sgwfnTableCheck_(sgwfnTables.getTable<SgwfnTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SGCWMIS")) {
                sgcwmisTableCheck_(sgcwmisTables.getTable<SgcwmisTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SORWMIS")) {
                sorwmisTableCheck_(sorwmisTables.getTable<SorwmisTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("SSFN")) {
                ssfnTableCheck_(ssfnTables.getTable<SsfnTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("MISC")) {
                miscTableCheck_(miscTables.getTable<MiscTable>(satnumIdx), satnumIdx+1);
            }
            if (deck->hasKeyword("MSFN")) {
                msfnTableCheck_(msfnTables.getTable<MsfnTable>(satnumIdx), satnumIdx+1);
            }
        }
    }





    void RelpermDiagnostics::swofTableCheck_(const Opm::SwofTable& swofTables, 
                                             const int satnumIdx)
    {
        const auto& sw = swofTables.getSwColumn();
        const auto& krw = swofTables.getKrwColumn();
        const auto& krow = swofTables.getKrowColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            const std::string msg = "-- Error:   In SWOF table SATNUM = "+ regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        //TODO check endpoint sw.back() == 1. - Sor.
        //Check krw column.
        if (krw.front() != 0.0) {
            const std::string msg = "-- Error:   In SWOF table SATNUM = " + regionIdx + ", first value of krw should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            const std::string msg = "-- Error:   In SWOF table SATNUM = " + regionIdx + ", krw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krow column.
        if (krow.front() > 1.0 || krow.back() < 0.0) {
            const std::string msg = "-- Error:   In SWOF table SATNUM = "+ regionIdx + ", krow should be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        ///TODO check if run with gas.
    }





    void RelpermDiagnostics::sgofTableCheck_(const Opm::SgofTable& sgofTables,
                                             const int satnumIdx)
    {
        const auto& sg = sgofTables.getSgColumn();
        const auto& krg = sgofTables.getKrgColumn();
        const auto& krog = sgofTables.getKrogColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sw column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGOF table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (sg.front() != 0.0) {
            const std::string msg = "-- Error:   In SGOF table SATNUM = " + regionIdx + ", first value of sg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        //TODO check endpoint sw.back() == 1. - Sor.
        //Check krw column.
        if (krg.front() != 0.0) {
            const std::string msg = "-- Error:   In SGOF table SATNUM = " + regionIdx + ", first value of krg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGOF table SATNUM = " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krow column.
        if (krog.front() > 1.0 || krog.back() < 0.0) {
            const std::string msg = "-- Error:   In SGOF table SATNUM = " + regionIdx + ", krog should be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        //TODO check if run with water.
    }

    void RelpermDiagnostics::slgofTableCheck_(const Opm::SlgofTable& slgofTables,
                                              const int satnumIdx) 
    {
        const auto& sl = slgofTables.getSlColumn();
        const auto& krg = slgofTables.getKrgColumn();
        const auto& krog = slgofTables.getKrogColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sl column.
        //TODO first value means sl = swco + sor
        if (sl.front() < 0.0 || sl.back() > 1.0) {
            const std::string msg = "-- Error:   In SLGOF table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (sl.back() != 1.0) {
            const std::string msg = "-- Error:   In SLGOF table SATNUM = " + regionIdx + ", last value of sl should be 1.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krg.front() > 1.0 || krg.back() < 0) {
            const std::string msg = "-- Error:   In SLGOF table SATNUM = " + regionIdx + ", krg shoule be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.back() != 0.0) {
            const std::string msg = "-- Error:   In SLGOF table SATNUM = " + regionIdx + ", last value of krg hould be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krog.front() < 0.0 || krog.back() > 1.0) {
            const std::string msg = "-- Error:   In SLGOF table SATNUM = " + regionIdx + ", krog shoule be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::swfnTableCheck_(const Opm::SwfnTable& swfnTables,
                                             const int satnumIdx)
    {
        const auto& sw = swfnTables.getSwColumn();
        const auto& krw = swfnTables.getKrwColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            const std::string msg = "-- Error:   In SWFN table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        
        //Check krw column.
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            const std::string msg = "-- Error:   In SWFN table SATNUM = " + regionIdx + ", krw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krw.front() != 0.0) {
            const std::string msg = "-- Error:   In SWFN table SATNUM = " + regionIdx + ", first value of krw should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }

    



    void RelpermDiagnostics::sgfnTableCheck_(const Opm::SgfnTable& sgfnTables,
                                             const int satnumIdx)
    {
        const auto& sg = sgfnTables.getSgColumn();
        const auto& krg = sgfnTables.getKrgColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGFN table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        
        //Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGFN table SATNUM = " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() != 0.0) {
            const std::string msg = "-- Error:   In SGFN table SATNUM = " + regionIdx + ", first value of krg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::sof3TableCheck_(const Opm::Sof3Table& sof3Tables,
                                             const int satnumIdx)
    {
        const auto& so = sof3Tables.getSoColumn();
        const auto& krow = sof3Tables.getKrowColumn();
        const auto& krog = sof3Tables.getKrogColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check so column.
        //TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krow column.
        if (krow.front() < 0.0 || krow.back() > 1.0) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", krow should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krow.front() != 0.0) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", first value of krow should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krog column.
        if (krog.front() < 0.0 || krog.back() > 1.0) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", krog should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krog.front() != 0.0) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", first value of krog should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    
        if (krog.back() != krow.back()) {
            const std::string msg = "-- Error:   In SOF3 table SATNUM = " + regionIdx + ", max value of krog and krow should be the same.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::sof2TableCheck_(const Opm::Sof2Table& sof2Tables,
                                             const int satnumIdx)
    {
        const auto& so = sof2Tables.getSoColumn();
        const auto& kro = sof2Tables.getKroColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check so column.
        //TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            const std::string msg = "-- Error:   In SOF2 table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krow column.
        if (kro.front() < 0.0 || kro.back() > 1.0) {
            const std::string msg = "-- Error:   In SOF2 table SATNUM = " + regionIdx + ", krow should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (kro.front() != 0.0) {
            const std::string msg = "-- Error:   In SOF2 table SATNUM = " + regionIdx + ", first value of krow should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::sgwfnTableCheck_(const Opm::SgwfnTable& sgwfnTables,
                                              const int satnumIdx)
    {
        const auto& sg = sgwfnTables.getSgColumn();
        const auto& krg = sgwfnTables.getKrgColumn();
        const auto& krgw = sgwfnTables.getKrgwColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGWFN table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            const std::string msg = "-- Error:   In SGWFN table SATNUM = " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() != 0.0) {
            const std::string msg = "-- Error:   In SGWFN table SATNUM = " + regionIdx + ", first value of krg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check krgw column.
        //TODO check saturation sw = 1. - sg
        if (krgw.front() > 1.0 || krgw.back() < 0.0) {
            const std::string msg = "-- Error:   In SGWFN table SATNUM = " + regionIdx + ", krgw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krgw.back() != 0.0) {
            const std::string msg = "-- Error:   In SGWFN table SATNUM = " + regionIdx + ", last value of krgw should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }



    void RelpermDiagnostics::sgcwmisTableCheck_(const Opm::SgcwmisTable& sgcwmisTables,
                                                const int satnumIdx)
    {
        const auto& sw = sgcwmisTables.getWaterSaturationColumn();
        const auto& sgc = sgcwmisTables.getMiscibleResidualGasColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            const std::string msg = "-- Error:   In SGCWMIS table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check critical gas column.
        if (sgc.front() < 0.0 || sgc.back() > 1.0) {
            const std::string msg = "-- Error:   In SGCWMIS table SATNUM = " + regionIdx + ", critical gas saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::sorwmisTableCheck_(const Opm::SorwmisTable& sorwmisTables,
                                                const int satnumIdx)
    {
        const auto& sw = sorwmisTables.getWaterSaturationColumn();
        const auto& sor = sorwmisTables.getMiscibleResidualOilColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            const std::string msg = "-- Error:   In SORWMIS table SATNUM = " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check critical oil column.
        if (sor.front() < 0.0 || sor.back() > 1.0) {
            const std::string msg = "-- Error:   In SORWMIS table SATNUM = " + regionIdx + ", critical oil saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }




    void RelpermDiagnostics::ssfnTableCheck_(const Opm::SsfnTable& ssfnTables,
                                             const int satnumIdx)
    {
        const auto& frac = ssfnTables.getSolventFractionColumn();
        const auto& krgm = ssfnTables.getGasRelPermMultiplierColumn();
        const auto& krsm = ssfnTables.getSolventRelPermMultiplierColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        //Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            const std::string msg = "-- Error:   In SSFN table SATNUM = " + regionIdx + ", phase fraction should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check gas relperm multiplier column.
        if (krgm.front() < 0.0 || krgm.back() > 1.0) {
            const std::string msg = "-- Error:   In SSFN table SATNUM = " + regionIdx + ", gas relative permeability multiplier should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check solvent relperm multiplier column.
        if (krsm.front() < 0.0 || krsm.back() > 1.0) {
            const std::string msg = "-- Error:   In SSFN table SATNUM = " + regionIdx + ", solvent relative permeability multiplier should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }






    void RelpermDiagnostics::miscTableCheck_(const Opm::MiscTable& miscTables,
                                             const int satnumIdx)
    {
        const auto& frac = miscTables.getSolventFractionColumn();
        const auto& misc = miscTables.getMiscibilityColumn();

        const std::string regionIdx = std::to_string(satnumIdx);
        //Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            const std::string msg = "-- Error:   In MISC table SATNUM = " + regionIdx + ", phase fraction should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check miscibility column.
        if (misc.front() < 0.0 || misc.back() > 1.0) {
            const std::string msg = "-- Error:   In MISC table SATNUM = " + regionIdx + ", miscibility should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::msfnTableCheck_(const Opm::MsfnTable& msfnTables,
                                             const int satnumIdx)
    {
        const auto& frac = msfnTables.getGasPhaseFractionColumn();
        const auto& krgsm = msfnTables.getGasSolventRelpermMultiplierColumn();
        const auto& krom = msfnTables.getOilRelpermMultiplierColumn();

        const std::string regionIdx = std::to_string(satnumIdx);
        //Check phase fraction column.
        if (frac.front() < 0.0 || frac.back() > 1.0) {
            const std::string msg = "-- Error:   In MSFN table SATNUM = " + regionIdx + ", total gas fraction should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check gas_solvent relperm multiplier column.
        if (krgsm.front() < 0.0 || krgsm.back() > 1.0) {
            const std::string msg = "-- Error:   In MSFN table SATNUM = " + regionIdx + ", gas+solvent relative permeability multiplier should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        //Check oil relperm multiplier column.
        if (krom.front() > 1.0 || krom.back() < 0.0) {
            const std::string msg = "-- Error:   In MSFN table SATNUM = " + regionIdx + ", oil relative permeability multiplier should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }





    void RelpermDiagnostics::unscaledEndPointsCheck_(DeckConstPtr deck,
                                                     EclipseStateConstPtr eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        const int numSatRegions = deck->getKeyword("TABDIMS").getRecord(0).getItem("NTSFUN").get< int >(0);
        unscaledEpsInfo_.resize(numSatRegions);

        const auto& tables = eclState->getTableManager();
        const TableContainer&  swofTables = tables.getSwofTables();
        const TableContainer&  sgofTables = tables.getSgofTables();
        const TableContainer& slgofTables = tables.getSlgofTables();
        const TableContainer&  sof3Tables = tables.getSof3Tables();

        // std::cout << "***************\nEnd-Points In all the Tables\n";
        for (int satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
             unscaledEpsInfo_[satnumIdx].extractUnscaled(deck, eclState, satnumIdx);
             const std::string regionIdx = std::to_string(satnumIdx + 1);
             ///Consistency check.
             if (unscaledEpsInfo_[satnumIdx].Sgu > (1. - unscaledEpsInfo_[satnumIdx].Swl)) {
                const std::string msg = "-- Warning: In saturation table SATNUM = " + regionIdx + ", Sgmax should not exceed 1-Swco.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.warning += 1;
             }
             if (unscaledEpsInfo_[satnumIdx].Sgl > (1. - unscaledEpsInfo_[satnumIdx].Swu)) {
                const std::string msg = "-- Warning: In saturation table SATNUM = " + regionIdx + ", Sgco should not exceed 1-Swmax.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.warning += 1;
             }

             //Krow(Sou) == Krog(Sou) for three-phase
             // means Krow(Swco) == Krog(Sgco)
             double krow_value = 1e20;
             double krog_value = 1e-20;
             if (fluidSystem_ == FluidSystem::BlackOil) {
                if (satFamily_ == SaturationFunctionFamily::FamilyI) {
                     if (!sgofTables.empty()) {
                         const auto& table = sgofTables.getTable<SgofTable>(satnumIdx);
                         krog_value = table.evaluate( "KROG" , unscaledEpsInfo_[satnumIdx].Sgl );
                     } else {
                         assert(!slgofTables.empty());
                         const auto& table = slgofTables.getTable<SlgofTable>(satnumIdx);
                         krog_value = table.evaluate( "KROG" , unscaledEpsInfo_[satnumIdx].Sgl );
                     }
                     {
                         const auto& table = swofTables.getTable<SwofTable>(satnumIdx);
                         krow_value = table.evaluate("KROW" , unscaledEpsInfo_[satnumIdx].Swl);
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
                     const std::string msg = "-- Warning: In saturation table SATNUM = " + regionIdx + ", Krow(Somax) should be equal to Krog(Somax).";
                     messages_.push_back(msg);
                     streamLog_->addMessage(Log::MessageType::Warning, msg);
                     counter_.warning += 1;
                 }
             }
             //Krw(Sw=0)=Krg(Sg=0)=Krow(So=0)=Krog(So=0)=0.
             //Mobile fluid requirements
            if (((unscaledEpsInfo_[satnumIdx].Sowcr + unscaledEpsInfo_[satnumIdx].Swcr)-1) >= 0) {
                const std::string msg = "-- Warning: In saturation table SATNUM = " + regionIdx + ", Sowcr + Swcr should be less than 1.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.warning += 1;
            }
            if (((unscaledEpsInfo_[satnumIdx].Sogcr + unscaledEpsInfo_[satnumIdx].Sgcr + unscaledEpsInfo_[satnumIdx].Swl) - 1 ) > 0) {
                const std::string msg = "-- Warning: In saturation table SATNUM = " + regionIdx + ", Sogcr + Sgcr + Swco should be less than 1.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.warning += 1;
            }
        }
    }






} //namespace Opm
