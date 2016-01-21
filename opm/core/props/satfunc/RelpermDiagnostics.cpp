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
#include <opm/core/utility/compressedToCartesian.hpp>

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


    void RelpermDiagnostics::diagnosis(Opm::EclipseStateConstPtr eclState,
                                       Opm::DeckConstPtr deck,
                                       const UnstructuredGrid& grid)
    {
        std::cout << "\n\n***************Relperm Diagnostics***************\n\n";
        phaseCheck_(deck);
        satFamilyCheck_(eclState);
        tableCheck_(eclState, deck);
        unscaledEndPointsCheck_(deck, eclState);
        scaledEndPointsCheck_(deck, eclState, grid);
        if (!messages_.empty()) {
            std::sort(messages_.begin(), messages_.end());
            auto it = std::unique(messages_.begin(), messages_.end());
            messages_.erase(it, messages_.end());
            std::cout << std::endl;
            for (const auto& x : messages_) {
                std::cout << "--" << x << std::endl;
            }
        }
        int limits = 0;
        if (!scaled_messages_.empty()) {
            std::cout << std::endl;
            for (const auto& x : scaled_messages_) {
                if (limits < 10) {
                    std::cout << "-- " << x << std::endl;
                    limits++;
                } else {
                    std::cout << "\nMore inconsistencies exist. Check saturation function input and LOGFILE!" << std::endl;
                    break;
                }
            }
        }

        std::string summary_msg = "\n\nError summary:" + 
            std::string("\nWarnings          " + std::to_string(counter_.warning)) +
            std::string("\nProblems          " + std::to_string(counter_.problem)) +
            std::string("\nErrors            " + std::to_string(counter_.error)) + 
            std::string("\nBugs              " + std::to_string(counter_.bug))+ "\n";
        streamLog_->addMessage(Log::MessageType::Info, summary_msg);
        std::cout << summary_msg << std::endl;
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

        if (hasWater && hasGas && !hasOil) {
            std::string msg = "System:  Water-Gas system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::WaterGas;
        }
        if (hasWater && hasOil && !hasGas) { 
            std::string msg = "System:  Oil-Water system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::OilWater; 
        }
        if (hasOil && hasGas && !hasWater) { 
            std::string msg = "System:  Oil-Gas system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::OilGas; 
        }
        if (hasOil && hasWater && hasGas) {
            std::string msg = "System:  Black-oil system.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
            fluidSystem_ = FluidSystem::BlackOil;
        }
    }





    void RelpermDiagnostics::satFamilyCheck_(Opm::EclipseStateConstPtr eclState)
    {
        const auto& tableManager = eclState->getTableManager();
        const TableContainer& swofTables = tableManager->getSwofTables();
        const TableContainer& slgofTables= tableManager->getSlgofTables();
        const TableContainer& sgofTables = tableManager->getSgofTables();
        const TableContainer& swfnTables = tableManager->getSwfnTables();
        const TableContainer& sgfnTables = tableManager->getSgfnTables();
        const TableContainer& sof3Tables = tableManager->getSof3Tables();
        const TableContainer& sof2Tables = tableManager->getSof2Tables();
        const TableContainer& sgwfnTables= tableManager->getSgwfnTables();

        bool family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
        bool family2 = !swfnTables.empty() && !sgfnTables.empty() && (!sof3Tables.empty() || !sof2Tables.empty()) && !sgwfnTables.empty();

        if (family1 && family2) {
            std::string msg = "Error:   Saturation families should not be mixed.\n Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (!family1 && !family2) {
            std::string msg = "Error,   Saturations function must be specified using either \n \
                             family 1 or family 2 keywords \n \
                             Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (family1 && !family2) {
            satFamily_ = SaturationFunctionFamily::FamilyI;
            std::string msg = "Relative permeability input format: Saturation Family I.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
        } 
        if (!family1 && family2) {
            satFamily_ = SaturationFunctionFamily::FamilyII;
            std::string msg = "Relative permeambility input format: Saturation Family II.";
            std::cout << msg << std::endl;
            streamLog_->addMessage(Log::MessageType::Info, msg);
        }
    }


 

    void RelpermDiagnostics::tableCheck_(EclipseStateConstPtr eclState, 
                                         DeckConstPtr deck)
    {
        int numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        std::string msg = "Number of saturation regions: " + std::to_string(numSatRegions) + "\n";
        std::cout << msg << std::endl;
        streamLog_->addMessage(Log::MessageType::Info, msg);
        const auto& tableManager = eclState->getTableManager();
        const TableContainer& swofTables = tableManager->getSwofTables();
        const TableContainer& slgofTables= tableManager->getSlgofTables();
        const TableContainer& sgofTables = tableManager->getSgofTables();
        const TableContainer& swfnTables = tableManager->getSwfnTables();
        const TableContainer& sgfnTables = tableManager->getSgfnTables();
        const TableContainer& sof3Tables = tableManager->getSof3Tables();
        const TableContainer& sof2Tables = tableManager->getSof2Tables();
        const TableContainer& sgwfnTables= tableManager->getSgwfnTables();

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
        }
    }





    void RelpermDiagnostics::swofTableCheck_(const Opm::SwofTable& swofTables, 
                                             const int satnumIdx)
    {
        const auto& sw = swofTables.getSwColumn();
        const auto& krw = swofTables.getKrwColumn();
        const auto& krow = swofTables.getKrowColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        ///Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            std::string msg = "Error:   In SWOF table region "+ regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        ///TODO check endpoint sw.back() == 1. - Sor.
        ///Check krw column.
        if (krw.front() != 0.0) {
            std::string msg = "Error:   In SWOF table region " + regionIdx + ", first value of krw should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            std::string msg = "Error:   In SWOF table region " + regionIdx + ", krw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krow column.
        if (krow.front() > 1.0 || krow.back() < 0.0) {
            std::string msg = "Error:   In SWOF table region "+ regionIdx + ", krow should be in range [0, 1].";
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
        ///Check sw column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            std::string msg = "Error:   In SGOF table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (sg.front() != 0.0) {
            std::string msg = "Error:   In SGOF table region " + regionIdx + ", first value of sg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        ///TODO check endpoint sw.back() == 1. - Sor.
        ///Check krw column.
        if (krg.front() != 0.0) {
            std::string msg = "Error:   In SGOF table region " + regionIdx + ", first value of krg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            std::string msg = "Error:   In SGOF table region " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krow column.
        if (krog.front() > 1.0 || krog.back() < 0.0) {
            std::string msg = "Error:   In SGOF table region " + regionIdx + ", krog should be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        ///TODO check if run with water.
    }

    void RelpermDiagnostics::slgofTableCheck_(const Opm::SlgofTable& slgofTables,
                                              const int satnumIdx) 
    {
        const auto& sl = slgofTables.getSlColumn();
        const auto& krg = slgofTables.getKrgColumn();
        const auto& krog = slgofTables.getKrogColumn();
        const std::string regionIdx = std::to_string(satnumIdx);
        ///Check sl column.
        ///TODO first value means sl = swco + sor
        if (sl.front() < 0.0 || sl.back() > 1.0) {
            std::string msg = "Error:   In SLGOF table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (sl.back() != 1.0) {
            std::string msg = "Error:   In SLGOF table region " + regionIdx + ", last value of sl should be 1.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krg.front() > 1.0 || krg.back() < 0) {
            std::string msg = "Error:   In SLGOF table region " + regionIdx + ", krg shoule be in range [0, 1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.back() != 0.0) {
            std::string msg = "Error:   In SLGOF table region " + regionIdx + ", last value of krg hould be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krog.front() < 0.0 || krog.back() > 1.0) {
            std::string msg = "Error:   In SLGOF table region " + regionIdx + ", krog shoule be in range [0, 1].";
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
        ///Check sw column.
        if (sw.front() < 0.0 || sw.back() > 1.0) {
            std::string msg = "Error:   In SWFN table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        
        ///Check krw column.
        if (krw.front() < 0.0 || krw.back() > 1.0) {
            std::string msg = "Error:   In SWFN table region " + regionIdx + ", krw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krw.front() != 0.0) {
            std::string msg = "Error:   In SWFN table region " + regionIdx + ", first value of krw should be 0.";
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
        ///Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            std::string msg = "Error:   In SGFN table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        
        ///Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            std::string msg = "Error:   In SGFN table region " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() != 0.0) {
            std::string msg = "Error:   In SGFN table region " + regionIdx + ", first value of krg should be 0.";
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
        ///Check so column.
        ///TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krow column.
        if (krow.front() < 0.0 || krow.back() > 1.0) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", krow should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krow.front() != 0.0) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", first value of krow should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krog column.
        if (krog.front() < 0.0 || krog.back() > 1.0) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", krog should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        if (krog.front() != 0.0) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", first value of krog should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    
        if (krog.back() != krow.back()) {
            std::string msg = "Error:   In SOF3 table region " + regionIdx + ", max value of krog and krow should be the same.";
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
        ///Check so column.
        ///TODO: The max so = 1 - Swco
        if (so.front() < 0.0 || so.back() > 1.0) {
            std::string msg = "Error:   In SOF2 table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krow column.
        if (kro.front() < 0.0 || kro.back() > 1.0) {
            std::string msg = "Error:   In SOF2 table region " + regionIdx + ", krow should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (kro.front() != 0.0) {
            std::string msg = "Error:   In SOF2 table region " + regionIdx + ", first value of krow should be 0.";
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
        ///Check sg column.
        if (sg.front() < 0.0 || sg.back() > 1.0) {
            std::string msg = "Error:   In SGWFN table region " + regionIdx + ", saturation should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krg column.
        if (krg.front() < 0.0 || krg.back() > 1.0) {
            std::string msg = "Error:   In SGWFN table region " + regionIdx + ", krg should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krg.front() != 0.0) {
            std::string msg = "Error:   In SGWFN table region " + regionIdx + ", first value of krg should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }

        ///Check krgw column.
        ///TODO check saturation sw = 1. - sg
        if (krgw.front() > 1.0 || krgw.back() < 0.0) {
            std::string msg = "Error:   In SGWFN table region " + regionIdx + ", krgw should be in range [0,1].";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
        if (krgw.back() != 0.0) {
            std::string msg = "Error:   In SGWFN table region " + regionIdx + ", last value of krgw should be 0.";
            messages_.push_back(msg);
            streamLog_->addMessage(Log::MessageType::Error, msg);
            counter_.error += 1;
        }
    }




    void RelpermDiagnostics::unscaledEndPointsCheck_(DeckConstPtr deck,
                                                     EclipseStateConstPtr eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        int numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        unscaledEpsInfo_.resize(numSatRegions);

        auto tables = eclState->getTableManager();
        const TableContainer&  swofTables = tables->getSwofTables();
        const TableContainer&  sgofTables = tables->getSgofTables();
        const TableContainer& slgofTables = tables->getSlgofTables();
        const TableContainer&  sof3Tables = tables->getSof3Tables();

        // std::cout << "***************\nEnd-Points In all the Tables\n";
        for (int satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
             unscaledEpsInfo_[satnumIdx].extractUnscaled(deck, eclState, satnumIdx);
  
             ///Consistency check.
             if (unscaledEpsInfo_[satnumIdx].Sgu > (1. - unscaledEpsInfo_[satnumIdx].Swl)) {
                 std::string msg = "Warning: In saturation region " + std::to_string(satnumIdx+1) + ", Sgmax should not exceed 1-Swco.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
             }
             if (unscaledEpsInfo_[satnumIdx].Sgl > (1. - unscaledEpsInfo_[satnumIdx].Swu)) {
                 std::string msg = "Warning: In saturation region " + std::to_string(satnumIdx+1) + ", Sgco should not exceed 1-Swmax.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
             }

             ///Krow(Sou) == Krog(Sou) for three-phase
             /// means Krow(Swco) == Krog(Sgco)
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
                     std::string msg = "Warning: In region " + std::to_string(satnumIdx+1) + ", Krow(sSomax) should equal Krog(Somax).";
                     messages_.push_back(msg);
                     streamLog_->addMessage(Log::MessageType::Warning, msg);
                     counter_.error += 1;
                 }
             }
             ///Krw(Sw=0)=Krg(Sg=0)=Krow(So=0)=Krog(So=0)=0.
             ///Mobile fluid requirements
            if (((unscaledEpsInfo_[satnumIdx].Sowcr + unscaledEpsInfo_[satnumIdx].Swcr)-1) >= 0) {
                std::string msg = "Warning: In saturation region " + std::to_string(satnumIdx+1) + ", Sowcr + Swcr should be less than 1.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
            }
            if (((unscaledEpsInfo_[satnumIdx].Sogcr + unscaledEpsInfo_[satnumIdx].Sgcr + unscaledEpsInfo_[satnumIdx].Swl) - 1 ) > 0) {
                std::string msg = "Warning: In saturation rgion " + std::to_string(satnumIdx+1) + ", Sogcr + Sgcr + Swco should be less than 1.";
                messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
            }
        }
    }





    void RelpermDiagnostics::scaledEndPointsCheck_(DeckConstPtr deck,
                                                   EclipseStateConstPtr eclState,
                                                   const UnstructuredGrid& grid)
    {
        const int nc = Opm::UgGridHelpers::numCells(grid);
        //std::vector<int> compressedToCartesianIdx(nc);
        const auto& global_cell = Opm::UgGridHelpers::globalCell(grid);
        const auto dims = Opm::UgGridHelpers::cartDims(grid);
        //        for (int cell = 0; cell < nc; ++cell) {
        //            if (global_cell) {
        //                compressedToCartesianIdx[cell] = global_cell[cell];
        //            } else {
        //                compressedToCartesianIdx[cell] = cell;
        //            }
        //        }
        const auto& compressedToCartesianIdx = Opm::compressedToCartesian(nc, global_cell);
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);       

        for (int c = 0; c < nc; ++c) {
            int cartIdx = compressedToCartesianIdx[c];
            std::array<int, 3> ijk;
            ijk[0] = cartIdx % dims[0];
            ijk[1] = (cartIdx / dims[0]) % dims[1];
            ijk[2] = cartIdx / dims[0] / dims[1];
            std::string cellIdx = "(" + std::to_string(ijk[0]) + ", " + 
                                   std::to_string(ijk[1]) + ", " +
                                   std::to_string(ijk[2]) + ")";
            scaledEpsInfo_[c].extractScaled(epsGridProperties, cartIdx);

            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl)) {
                std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SGU exceed 1.0 - SWL";
                scaled_messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
            }
            
            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu)) {
                std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SGL exceed 1.0 - SWU";
                scaled_messages_.push_back(msg);
                streamLog_->addMessage(Log::MessageType::Warning, msg);
                counter_.error += 1;
            }

            if (deck->hasKeyword("SCALECRS") && fluidSystem_ == FluidSystem::BlackOil) {
                // Mobilility check.
                if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= 1.0) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SOWCR + SWCR exceed 1.0";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.error += 1;
                }

                if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= 1.0) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SOGCR + SGCR + SWL exceed 1.0";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.error += 1;
                }
            }
            ///Following rules come from NEXUS.
            if (fluidSystem_ != FluidSystem::WaterGas) {
                if (scaledEpsInfo_[c].Swl > scaledEpsInfo_[c].Swcr) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SWL > SWCR";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }

                if (scaledEpsInfo_[c].Swcr > scaledEpsInfo_[c].Sowcr) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SWCR > SOWCR";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Opm::Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }
            
                if (scaledEpsInfo_[c].Sowcr > scaledEpsInfo_[c].Swu) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SOWCR > SWU";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }
            }

            if (fluidSystem_ != FluidSystem::OilWater) {
                if (scaledEpsInfo_[c].Sgl > scaledEpsInfo_[c].Sgcr) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SGL > SGCR";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }
            }

            if (fluidSystem_ != FluidSystem::BlackOil) {
                if (scaledEpsInfo_[c].Sgcr > scaledEpsInfo_[c].Sogcr) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SGCR > SOGCR";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }

                if (scaledEpsInfo_[c].Sogcr > scaledEpsInfo_[c].Sgu) {
                    std::string msg = "Warning: For scaled endpoints input, cell" + cellIdx + " SOGCR > SGU";
                    scaled_messages_.push_back(msg);
                    streamLog_->addMessage(Log::MessageType::Warning, msg);
                    counter_.warning += 1;
                }
            }
        } 
    }

} //namespace Opm
