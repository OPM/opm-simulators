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

#ifndef OPM_RELPERMDIAGNOSTICS_HEADER_INCLUDED
#define OPM_RELPERMDIAGNOSTICS_HEADER_INCLUDED

#include <vector>
#include <utility>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

namespace Opm {
    enum FluidSystem {
        OilWater,
        OilGas,
        WaterGas,
        BlackOil
    };

    enum SaturationFunctionFamily {
        FamilyI,
        FamilyII,
        NoFamily
    };


    ///This class is intend to be a relpmer diganostics, to detect
    ///wrong input of relperm table and endpoints.
    class RelpermDiagnostics 
    {
    public:
        RelpermDiagnostics();
        RelpermDiagnostics(EclipseStateConstPtr eclstate);

        void diagnosis(EclipseStateConstPtr eclState,
                       DeckConstPtr deck,
                       const UnstructuredGrid& grid);

    private:
        FluidSystem fluidSystem_;
        
        SaturationFunctionFamily satFamily_;

        std::vector<Opm::EclEpsScalingPointsInfo<double> > unscaledEpsInfo_;
        std::vector<Opm::EclEpsScalingPointsInfo<double> > scaledEpsInfo_;

        std::vector<std::string> messager_;
        
        ///Display all the tables.
        void keywordsDisplay_(EclipseStateConstPtr eclState);
        
        ///Check the phase that used.
        void phaseCheck_(EclipseStateConstPtr eclState, 
                         DeckConstPtr deck);

        
        ///Check saturation family I and II.
        void satFamilyCheck_(EclipseStateConstPtr eclState);
 
        ///Check saturation tables.
        void tableCheck_(EclipseStateConstPtr eclState, 
                         DeckConstPtr deck);

        ///Check endpoints in the saturation tables.
        void unscaledEndPointsCheck_(DeckConstPtr deck,
                                     EclipseStateConstPtr eclState);

        void scaledEndPointsCheck_(DeckConstPtr deck,
                                   EclipseStateConstPtr eclState,
                                   const UnstructuredGrid& grid);

        ///For every table, need to deal with case by case.
        void swofTableCheck_(const Opm::SwofTable& swofTables);
        void sgofTableCheck_(const Opm::SgofTable& sgofTables);
        void slgofTableCheck_(const Opm::SlgofTable& slgofTables);
        void swfnTableCheck_(const Opm::SwfnTable& swfnTables);
        void sgfnTableCheck_(const Opm::SgfnTable& sgfnTables);
        void sof3TableCheck_(const Opm::Sof3Table& sof3Tables);
        void sof2TableCheck_(const Opm::Sof2Table& sof2Tables);
        void sgwfnTableCheck_(const Opm::SgwfnTable& sgwfnTables);
    };

    RelpermDiagnostics::RelpermDiagnostics()
    {}

    RelpermDiagnostics::RelpermDiagnostics(Opm::EclipseStateConstPtr eclState)
    {}


    void RelpermDiagnostics::diagnosis(Opm::EclipseStateConstPtr eclState,
                                       Opm::DeckConstPtr deck,
                                       const UnstructuredGrid& grid)
    {
        std::cout << "***************Relperm Diagnostics***************\n";
        phaseCheck_(eclState, deck);
        satFamilyCheck_(eclState);
        tableCheck_(eclState, deck);
        unscaledEndPointsCheck_(deck, eclState);
        scaledEndPointsCheck_(deck, eclState, grid);
        if (!messager_.empty()) {
            int counter = 1;
            std::cout << "***************\nProblem found:\n";
            for (const auto& x : messager_) {
                std::cout << counter << ".  " << x << std::endl;
                counter++;
            }
        } else {
            std::cout << "****************\nConsistency check all passed!" << std::endl;
        }
        std::cout << "********************************************************\n";
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
            std::string s = "Saturation families should not be mixed.\n Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messager_.push_back(s);
        }

        if (!family1 && !family2) {
            std::string s = "Saturations function must be specified using either \n \
                             family 1 or family 2 keywords \n \
                             Use either SGOF and SWOF or SGFN, SWFN and SOF3.";
            messager_.push_back(s);
        }

        if (family1 && !family2) {
            satFamily_ = SaturationFunctionFamily::FamilyI;
            std::cout << "relperm: Saturation Family I." << std::endl;
        } 
        if (!family1 && family2) {
            satFamily_ = SaturationFunctionFamily::FamilyII;
            std::cout << "relperm: Saturation Family II." << std::endl;
        }
    }


    void RelpermDiagnostics::phaseCheck_(EclipseStateConstPtr eclState,
                                         DeckConstPtr deck)
    {
        bool hasWater = deck->hasKeyword("WATER");
        bool hasGas = deck->hasKeyword("GAS");
        bool hasOil = deck->hasKeyword("OIL");

        if (hasWater && hasGas && !hasOil) {
            std::cout << "System:  Water-Gas system." << std::endl;
            fluidSystem_ = FluidSystem::WaterGas;
        }
        if (hasWater && hasOil && !hasGas) { 
            std::cout << "System:  Oil-Water system." << std::endl;
            fluidSystem_ = FluidSystem::OilWater; 
        }
        if (hasOil && hasGas && !hasWater) { 
            std::cout << "System:  Oil-Gas system." << std::endl;
            fluidSystem_ = FluidSystem::OilGas; 
        }
        if (hasOil && hasWater && hasGas) {
            std::cout << "System:  Black-oil system." << std::endl;
            fluidSystem_ = FluidSystem::BlackOil;
        }
    } 

    void RelpermDiagnostics::tableCheck_(EclipseStateConstPtr eclState, 
                                         DeckConstPtr deck)
    {
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));
        const auto& tableManager = eclState->getTableManager();
        const TableContainer& swofTables = tableManager->getSwofTables();
        const TableContainer& slgofTables= tableManager->getSlgofTables();
        const TableContainer& sgofTables = tableManager->getSgofTables();
        const TableContainer& swfnTables = tableManager->getSwfnTables();
        const TableContainer& sgfnTables = tableManager->getSgfnTables();
        const TableContainer& sof3Tables = tableManager->getSof3Tables();
        const TableContainer& sof2Tables = tableManager->getSof2Tables();
        const TableContainer& sgwfnTables= tableManager->getSgwfnTables();

        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            if (deck->hasKeyword("SWOF")) {
                swofTableCheck_(swofTables.getTable<SwofTable>(satnumIdx));
            }
            if (deck->hasKeyword("SGOF")) {
                sgofTableCheck_(sgofTables.getTable<SgofTable>(satnumIdx));
            }
            if (deck->hasKeyword("SLGOF")) {
                slgofTableCheck_(slgofTables.getTable<SlgofTable>(satnumIdx));
            }
            if (deck->hasKeyword("SWFN")) {
                swfnTableCheck_(swfnTables.getTable<SwfnTable>(satnumIdx));
            }
            if (deck->hasKeyword("SGFN")) {
                sgfnTableCheck_(sgfnTables.getTable<SgfnTable>(satnumIdx));
            }
            if (deck->hasKeyword("SOF3")) {
                sof3TableCheck_(sof3Tables.getTable<Sof3Table>(satnumIdx));
            }
            if (deck->hasKeyword("SOF2")) {
                sof2TableCheck_(sof2Tables.getTable<Sof2Table>(satnumIdx));
            }
            if (deck->hasKeyword("SGWFN")) {
                sgwfnTableCheck_(sgwfnTables.getTable<SgwfnTable>(satnumIdx));
            }
        }
    }

    void RelpermDiagnostics::swofTableCheck_(const Opm::SwofTable& swofTables)
    {
        const auto& sw = swofTables.getSwColumn();
        const auto& krw = swofTables.getKrwColumn();
        const auto& krow = swofTables.getKrowColumn();
        ///Check sw column.
        if (sw.front()< 0 || sw.back() > 1) {
            std::string s = "In SWOF table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }
        ///TODO check endpoint sw.back() == 1. - Sor.
        ///Check krw column.
        if (krw.front() != 0) {
            std::string s = "In SWOF table, first value of krw should be 0.";
            messager_.push_back(s);
        }
        if (krw.front() < 0 || krw.back() > 1) {
            std::string s = "In SWOF table, krw should be in range [0,1].";
            messager_.push_back(s);
        }

        ///Check krow column.
        if (krow.front() > 1 || krow.back() < 0) {
            std::string s = "In SWOF table, krow should be in range [0, 1].";
            messager_.push_back(s);
        }
        ///TODO check if run with gas.
    }


    void RelpermDiagnostics::sgofTableCheck_(const Opm::SgofTable& sgofTables)
    {
        const auto& sg = sgofTables.getSgColumn();
        const auto& krg = sgofTables.getKrgColumn();
        const auto& krog = sgofTables.getKrogColumn();
                ///Check sw column.
        if (sg.front()< 0 || sg.back() > 1) {
            std::string s = "In SGOF table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }
        if (sg.front() != 0) {
            std::string s = "In SGOF table, first value of sg should be 0.";
            messager_.push_back(s);
        }
        ///TODO check endpoint sw.back() == 1. - Sor.
        ///Check krw column.
        if (krg.front() != 0) {
            std::string s = "In SGOF table, first value of krg should be 0.";
            messager_.push_back(s);
        }
        if (krg.front() < 0 || krg.back() > 1) {
            std::string s = "In SGOF table, krg should be in range [0,1].";
            messager_.push_back(s);
        }

        ///Check krow column.
        if (krog.front() > 1 || krog.back() < 0) {
            std::string s = "In SGOF table, krog should be in range [0, 1].";
            messager_.push_back(s);
        }
        ///TODO check if run with water.
    }

    void RelpermDiagnostics::slgofTableCheck_(const Opm::SlgofTable& slgofTables) 
    {
        const auto& sl = slgofTables.getSlColumn();
        const auto& krg = slgofTables.getKrgColumn();
        const auto& krog = slgofTables.getKrogColumn();

        ///Check sl column.
        ///TODO first value means sl = swco + sor
        if (sl.front()< 0 || sl.back() > 1) {
            std::string s = "In SLGOF table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }
        if (sl.back() != 1) {
            std::string s = "In SLGOF table, last value of sl should be 1.";
            messager_.push_back(s);
        }

        if (krg.front() > 1 || krg.back() < 0) {
            std::string s = "In SLGOF table, krg shoule be in range [0, 1].";
            messager_.push_back(s);
        }
        if (krg.back() != 0) {
            std::string s = "In SLGOF table, last value of krg hould be 0.";
            messager_.push_back(s);
        }

        if (krog.front() < 0 || krog.back() > 1) {
            std::string s = "In SLGOF table, krog shoule be in range [0, 1].";
            messager_.push_back(s);
        }
    }

    void RelpermDiagnostics::swfnTableCheck_(const Opm::SwfnTable& swfnTables)
    {
        const auto& sw = swfnTables.getSwColumn();
        const auto& krw = swfnTables.getKrwColumn();
        
        ///Check sw column.
        if (sw.front() < 0 || sw.back() > 1) {
            std::string s = "In SWFN table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }
        
        ///Check krw column.
        if (krw.front() < 0 || krw.back() > 1) {
            std::string s = "In SWFN table, krw should be in range [0,1].";
            messager_.push_back(s);
        }

        if (krw.front() != 0) {
            std::string s = "In SWFN table, first value of krw should be 0.";
            messager_.push_back(s);
        }
    }

    
    void RelpermDiagnostics::sgfnTableCheck_(const Opm::SgfnTable& sgfnTables)
    {
        const auto& sg = sgfnTables.getSgColumn();
        const auto& krg = sgfnTables.getKrgColumn();
        
        ///Check sg column.
        if (sg.front() < 0 || sg.back() > 1) {
            std::string s = "In SGFN table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }
        
        ///Check krg column.
        if (krg.front() < 0 || krg.back() > 1) {
            std::string s = "In SGFN table, krg should be in range [0,1].";
            messager_.push_back(s);
        }
        if (krg.front() != 0) {
            std::string s = "In SGFN table, first value of krg should be 0.";
            messager_.push_back(s);
        }
    }

    void RelpermDiagnostics::sof3TableCheck_(const Opm::Sof3Table& sof3Tables)
    {
        const auto& so = sof3Tables.getSoColumn();
        const auto& krow = sof3Tables.getKrowColumn();
        const auto& krog = sof3Tables.getKrogColumn();

        ///Check so column.
        ///TODO: The max so = 1 - Swco
        if (so.front() < 0 || so.back() > 1) {
            std::string s = "In SOF3 table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }

        ///Check krow column.
        if (krow.front() < 0 || krow.back() > 1) {
            std::string s = "In SOF3 table, krow should be in range [0,1].";
            messager_.push_back(s);
        }
        if (krow.front() != 0) {
            std::string s = "In SOF3 table, first value of krow should be 0.";
            messager_.push_back(s);
        }

        ///Check krog column.
        if (krog.front() < 0 || krog.back() > 1) {
            std::string s = "In SOF3 table, krog should be in range [0,1].";
            messager_.push_back(s);
        }

        if (krog.front() != 0) {
            std::string s = "In SOF3 table, first value of krog should be 0.";
            messager_.push_back(s);
        }
    
        if (krog.back() != krow.back()) {
            std::string s = "In SOF3 table, max value of krog and krow should be the same.";
            messager_.push_back(s);
        }
    }


    void RelpermDiagnostics::sof2TableCheck_(const Opm::Sof2Table& sof2Tables)
    {
        const auto& so = sof2Tables.getSoColumn();
        const auto& kro = sof2Tables.getKroColumn();

        ///Check so column.
        ///TODO: The max so = 1 - Swco
        if (so.front() < 0 || so.back() > 1) {
            std::string s = "In SOF2 table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }

        ///Check krow column.
        if (kro.front() < 0 || kro.back() > 1) {
            std::string s = "In SOF2 table, krow should be in range [0,1].";
            messager_.push_back(s);
        }
        if (kro.front() != 0) {
            std::string s = "In SOF2 table, first value of krow should be 0.";
            messager_.push_back(s);
        }
    }


    void RelpermDiagnostics::sgwfnTableCheck_(const Opm::SgwfnTable& sgwfnTables)
    {
        const auto& sg = sgwfnTables.getSgColumn();
        const auto& krg = sgwfnTables.getKrgColumn();
        const auto& krgw = sgwfnTables.getKrgwColumn();
        
        ///Check sg column.
        if (sg.front() < 0 || sg.back() > 1) {
            std::string s = "In SGWFN table, saturation should be in range [0,1].";
            messager_.push_back(s);
        }

        ///Check krg column.
        if (krg.front() < 0 || krg.back() > 1) {
            std::string s = "In SGWFN table, krg should be in range [0,1].";
            messager_.push_back(s);
        }
        if (krg.front() != 0) {
            std::string s = "In SGWFN table, first value of krg should be 0.";
            messager_.push_back(s);
        }

        ///Check krgw column.
        ///TODO check saturation sw = 1. - sg
        if (krgw.front() > 1 || krgw.back() < 0) {
            std::string s = "In SGWFN table, krgw should be in range [0,1].";
            messager_.push_back(s);
        }
        if (krgw.back() != 0) {
            std::string s = "In SGWFN table, last value of krgw should be 0.";
            messager_.push_back(s);
        }
    }




    void RelpermDiagnostics::unscaledEndPointsCheck_(DeckConstPtr deck,
                                                     EclipseStateConstPtr eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));
        unscaledEpsInfo_.resize(numSatRegions);

        auto tables = eclState->getTableManager();
        const TableContainer&  swofTables = tables->getSwofTables();
        const TableContainer&  sgofTables = tables->getSgofTables();
        const TableContainer& slgofTables = tables->getSlgofTables();
        const TableContainer&  sof3Tables = tables->getSof3Tables();

        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
             unscaledEpsInfo_[satnumIdx].extractUnscaled(deck, eclState, satnumIdx);
             std::cout << "***************\nEnd-Points In all the Tables\n";
             unscaledEpsInfo_[satnumIdx].print();
             ///Consistency check.
             if (unscaledEpsInfo_[satnumIdx].Sgu > (1. - unscaledEpsInfo_[satnumIdx].Swl)) {
                messager_.push_back("Sgmax should not exceed 1-Swco.");
             }
             if (unscaledEpsInfo_[satnumIdx].Sgl > (1. - unscaledEpsInfo_[satnumIdx].Swu)) {
                messager_.push_back("Sgco should not exceed 1-Swmax.");
             }

             ///Krow(Sou) == Krog(Sou) for three-phase
             /// means Krow(Swco) == Krog(Sgco)
             double krow_value = 1e20;
             double krog_value = 1e-20;
             if (fluidSystem_ = FluidSystem::BlackOil) {
                if (satFamily_ == SaturationFunctionFamily::FamilyI) {
                     if (!sgofTables.empty()) {
                         auto sg = sgofTables.getTable<SgofTable>(satnumIdx).getSgColumn();
                         auto krog = sgofTables.getTable<SgofTable>(satnumIdx).getKrogColumn();
                         krog_value=Opm::linearInterpolation(sg, krog,unscaledEpsInfo_[satnumIdx].Sgl);
                     } else {
                         assert(!slgofTables.empty());
                         auto sl = slgofTables.getTable<SlgofTable>(satnumIdx).getSlColumn();
                         auto krog = slgofTables.getTable<SlgofTable>(satnumIdx).getKrogColumn();
                         krog_value=Opm::linearInterpolation(sl, krog, unscaledEpsInfo_[satnumIdx].Sgl);                        
                     }
                     auto sw = swofTables.getTable<SwofTable>(satnumIdx).getSwColumn();
                     auto krow = swofTables.getTable<SwofTable>(satnumIdx).getKrowColumn();
                     krow_value = Opm::linearInterpolation(sw, krow,unscaledEpsInfo_[satnumIdx].Swl);
                 }
                 if (satFamily_ == SaturationFunctionFamily::FamilyII) {
                     assert(!sof3Table.empty());
                     const double Sou = 1.- unscaledEpsInfo_[satnumIdx].Swl - unscaledEpsInfo_[satnumIdx].Sgl;
                     auto so = sof3Tables.getTable<Sof3Table>(satnumIdx).getSoColumn();
                     auto krow = sof3Tables.getTable<Sof3Table>(satnumIdx).getKrowColumn();
                     auto krog = sof3Tables.getTable<Sof3Table>(satnumIdx).getKrogColumn();
                     krow_value = Opm::linearInterpolation(so, krow, Sou);
                     krog_value = Opm::linearInterpolation(so, krog, Sou);
                 }
                 if (krow_value != krog_value) {
                     messager_.push_back("Krow(sSomax) should equal Krog(Somax).");
                 }
             }
             ///Krw(Sw=0)=Krg(Sg=0)=Krow(So=0)=Krog(So=0)=0.
             ///Mobile fluid requirements
            if (((unscaledEpsInfo_[satnumIdx].Sowcr + unscaledEpsInfo_[satnumIdx].Swcr)-1) >= 0) {
                messager_.push_back("Sowcr + Swcr should less than 1.");
            }
            if (((unscaledEpsInfo_[satnumIdx].Sogcr + unscaledEpsInfo_[satnumIdx].Sgcr + unscaledEpsInfo_[satnumIdx].Swl) - 1 ) > 0) {
                messager_.push_back("Sogcr + Sgcr + Swco should less than 1.");
            }
        }
    }



    void RelpermDiagnostics::scaledEndPointsCheck_(DeckConstPtr deck,
                                                   EclipseStateConstPtr eclState,
                                                   const UnstructuredGrid& grid)
    {
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));

        const int nc = Opm::UgGridHelpers::numCells(grid);
        std::vector<int> compressedToCartesianIdx(nc);
        const auto& global_cell = Opm::UgGridHelpers::globalCell(grid);
        for (int cell = 0; cell < nc; ++cell) {
            if (global_cell) {
                compressedToCartesianIdx[cell] = global_cell[cell];
            } else {
                compressedToCartesianIdx[cell] = cell;
            }
        }
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);       

        for (int c = 0; c < nc; ++c) {
            int cartIdx = compressedToCartesianIdx[c];
            scaledEpsInfo_[c].extractScaled(epsGridProperties, cartIdx);
            //std::cout << "Scaled endpoints info:\n";
            //scaledEpsInfo_[c].print();
            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl)) {
                std::string msg = "In cell:" + std::to_string(c) + "SGU exceed 1.0 - SWL";
                messager_.push_back(msg);
            }
            
            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu)) {
                std::string msg = "In cell: " + std::to_string(c) + "SGL exceed 1.0 - SWU";
                messager_.push_back(msg);
            }
        
            // Mobilility check.
            if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= 1.0) {
                std::string msg = "In cell: " + std::to_string(c) + "SOWCR + SWCR exceed 1.0";
                messager_.push_back(msg);
            }

            if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= 1.0) {
                std::string msg = "In cell: " + std::to_string(c) + "SOGCR + SGCR + SWL";
                messager_.push_back(msg);
            }
        } 
    }

} //namespace Opm

#endif // OPM_RELPERMDIAGNOSTICS_HEADER_INCLUDED
