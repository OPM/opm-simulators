/*
  Copyright 2016 Statoil ASA.

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

#ifndef OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED
#define OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED

#include <vector>
#include <utility>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/core/utility/compressedToCartesian.hpp>

namespace Opm {

    template <class GridT>
    void RelpermDiagnostics::diagnosis(Opm::EclipseStateConstPtr eclState,
                                       Opm::DeckConstPtr deck,
                                       const GridT& grid)
    {
        std::cout << "\n\n***************Saturation Functions Diagnostics***************\n\n";
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
                std::cout << x << std::endl;
            }
        }
        int limits = 0;
        if (!scaled_messages_.empty()) {
            std::cout << std::endl;
            for (const auto& x : scaled_messages_) {
                if (limits < 10) {
                    std::cout << x << std::endl;
                    limits++;
                } else {
                    std::cout << "\nMore inconsistencies exist. Check saturation function input and LOGFILE!" << std::endl;
                    break;
                }
            }
        }

        const std::string summary_msg = "\n\nError summary:" + 
            std::string("\nWarnings          " + std::to_string(counter_.warning)) +
            std::string("\nProblems          " + std::to_string(counter_.problem)) +
            std::string("\nErrors            " + std::to_string(counter_.error)) + 
            std::string("\nBugs              " + std::to_string(counter_.bug))+ "\n";
        std::cout << summary_msg << std::endl;
    }

    template <class GridT>
    void RelpermDiagnostics::scaledEndPointsCheck_(DeckConstPtr deck,
                                                   EclipseStateConstPtr eclState,
                                                   const GridT& grid)
    {
        const int nc = Opm::UgGridHelpers::numCells(grid);
        const auto& global_cell = Opm::UgGridHelpers::globalCell(grid);
        const auto dims = Opm::UgGridHelpers::cartDims(grid);
        const auto& compressedToCartesianIdx = Opm::compressedToCartesian(nc, global_cell);
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);       
        const auto satnum = eclState->getIntGridProperty("SATNUM");

        for (int c = 0; c < nc; ++c) {
            const int cartIdx = compressedToCartesianIdx[c];
            const std::string satnumIdx = std::to_string(satnum->iget(cartIdx));
            std::array<int, 3> ijk;
            ijk[0] = cartIdx % dims[0];
            ijk[1] = (cartIdx / dims[0]) % dims[1];
            ijk[2] = cartIdx / dims[0] / dims[1];
            const std::string cellIdx = "(" + std::to_string(ijk[0]) + ", " + 
                                   std::to_string(ijk[1]) + ", " +
                                   std::to_string(ijk[2]) + ")";
            scaledEpsInfo_[c].extractScaled(epsGridProperties, cartIdx);

            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl)) {
                const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGU exceed 1.0 - SWL";
                scaled_messages_.push_back(msg);
                OpmLog::warning(msg);
                counter_.warning += 1;
            }
            
            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu)) {
                const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGL exceed 1.0 - SWU";
                scaled_messages_.push_back(msg);
                OpmLog::warning(msg);
                counter_.warning += 1;
            }

            if (deck->hasKeyword("SCALECRS") && fluidSystem_ == FluidSystem::BlackOil) {
                // Mobilility check.
                if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= 1.0) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOWCR + SWCR exceed 1.0";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }

                if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= 1.0) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOGCR + SGCR + SWL exceed 1.0";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.error += 1;
                }
            }
            ///Following rules come from NEXUS.
            if (fluidSystem_ != FluidSystem::WaterGas) {
                if (scaledEpsInfo_[c].Swl > scaledEpsInfo_[c].Swcr) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SWL > SWCR";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }

                if (scaledEpsInfo_[c].Swcr > scaledEpsInfo_[c].Sowcr) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SWCR > SOWCR";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }
            
                if (scaledEpsInfo_[c].Sowcr > scaledEpsInfo_[c].Swu) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOWCR > SWU";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }
            }

            if (fluidSystem_ != FluidSystem::OilWater) {
                if (scaledEpsInfo_[c].Sgl > scaledEpsInfo_[c].Sgcr) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGL > SGCR";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }
            }

            if (fluidSystem_ != FluidSystem::BlackOil) {
                if (scaledEpsInfo_[c].Sgcr > scaledEpsInfo_[c].Sogcr) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGCR > SOGCR";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }

                if (scaledEpsInfo_[c].Sogcr > scaledEpsInfo_[c].Sgu) {
                    const std::string msg = "-- Warning: For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOGCR > SGU";
                    scaled_messages_.push_back(msg);
                    OpmLog::warning(msg);
                    counter_.warning += 1;
                }
            }
        } 
    }

} //namespace Opm

#endif // OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED
