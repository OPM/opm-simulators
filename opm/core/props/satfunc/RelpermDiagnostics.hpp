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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/OpmLog/StreamLog.hpp>
#include <opm/parser/eclipse/OpmLog/OpmLog.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

namespace Opm {

    class Sof2Table;
    class SgwfnTable;

    ///This class is intend to be a relpmer diganostics, to detect
    ///wrong input of relperm table and endpoints.
    class RelpermDiagnostics 
    {
    public:

        ///Constructor for OpmLog.
        explicit RelpermDiagnostics(std::string& logFile);

        ///This function is used to diagnosis relperm in
        ///eclipse data file. Errors and warings will be 
        ///output if they're found.
        ///\param[in] eclState  eclipse state.
        ///\param[in] deck      ecliplise data file.
        ///\param[in] grid      unstructured grid.
        template <class GridT>
        void diagnosis(EclipseStateConstPtr eclState,
                       DeckConstPtr deck,
                       const GridT& grid);

        ///return streamLog
        std::shared_ptr<Opm::StreamLog> getOpmLog() const;
        std::vector<std::string> getMessages() const;

    private:
        enum FluidSystem {
            OilWater,
            OilGas,
            WaterGas,
            BlackOil,
            Solvent
        };
        
        FluidSystem fluidSystem_;

        enum SaturationFunctionFamily {
            FamilyI,
            FamilyII,
            NoFamily
        };
  
        SaturationFunctionFamily satFamily_;

        struct Counter {
            Counter();
            int error;
            int warning;
            int problem;
            int bug;
        };

        Counter counter_;
        
        std::vector<Opm::EclEpsScalingPointsInfo<double> > unscaledEpsInfo_;
        std::vector<Opm::EclEpsScalingPointsInfo<double> > scaledEpsInfo_;

        std::vector<std::string> messages_;
        ///Store scaled information.
        std::vector<std::string> scaled_messages_;

        ///Use OpmLog
        std::shared_ptr<Opm::StreamLog> streamLog_;

        ///Check the phase that used.
        void phaseCheck_(DeckConstPtr deck);

        ///Check saturation family I and II.
        void satFamilyCheck_(EclipseStateConstPtr eclState);
 
        ///Check saturation tables.
        void tableCheck_(EclipseStateConstPtr eclState, 
                         DeckConstPtr deck);

        ///Check endpoints in the saturation tables.
        void unscaledEndPointsCheck_(DeckConstPtr deck,
                                     EclipseStateConstPtr eclState);

        template <class GridT>
        void scaledEndPointsCheck_(DeckConstPtr deck,
                                   EclipseStateConstPtr eclState,
                                   const GridT& grid);

        ///For every table, need to deal with case by case.
        void swofTableCheck_(const Opm::SwofTable& swofTables,
                             const int satnumIdx);
        void sgofTableCheck_(const Opm::SgofTable& sgofTables,
                             const int satnumIdx);
        void slgofTableCheck_(const Opm::SlgofTable& slgofTables,
                              const int satnumIdx);
        void swfnTableCheck_(const Opm::SwfnTable& swfnTables,
                             const int satnumIdx);
        void sgfnTableCheck_(const Opm::SgfnTable& sgfnTables,
                             const int satnumIdx);
        void sof3TableCheck_(const Opm::Sof3Table& sof3Tables,
                             const int satnumIdx);
        void sof2TableCheck_(const Opm::Sof2Table& sof2Tables,
                             const int satnumIdx);
        void sgwfnTableCheck_(const Opm::SgwfnTable& sgwfnTables,
                              const int satnumIdx);
    };

} //namespace Opm

#include <opm/core/props/satfunc/RelpermDiagnostics_impl.hpp>

#endif // OPM_RELPERMDIAGNOSTICS_HEADER_INCLUDED
