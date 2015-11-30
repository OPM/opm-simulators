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

#include "config.h"
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

} //namespace Opm

#endif // OPM_RELPERMDIAGNOSTICS_HEADER_INCLUDED
