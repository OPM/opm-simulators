/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_THRESHOLDPRESSURES_HEADER_INCLUDED
#define OPM_THRESHOLDPRESSURES_HEADER_INCLUDED

#include <vector>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>


namespace Opm
{

    /// \brief Get a vector of pressure thresholds from EclipseState.
    /// This function looks at EQLOPTS, THPRES and EQLNUM to determine
    /// pressure thresholds.  It does not consider the case where the
    /// threshold values are defaulted, in which case they should be
    /// determined from the initial, equilibrated simulation state.
    /// \tparam    Grid           Type of grid object (UnstructuredGrid or CpGrid).
    /// \param[in] deck           Input deck, EQLOPTS and THPRES are accessed from it.
    /// \param[in] eclipseState   Processed eclipse state, EQLNUM is accessed from it.
    /// \param[in] grid           The grid to which the thresholds apply.
    /// \return                   A vector of pressure thresholds, one
    ///                           for each face in the grid. A value
    ///                           of zero means no threshold for that
    ///                           particular face. An empty vector is
    ///                           returned if there is no THPRES
    ///                           feature used in the deck.



    template <class Grid>
    std::vector<double> thresholdPressures(const ParseMode& parseMode ,EclipseStateConstPtr eclipseState, const Grid& grid)
    {
        SimulationConfigConstPtr simulationConfig = eclipseState->getSimulationConfig();
        std::vector<double> thpres_vals;
        if (simulationConfig->hasThresholdPressure()) {
            std::shared_ptr<const ThresholdPressure> thresholdPressure = simulationConfig->getThresholdPressure();
            std::shared_ptr<GridProperty<int>> eqlnum = eclipseState->getIntGridProperty("EQLNUM");
            auto eqlnumData = eqlnum->getData();

            // Set values for each face.
            const int num_faces = UgGridHelpers::numFaces(grid);
            thpres_vals.resize(num_faces, 0.0);
            const int* gc = UgGridHelpers::globalCell(grid);
            auto fc = UgGridHelpers::faceCells(grid);
            for (int face = 0; face < num_faces; ++face) {
                const int c1 = fc(face, 0);
                const int c2 = fc(face, 1);
                if (c1 < 0 || c2 < 0) {
                    // Boundary face, skip this.
                    continue;
                }
                const int gc1 = (gc == 0) ? c1 : gc[c1];
                const int gc2 = (gc == 0) ? c2 : gc[c2];
                const int eq1 = eqlnumData[gc1];
                const int eq2 = eqlnumData[gc2];

                if (thresholdPressure->hasRegionBarrier(eq1,eq2)) {
                    if (thresholdPressure->hasThresholdPressure(eq1,eq2)) {
                        thpres_vals[face] = thresholdPressure->getThresholdPressure(eq1,eq2);
                    } else {
                        std::string msg = "Initializing the THPRES pressure values from the initial state is not supported - using 0.0";
                        parseMode.handleError( ParseMode::UNSUPPORTED_INITIAL_THPRES , msg );
                    }
                }
            }
        }
        return thpres_vals;
    }
}

#endif // OPM_THRESHOLDPRESSURES_HEADER_INCLUDED
