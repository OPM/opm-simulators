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

namespace Opm
{

    /// \brief Make a vector of pressure thresholds from deck input.
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
    std::vector<double> thresholdPressures(DeckConstPtr deck, EclipseStateConstPtr eclipseState, const Grid& grid)
    {
        // Is the EQLOPTS option THPRES activated?
        bool has_thpres_option = false;
        if (deck->hasKeyword("EQLOPTS")) {
            auto eqlopts = deck->getKeyword("EQLOPTS");
            auto rec = eqlopts->getRecord(0);
            for (size_t i = 0; i < rec->size(); ++i) {
                if (rec->getItem(i)->getString(0) == "THPRES") {
                    has_thpres_option = true;
                } else if (rec->getItem(i)->getString(0) == "IRREVERS") {
                    OPM_THROW(std::runtime_error, "Cannot use IRREVERS version of THPRES option, not implemented.");
                }
            }
        }

        // Do we have the THPRES keyword?
        // Check for consistency.
        const bool has_thpres_keyword = deck->hasKeyword("THPRES");
        if (has_thpres_keyword != has_thpres_option) {
            OPM_THROW(std::runtime_error, "Invalid deck, the THPRES keyword must be present and the THPRES "
                      "option of EQLOPTS must be used for the threshold pressure feature.");
        }

        if (has_thpres_option) {
            // Check for EQLNUM
            if (!eclipseState->hasIntGridProperty("EQLNUM")) {
                OPM_THROW(std::runtime_error, "Invalid deck, EQLNUM must be set in order to use THPRES.");
            }

            // Create threshold pressure table.
            auto thpres = deck->getKeyword("THPRES");
            auto eqlnum = eclipseState->getIntGridProperty("EQLNUM")->getData();
            const int max_eqlnum = *max_element(eqlnum.begin(), eqlnum.end());
            std::vector<double> thp_table(max_eqlnum * max_eqlnum, 0.0);
            const int num_records = thpres->size();
            for (int rec_ix = 0; rec_ix < num_records; ++rec_ix) {
                auto rec = thpres->getRecord(rec_ix);
                const int r1 = rec->getItem("REGION1")->getInt(0) - 1;
                const int r2 = rec->getItem("REGION2")->getInt(0) - 1;
                const double p = rec->getItem("THPRES")->getSIDouble(0);
                if (r1 >= max_eqlnum || r2 >= max_eqlnum) {
                    OPM_THROW(std::runtime_error, "Too high region numbers in THPRES keyword.");
                }
                thp_table[r1 + max_eqlnum*r2] = p;
                thp_table[r2 + max_eqlnum*r1] = p;
            }

            // Set values for each face.
            const int num_faces = UgGridHelpers::numFaces(grid);
            std::vector<double> thpres_vals(num_faces, 0.0);
            const int* gc = UgGridHelpers::globalCell(grid);
            auto fc = UgGridHelpers::faceCells(grid);
            for (int face = 0; face < num_faces; ++face) {
                const int c1 = fc(face, 0);
                const int c2 = fc(face, 1);
                if (c1 < 0 || c2 < 0) {
                    // Boundary face, skip this.
                    continue;
                }
                const int eq1 = eqlnum[gc[c1]] - 1;
                const int eq2 = eqlnum[gc[c2]] - 1;
                thpres_vals[face] = thp_table[eq1 + max_eqlnum*eq2];
            }
            return thpres_vals;
        } else {
            return std::vector<double>();
        }
    }

}

#endif // OPM_THRESHOLDPRESSURES_HEADER_INCLUDED
