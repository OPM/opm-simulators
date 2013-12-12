/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/polymer/PolymerInflow.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/wells.h>
#include <map>

namespace Opm
{

    // ---------- Methods of PolymerInflowBasic ----------

    /// Constructor.
    /// @param[in]  starttime  Start time of injection in seconds.
    /// @param[in]  endtime    End time of injection in seconds.
    /// @param[in]  amount     Amount to be injected per second.
    PolymerInflowBasic::PolymerInflowBasic(const double starttime,
                                           const double endtime,
                                           const double amount)
        : stime_(starttime), etime_(endtime), amount_(amount)
    {
    }

    void PolymerInflowBasic::getInflowValues(const double step_start,
                                             const double step_end,
                                             std::vector<double>& poly_inflow_c) const
    {
        const double eps = 1e-5*(step_end - step_start);
        if (step_start + eps >= stime_ && step_end - eps <= etime_) {
            std::fill(poly_inflow_c.begin(), poly_inflow_c.end(), amount_);
        } else if (step_start + eps <= etime_ && step_end - eps >= stime_) {
            OPM_MESSAGE("Warning: polymer injection set to change inside timestep. Using value at start of step.");
            std::fill(poly_inflow_c.begin(), poly_inflow_c.end(), amount_);
        } else {
            std::fill(poly_inflow_c.begin(), poly_inflow_c.end(), 0.0);
        }
    }


    // ---------- Methods of PolymerInflowFromDeck ----------



    /// Constructor.
    /// @param[in]  deck     Input deck expected to contain WPOLYMER.
    PolymerInflowFromDeck::PolymerInflowFromDeck(const EclipseGridParser& deck,
                                                 const Wells& wells,
                                                 const int num_cells)
        : sparse_inflow_(num_cells)
    {
        if (!deck.hasField("WPOLYMER")) {
            OPM_MESSAGE("PolymerInflowFromDeck initialized without WPOLYMER in current epoch.");
            return;
        }

        // Extract concentrations and put into cell->concentration map.
        const std::vector<WpolymerLine>& wpl = deck.getWPOLYMER().wpolymer_;
        const int num_wpl = wpl.size();
        std::map<int, double> perfcell_conc;
        for (int i = 0; i < num_wpl; ++i) {
            // Only use well name and polymer concentration.
            // That is, we ignore salt concentration and group
            // names.
            int wix = 0;
            for (; wix < wells.number_of_wells; ++wix) {
                if (wpl[i].well_ == wells.name[wix]) {
                    break;
                }
            }
            if (wix == wells.number_of_wells) {
                OPM_THROW(std::runtime_error, "Could not find a match for well " << wpl[i].well_ << " from WPOLYMER.");
            }
            for (int j = wells.well_connpos[wix]; j < wells.well_connpos[wix+1]; ++j) {
                const int perf_cell = wells.well_cells[j];
                perfcell_conc[perf_cell] = wpl[i].polymer_concentration_;
            }
        }

        // Build sparse vector from map.
        std::map<int, double>::const_iterator it = perfcell_conc.begin();
        for (; it != perfcell_conc.end(); ++it) {
            sparse_inflow_.addElement(it->second, it->first);
        }
    }


    void PolymerInflowFromDeck::getInflowValues(const double /*step_start*/,
                                                const double /*step_end*/,
                                                std::vector<double>& poly_inflow_c) const
    {
        // This method does not depend on the given time,
        // instead one would have a new epoch (and create a new
        // instance) for each change in WPOLYMER.
        std::fill(poly_inflow_c.begin(), poly_inflow_c.end(), 0.0);
        const int nnz = sparse_inflow_.nonzeroSize();
        for (int i = 0; i < nnz; ++i) {
            poly_inflow_c[sparse_inflow_.nonzeroIndex(i)] = sparse_inflow_.nonzeroElement(i) ;
        }
    }


} // namespace Opm
