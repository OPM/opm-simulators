/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014-2016 IRIS AS
  Copyright 2015 Andreas Lauser

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
#ifndef OPM_SIM_FIBO_DETAILS_HPP
#define OPM_SIM_FIBO_DETAILS_HPP

#include <utility>
#include <algorithm>
#include <locale>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/core/well_controls.h>

namespace Opm
{
    namespace SimFIBODetails {
        typedef std::unordered_map<std::string, const Well* > WellMap;

        inline WellMap
        mapWells(const std::vector< const Well* >& wells)
        {
            WellMap wmap;

            for (std::vector< const Well* >::const_iterator
                     w = wells.begin(), e = wells.end();
                 w != e; ++w)
            {
                wmap.insert(std::make_pair((*w)->name(), *w));
            }

            return wmap;
        }

        inline int
        resv_control(const WellControls* ctrl)
        {
            int i, n = well_controls_get_num(ctrl);

            bool match = false;
            for (i = 0; (! match) && (i < n); ++i) {
                match = well_controls_iget_type(ctrl, i) == RESERVOIR_RATE;
            }

            if (! match) { i = 0; }

            return i - 1; // -1 if no match, undo final "++" otherwise
        }

        inline bool
        is_resv(const Wells& wells,
                const int    w)
        {
            return (0 <= resv_control(wells.ctrls[w]));
        }

        inline std::vector<int>
        resvWells(const Wells*      wells)
        {
            std::vector<int> resv_wells;
            if( wells )
            {
                for (int w = 0, nw = wells->number_of_wells; w < nw; ++w) {
                    if ( is_resv(*wells, w) ) {
                        resv_wells.push_back(w);
                    }
                }
            }

            return resv_wells;
        }

        inline void
        historyRates(const PhaseUsage&               pu,
                     const WellProductionProperties& p,
                     std::vector<double>&            rates)
        {
            assert (! p.predictionMode);
            assert (rates.size() ==
                    std::vector<double>::size_type(pu.num_phases));

            if (pu.phase_used[ BlackoilPhases::Aqua ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Aqua ];

                rates[i] = p.WaterRate;
            }

            if (pu.phase_used[ BlackoilPhases::Liquid ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Liquid ];

                rates[i] = p.OilRate;
            }

            if (pu.phase_used[ BlackoilPhases::Vapour ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Vapour ];

                rates[i] = p.GasRate;
            }
        }
    } // namespace SimFIBODetails
} // namespace Opm

#endif // OPM_SIM_FIBO_DETAILS_HPP
