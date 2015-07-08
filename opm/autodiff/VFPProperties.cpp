/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#include "config.h"

#include <opm/autodiff/VFPProperties.hpp>

#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <algorithm>

namespace Opm {


VFPProperties::VFPProperties() {
}

VFPProperties::VFPProperties(const VFPInjTable* inj_table, const VFPProdTable* prod_table) {
    if (inj_table != NULL) {
        m_inj_tables[inj_table->getTableNum()] = inj_table;
    }
    if (prod_table != NULL) {
        m_prod_tables[prod_table->getTableNum()] = prod_table;
    }
}

VFPProperties::VFPProperties(const std::map<int, VFPInjTable>& inj_tables,
                             const std::map<int, VFPProdTable>& prod_tables) {
    init(inj_tables);
    init(prod_tables);
}


VFPProperties::VFPProperties(const std::map<int, VFPInjTable>& inj_tables) {
    init(inj_tables);
}


VFPProperties::VFPProperties(const std::map<int, VFPProdTable>& prod_tables) {
    init(prod_tables);
}

void VFPProperties::init(const std::map<int, VFPInjTable>& inj_tables) {
    //Populate injection table pointers.
    for (const auto& table : inj_tables) {
        m_inj_tables[table.first] = &table.second;
    }
}

void VFPProperties::init(const std::map<int, VFPProdTable>& prod_tables) {
    //Populate production table pointers
    for (const auto& table : prod_tables) {
        m_prod_tables[table.first] = &table.second;
    }
}

VFPProperties::ADB::V VFPProperties::prod_bhp(int table_id,
        const Wells& wells,
        const ADB::V& qs,
        const ADB::V& thp,
        const ADB::V& alq) const {
    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(np == 3);
    const ADB::V& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB::V& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB::V& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    return prod_bhp(table_id, w, o, g, thp, alq);
}

VFPProperties::ADB::V VFPProperties::prod_bhp(int table_id,
        const ADB::V& aqua,
        const ADB::V& liquid,
        const ADB::V& vapour,
        const ADB::V& thp,
        const ADB::V& alq) const {
    const int nw = thp.size();

    assert(aqua.size()   == nw);
    assert(liquid.size() == nw);
    assert(vapour.size() == nw);
    assert(thp.size()    == nw);
    assert(alq.size()    == nw);

    //Compute the BHP for each well independently
    ADB::V bhp_vals;
    bhp_vals.resize(nw);
    for (int i=0; i<nw; ++i) {
        bhp_vals[i] = prod_bhp(table_id, aqua[i], liquid[i], vapour[i], thp[i], alq[i]);
    }
    return bhp_vals;
}

double VFPProperties::prod_bhp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp,
        const double& alq) const {
    const VFPProdTable* table = getProdTable(table_id);

    //Find interpolation variables
    double flo = getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = getWFR(aqua, liquid, vapour, table->getWFRType());
    double gfr = getGFR(aqua, liquid, vapour, table->getGFRType());

    //First, find the values to interpolate between
    auto flo_i = find_interp_data(flo, table->getFloAxis());
    auto thp_i = find_interp_data(thp, table->getTHPAxis());
    auto wfr_i = find_interp_data(wfr, table->getWFRAxis());
    auto gfr_i = find_interp_data(gfr, table->getGFRAxis());
    auto alq_i = find_interp_data(alq, table->getALQAxis());

    //Then perform the interpolation itself
    return interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);
}


double VFPProperties::prod_thp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& bhp,
        const double& alq) const {
    const VFPProdTable* table = getProdTable(table_id);
    const VFPProdTable::array_type& data = table->getTable();

    double thp = -1e100;

    //Find interpolation variables
    double flo = getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = getWFR(aqua, liquid, vapour, table->getWFRType());
    double gfr = getGFR(aqua, liquid, vapour, table->getGFRType());

    /**
     * Get THP axis, assume that it is sorted
     */
    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();
    assert(std::is_sorted(thp_array.begin(), thp_array.end()));

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    auto flo_i = find_interp_data(flo, table->getFloAxis());
    auto wfr_i = find_interp_data(wfr, table->getWFRAxis());
    auto gfr_i = find_interp_data(gfr, table->getGFRAxis());
    auto alq_i = find_interp_data(alq, table->getALQAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = find_interp_data(thp_array[i], thp_array);
        bhp_array[i] = interpolate(data, flo_i, thp_i, wfr_i, gfr_i, alq_i);
    }

    /**
     * Our *interpolated* bhp_array will be montoic increasing for increasing
     * THP if our input BHP values are monotonic increasing for increasing
     * THP values. However, if we have to *extrapolate* along any of the other
     * axes, this guarantee holds no more, and bhp_array may be "random"
     */
    if (std::is_sorted(bhp_array.begin(), bhp_array.end())) {
        //Target bhp less than all values in array, extrapolate
        if (bhp <= bhp_array[0]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[0];
            const double& x1 = thp_array[1];
            const double& y0 = bhp_array[0];
            const double& y1 = bhp_array[1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
        //Target bhp within table ranges, interpolate
        else {
            //Loop over the values and find min(bhp_array(thp)) == bhp
            //so that we maximize the rate.

            //Find i so that bhp_array[i-1] <= bhp <= bhp_array[i];
            //Assuming a small number of values in bhp_array, this should be quite
            //efficient. Other strategies might be bisection, etc.
            int i=0;
            bool found = false;
            for (; i<nthp-1; ++i) {
                const double& y0 = bhp_array[i  ];
                const double& y1 = bhp_array[i+1];

                if (y0 < bhp && bhp <= y1) {
                    found = true;
                    break;
                }
            }
            //Canary in a coal mine: shouldn't really be required
            assert(found == true);

            const double& x0 = thp_array[i  ];
            const double& x1 = thp_array[i+1];
            const double& y0 = bhp_array[i  ];
            const double& y1 = bhp_array[i+1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
    }
    //bhp_array not sorted, raw search.
    else {
        //Find i so that bhp_array[i-1] <= bhp <= bhp_array[i];
        //Since the BHP values might not be sorted, first search within
        //our interpolation values, and then try to extrapolate.
        int i=0;
        bool found = false;
        for (; i<nthp-1; ++i) {
            const double& y0 = bhp_array[i  ];
            const double& y1 = bhp_array[i+1];

            if (y0 < bhp && bhp <= y1) {
                found = true;
                break;
            }
        }
        if (found) {
            const double& x0 = thp_array[i  ];
            const double& x1 = thp_array[i+1];
            const double& y0 = bhp_array[i  ];
            const double& y1 = bhp_array[i+1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
        else if (bhp <= bhp_array[0]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[0];
            const double& x1 = thp_array[1];
            const double& y0 = bhp_array[0];
            const double& y1 = bhp_array[1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
            thp = find_x(x0, x1, y0, y1, bhp);
        }
        else {
            OPM_THROW(std::logic_error, "Programmer error: Unable to find THP in THP array");
        }
    }

    return thp;
}



const VFPProdTable* VFPProperties::getProdTable(int table_id) const {
    auto entry = m_prod_tables.find(table_id);
    if (entry == m_prod_tables.end()) {
        OPM_THROW(std::invalid_argument, "Nonexistent table " << table_id << " referenced.");
    }
    else {
        return entry->second;
    }
}

VFPProperties::InterpData VFPProperties::find_interp_data(const double& value, const std::vector<double>& values) {
    InterpData retval;

    //First element greater than or equal to value
    //Start with the second element, so that floor_iter does not go out of range
    //Don't access out-of-range, therefore values.end()-1
    auto ceil_iter = std::lower_bound(values.begin()+1, values.end()-1, value);

    //Find last element smaller than range
    auto floor_iter = ceil_iter-1;

    //Find the indices
    const int a = floor_iter - values.begin();
    const int b = ceil_iter - values.begin();
    const int max_size = std::max(static_cast<int>(values.size()) - 1, 0);

    //Clamp indices to range of vector
    retval.ind_[0] = a;
    retval.ind_[1] = std::min(b, max_size);

    //Find interpolation ratio
    double dist = (*ceil_iter - *floor_iter);
    assert(dist >= 0.0);
    if (dist > 0.0) {
        //Possible source for floating point error here if value and floor are large,
        //but very close to each other
        retval.factor_ = (value-*floor_iter) / dist;
    }
    else {
        retval.factor_ = 1.0;
    }

    return retval;
}


#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif

double VFPProperties::interpolate(const VFPProdTable::array_type& array,
        const InterpData& flo_i,
        const InterpData& thp_i,
        const InterpData& wfr_i,
        const InterpData& gfr_i,
        const InterpData& alq_i) {
    double nn[2][2][2][2][2];

    //Pick out nearest neighbors (nn) to our evaluation point
    //This is not really required, but performance-wise it may pay off, since the 32-elements
    //we copy to (nn) will fit better in cache than the full original table for the
    //interpolation below.
    //The following ladder of for loops will presumably be unrolled by a reasonable compiler.
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                for (int a=0; a<=1; ++a) {
                    for (int f=0; f<=1; ++f) {
                        //Shorthands for indexing
                        const int ti = thp_i.ind_[t];
                        const int wi = wfr_i.ind_[w];
                        const int gi = gfr_i.ind_[g];
                        const int ai = alq_i.ind_[a];
                        const int fi = flo_i.ind_[f];

                        //Copy element
                        nn[t][w][g][a][f] = array[ti][wi][gi][ai][fi];
                    }
                }
            }
        }
    }

    //Remove dimensions one by one
    // Example: going from 3D to 2D to 1D, we start by interpolating along
    // the z axis first, leaving a 2D problem. Then interpolating along the y
    // axis, leaving a 1D, problem, etc.
    double tf = flo_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                for (int a=0; a<=1; ++a) {
                    nn[t][w][g][a][0] = (1.0-tf)*nn[t][w][g][a][0] + tf*nn[t][w][g][a][1];
                }
            }
        }
    }

    tf = alq_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                nn[t][w][g][0][0] = (1.0-tf)*nn[t][w][g][0][0] + tf*nn[t][w][g][1][0];
            }
        }
    }

    tf = gfr_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            nn[t][w][0][0][0] = (1.0-tf)*nn[t][w][0][0][0] + tf*nn[t][w][1][0][0];
        }
    }

    tf = wfr_i.factor_;
    for (int t=0; t<=1; ++t) {
        nn[t][0][0][0][0] = (1.0-tf)*nn[t][0][0][0][0] + tf*nn[t][1][0][0][0];
    }

    tf = thp_i.factor_;
    return (1.0-tf)*nn[0][0][0][0][0] + tf*nn[1][0][0][0][0];
}

#ifdef __GNUC__
#pragma GCC pop_options //unroll loops
#endif


double VFPProperties::find_x(const double& x0,
        const double& x1,
        const double& y0,
        const double& y1,
        const double& y) {
    const double dx = x1 - x0;
    const double dy = y1 - y0;

    /**
     *       y = y0 + (dy / dx) * (x - x0)
     *   =>  x = x0 + (y - y0) * (dx / dy)
     *
     * If dy is zero, use x1 as the value.
     */

    double x = 0.0;

    if (dy != 0.0) {
        x = x0 + (y-y0) * (dx/dy);
    }
    else {
        x = x1;
    }

    return x;
}

} //Namespace
