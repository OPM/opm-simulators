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
#include <opm/core/utility/ErrorMacros.hpp>

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


double VFPProperties::prod_bhp(int table, const double& flo, const double& thp, const double& wfr, const double& gfr, const double& alq) {
    if (m_prod_tables.find(table) == m_prod_tables.end()) {
        OPM_THROW(std::invalid_argument, "Nonexistent table " << table << " referenced.");
    }
    const auto* tab = m_prod_tables[table];
    //First, find the values to interpolate between
    auto flo_i = find_interp_data(flo, tab->getFloAxis());
    auto thp_i = find_interp_data(thp, tab->getTHPAxis());
    auto wfr_i = find_interp_data(wfr, tab->getWFRAxis());
    auto gfr_i = find_interp_data(gfr, tab->getGFRAxis());
    auto alq_i = find_interp_data(alq, tab->getALQAxis());

    //Then perform the interpolation itself
    return interpolate(tab->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);
}

VFPProperties::ADB VFPProperties::prod_bhp(int table, const ADB& flo, const ADB& thp, const ADB& wfr, const ADB& gfr, const ADB& alq) {
    const ADB::V& f_v = flo.value();
    const ADB::V& t_v = thp.value();
    const ADB::V& w_v = wfr.value();
    const ADB::V& g_v = gfr.value();
    const ADB::V& a_v = alq.value();

    const int nw = f_v.size();

    //Compute the BHP for each well independently
    ADB::V bhp_vals;
    bhp_vals.resize(nw);
    for (int i=0; i<nw; ++i) {
        bhp_vals[i] = prod_bhp(table, f_v[i], t_v[i], w_v[i], g_v[i], a_v[i]);
    }
    //Create an ADB constant value.
    return ADB::constant(bhp_vals);
}



VFPProperties::ADB VFPProperties::prod_bhp(int table, const Wells& wells, const ADB& qs, const ADB& thp, const ADB& alq) {
    if (m_prod_tables.find(table) == m_prod_tables.end()) {
        OPM_THROW(std::invalid_argument, "Nonexistant table " << table << " referenced.");
    }

    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(np == 3);
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    const auto* tab = m_prod_tables[table];
    ADB flo = getFlo(w, o, g, tab->getFloType());
    ADB wfr = getWFR(w, o, g, tab->getWFRType());
    ADB gfr = getGFR(w, o, g, tab->getGFRType());

    //TODO: Check ALQ type here?

    return prod_bhp(table, flo, thp, wfr, gfr, alq);
}


VFPProperties::ADB VFPProperties::getFlo(const ADB& aqua, const ADB& liquid, const ADB& vapour,
                                         const VFPProdTable::FLO_TYPE& type) {
    switch (type) {
        case VFPProdTable::FLO_OIL:
            //Oil = liquid phase
            return liquid;
        case VFPProdTable::FLO_LIQ:
            //Liquid = aqua + liquid phases
            return aqua + liquid;
        case VFPProdTable::FLO_GAS:
            //Gas = vapor phase
            return vapour;
        case VFPProdTable::FLO_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid FLO_TYPE: '" << type << "'");
            return ADB::null();
    }
}

VFPProperties::ADB VFPProperties::getWFR(const ADB& aqua, const ADB& liquid, const ADB& vapour,
                                         const VFPProdTable::WFR_TYPE& type) {
    switch(type) {
        case VFPProdTable::WFR_WOR:
            //Water-oil ratio = water / oil
            return aqua / liquid;
        case VFPProdTable::WFR_WCT:
            //Water cut = water / (water + oil + gas)
            return aqua / (aqua + liquid + vapour);
        case VFPProdTable::WFR_WGR:
            //Water-gas ratio = water / gas
            return aqua / vapour;
        case VFPProdTable::WFR_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid WFR_TYPE: '" << type << "'");
            return ADB::null();
    }
}


VFPProperties::ADB VFPProperties::getGFR(const ADB& aqua, const ADB& liquid, const ADB& vapour,
                                         const VFPProdTable::GFR_TYPE& type) {
    switch(type) {
        case VFPProdTable::GFR_GOR:
            // Gas-oil ratio = gas / oil
            return vapour / liquid;
        case VFPProdTable::GFR_GLR:
            // Gas-liquid ratio = gas / (oil + water)
            return vapour / (liquid + aqua);
        case VFPProdTable::GFR_OGR:
            // Oil-gas ratio = oil / gas
            return liquid / vapour;
        case VFPProdTable::GFR_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid GFR_TYPE: '" << type << "'");
            return ADB::null();
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


} //Namespace
