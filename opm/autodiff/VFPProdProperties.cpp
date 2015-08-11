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

#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/VFPProdTable.hpp>
#include <opm/core/props/BlackoilPhases.hpp>



namespace Opm {






VFPProdProperties::VFPProdProperties() {

}

VFPProdProperties::VFPProdProperties(const VFPProdTable* table){
    m_tables[table->getTableNum()] = table;
}


VFPProdProperties::VFPProdProperties(const std::map<int, VFPProdTable>& tables) {
    init(tables);
}





void VFPProdProperties::init(const std::map<int, VFPProdTable>& prod_tables) {
    //Populate production table pointers
    for (const auto& table : prod_tables) {
        m_tables[table.first] = &table.second;
    }
}

VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
        const Wells& wells,
        const ADB& qs,
        const ADB& thp,
        const ADB& alq) const {
    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(np == 3);
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    return bhp(table_id, w, o, g, thp, alq);
}

namespace detail {
    /**
     * Returns the type variable for FLO/GFR/WFR
     */
    template <typename TYPE>
    TYPE getType(const VFPProdTable* table);

    template <>
    VFPProdTable::FLO_TYPE getType(const VFPProdTable* table) {
        return table->getFloType();
    }

    template <>
    VFPProdTable::WFR_TYPE getType(const VFPProdTable* table) {
        return table->getWFRType();
    }

    template <>
    VFPProdTable::GFR_TYPE getType(const VFPProdTable* table) {
        return table->getGFRType();
    }

    /**
     * Returns the actual ADB for the type of FLO/GFR/WFR type
     */
    template <typename TYPE>
    VFPProdProperties::ADB getValue(
            const VFPProdProperties::ADB& aqua,
            const VFPProdProperties::ADB& liquid,
            const VFPProdProperties::ADB& vapour, TYPE type);

    template <>
    VFPProdProperties::ADB getValue(
            const VFPProdProperties::ADB& aqua,
            const VFPProdProperties::ADB& liquid,
            const VFPProdProperties::ADB& vapour,
            VFPProdTable::FLO_TYPE type) {
        return VFPProdProperties::getFlo(aqua, liquid, vapour, type);
    }

    template <>
    VFPProdProperties::ADB getValue(
            const VFPProdProperties::ADB& aqua,
            const VFPProdProperties::ADB& liquid,
            const VFPProdProperties::ADB& vapour,
            VFPProdTable::WFR_TYPE type) {
        return VFPProdProperties::getWFR(aqua, liquid, vapour, type);
    }

    template <>
    VFPProdProperties::ADB getValue(
            const VFPProdProperties::ADB& aqua,
            const VFPProdProperties::ADB& liquid,
            const VFPProdProperties::ADB& vapour,
            VFPProdTable::GFR_TYPE type) {
        return VFPProdProperties::getGFR(aqua, liquid, vapour, type);
    }

    /**
     * Given m wells and n types of VFP variables (e.g., FLO = {FLO_OIL, FLO_LIQ}
     * this function combines the n types of ADB objects, so that each of the
     * m wells gets the right ADB.
     */
    template <typename TYPE>
    VFPProdProperties::ADB gather_vars(const std::vector<const VFPProdTable*>& well_tables,
            const VFPProdProperties::ADB& aqua,
            const VFPProdProperties::ADB& liquid,
            const VFPProdProperties::ADB& vapour) {

        typedef VFPProdProperties::ADB ADB;

        const int num_wells = static_cast<int>(well_tables.size());
        assert(aqua.size() == num_wells);
        assert(liquid.size() == num_wells);
        assert(vapour.size() == num_wells);

        //Caching variable for flo/wfr/gfr
        std::map<TYPE, ADB> map;

        //Indexing variable used when combining the different ADB types
        std::map<TYPE, std::vector<int> > elems;

        //Compute all of the different ADB types,
        //and record which wells use which types
        for (int i=0; i<num_wells; ++i) {
            const VFPProdTable* table = well_tables[i];

            //Only do something if this well is under THP control
            if (table != NULL) {
                TYPE type = getType<TYPE>(table);

                //"Caching" of flo_type etc: Only calculate used types
                //Create type if it does not exist
                if (map.find(type) == map.end()) {
                    map.insert(std::pair<TYPE, ADB>(
                            type,
                            detail::getValue<TYPE>(aqua, liquid, vapour, type)
                            ));
                }

                //Add the index for assembly later in gather_vars
                elems[type].push_back(i);
            }
        }

        //Loop over all types of ADB variables, and combine them
        //so that each well gets the proper variable
        ADB retval = ADB::constant(ADB::V::Zero(num_wells));
        for (const auto& entry : elems) {
            const auto& key = entry.first;
            const auto& value = entry.second;

            //Get the ADB for this type of variable
            assert(map.find(key) != map.end());
            const ADB& values = map.find(key)->second;

            //Get indices to all elements that should use this ADB
            const std::vector<int>& elems = value;

            //Add these elements to retval
            retval = retval + superset(subset(values, elems), elems, values.size());
        }

        return retval;
    }

    void extendBlockPattern(const VFPProdProperties::ADB& x, std::vector<int>& block_pattern) {
        std::vector<int> x_block_pattern = x.blockPattern();

        if (x_block_pattern.empty()) {
            return;
        }
        else {
            if (block_pattern.empty()) {
                block_pattern = x_block_pattern;
                return;
            }
            else {
                if (x_block_pattern != block_pattern) {
                    OPM_THROW(std::logic_error, "Block patterns do not match");
                }
            }
        }
    }

    std::vector<int> commonBlockPattern(
            const VFPProdProperties::ADB& x1,
            const VFPProdProperties::ADB& x2,
            const VFPProdProperties::ADB& x3,
            const VFPProdProperties::ADB& x4,
            const VFPProdProperties::ADB& x5) {
        std::vector<int> block_pattern;

        extendBlockPattern(x1, block_pattern);
        extendBlockPattern(x2, block_pattern);
        extendBlockPattern(x3, block_pattern);
        extendBlockPattern(x4, block_pattern);
        extendBlockPattern(x5, block_pattern);

        return block_pattern;
    }

} //Namespace

VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        const ADB& thp,
        const ADB& alq) const {
    const int nw = thp.size();

    std::vector<int> block_pattern = detail::commonBlockPattern(aqua, liquid, vapour, thp, alq);

    assert(static_cast<int>(table_id.size()) == nw);
    assert(aqua.size()     == nw);
    assert(liquid.size()   == nw);
    assert(vapour.size()   == nw);
    assert(thp.size()      == nw);
    assert(alq.size()      == nw);

    //Allocate data for bhp's and partial derivatives
    ADB::V value = ADB::V::Zero(nw);
    ADB::V dthp = ADB::V::Zero(nw);
    ADB::V dwfr = ADB::V::Zero(nw);
    ADB::V dgfr = ADB::V::Zero(nw);
    ADB::V dalq = ADB::V::Zero(nw);
    ADB::V dflo = ADB::V::Zero(nw);

    //Get the table for each well
    std::vector<const VFPProdTable*> well_tables(nw, NULL);
    for (int i=0; i<nw; ++i) {
        if (table_id[i] >= 0) {
            well_tables[i] = getProdTable(table_id[i]);
        }
    }

    //Get the right FLO/GFR/WFR variable for each well as a single ADB
    const ADB flo = detail::gather_vars<VFPProdTable::FLO_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB wfr = detail::gather_vars<VFPProdTable::WFR_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB gfr = detail::gather_vars<VFPProdTable::GFR_TYPE>(well_tables, aqua, liquid, vapour);

    //Compute the BHP for each well independently
    for (int i=0; i<nw; ++i) {
        const VFPProdTable* table = well_tables[i];
        if (table != NULL) {
            //First, find the values to interpolate between
            auto flo_i = find_interp_data(flo.value()[i], table->getFloAxis());
            auto thp_i = find_interp_data(thp.value()[i], table->getTHPAxis());
            auto wfr_i = find_interp_data(wfr.value()[i], table->getWFRAxis());
            auto gfr_i = find_interp_data(gfr.value()[i], table->getGFRAxis());
            auto alq_i = find_interp_data(alq.value()[i], table->getALQAxis());

            adb_like bhp_val = interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

            value[i] = bhp_val.value;
            dthp[i] = bhp_val.dthp;
            dwfr[i] = bhp_val.dwfr;
            dgfr[i] = bhp_val.dgfr;
            dalq[i] = bhp_val.dalq;
            dflo[i] = bhp_val.dflo;
        }
        else {
            value[i] = -1e100; //Signal that this value has not been calculated properly, due to "missing" table
        }
    }

    //Create diagonal matrices from ADB::Vs
    ADB::M dthp_diag = spdiag(dthp);
    ADB::M dwfr_diag = spdiag(dwfr);
    ADB::M dgfr_diag = spdiag(dgfr);
    ADB::M dalq_diag = spdiag(dalq);
    ADB::M dflo_diag = spdiag(dflo);

    //Calculate the Jacobians
    const int num_blocks = block_pattern.size();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        //Could have used fastSparseProduct and temporary variables
        //but may not save too much on that.
        jacs[block] = ADB::M(nw, block_pattern[block]);

        if (!thp.derivative().empty()) {
            jacs[block] += dthp_diag * thp.derivative()[block];
        }
        if (!wfr.derivative().empty()) {
            jacs[block] += dwfr_diag * wfr.derivative()[block];
        }
        if (!gfr.derivative().empty()) {
            jacs[block] += dgfr_diag * gfr.derivative()[block];
        }
        if (!alq.derivative().empty()) {
            jacs[block] += dalq_diag * alq.derivative()[block];
        }
        if (!flo.derivative().empty()) {
            jacs[block] += dflo_diag * flo.derivative()[block];
        }
    }

    ADB retval = ADB::function(std::move(value), std::move(jacs));
    return retval;
}

VFPProdProperties::adb_like VFPProdProperties::bhp(int table_id,
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
    adb_like retval = interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);
    return retval;
}


double VFPProdProperties::thp(int table_id,
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
        bhp_array[i] = interpolate(data, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
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



const VFPProdTable* VFPProdProperties::getProdTable(int table_id) const {
    auto entry = m_tables.find(table_id);
    if (entry == m_tables.end()) {
        OPM_THROW(std::invalid_argument, "Nonexistent table " << table_id << " referenced.");
    }
    else {
        return entry->second;
    }
}

VFPProdProperties::InterpData VFPProdProperties::find_interp_data(const double& value, const std::vector<double>& values) {
    InterpData retval;

    //If we only have one value in our vector, return that
    if (values.size() == 1) {
        retval.ind_[0] = 0;
        retval.ind_[1] = 0;
        retval.inv_dist_ = 0.0;
        retval.factor_ = 0.0;
    }
    // Else search in the vector
    else {
        //First element greater than or equal to value
        //Start with the second element, so that floor_iter does not go out of range
        //Don't access out-of-range, therefore values.end()-1
        auto ceil_iter = std::lower_bound(values.begin()+1, values.end()-1, value);

        //Find last element smaller than range
        auto floor_iter = ceil_iter-1;

        //Find the indices
        retval.ind_[0] = floor_iter - values.begin();
        retval.ind_[1] = ceil_iter - values.begin();

        //Find interpolation ratio
        double dist = (*ceil_iter - *floor_iter);
        if (std::abs(dist) > 0.0) {
            //Possible source for floating point error here if value and floor are large,
            //but very close to each other
            retval.inv_dist_ = 1.0 / dist;
            retval.factor_ = (value-*floor_iter) * retval.inv_dist_;
        }
        else {
            retval.inv_dist_ = 0.0;
            retval.factor_ = 0.0;
        }
    }

    return retval;
}


#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif

VFPProdProperties::adb_like VFPProdProperties::interpolate(
        const VFPProdTable::array_type& array,
        const InterpData& flo_i,
        const InterpData& thp_i,
        const InterpData& wfr_i,
        const InterpData& gfr_i,
        const InterpData& alq_i) {

    //Values and derivatives in a 5D hypercube
    adb_like nn[2][2][2][2][2];


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
                        nn[t][w][g][a][f].value = array[ti][wi][gi][ai][fi];
                    }
                }
            }
        }
    }

    //Calculate derivatives
    //Note that the derivative of the two end points of a line aligned with the
    //"axis of the derivative" are equal
    for (int i=0; i<=1; ++i) {
        for (int j=0; j<=1; ++j) {
            for (int k=0; k<=1; ++k) {
                for (int l=0; l<=1; ++l) {
                    nn[0][i][j][k][l].dthp = (nn[1][i][j][k][l].value - nn[0][i][j][k][l].value) * thp_i.inv_dist_;
                    nn[i][0][j][k][l].dwfr = (nn[i][1][j][k][l].value - nn[i][0][j][k][l].value) * wfr_i.inv_dist_;
                    nn[i][j][0][k][l].dgfr = (nn[i][j][1][k][l].value - nn[i][j][0][k][l].value) * gfr_i.inv_dist_;
                    nn[i][j][k][0][l].dalq = (nn[i][j][k][1][l].value - nn[i][j][k][0][l].value) * alq_i.inv_dist_;
                    nn[i][j][k][l][0].dflo = (nn[i][j][k][l][1].value - nn[i][j][k][l][0].value) * flo_i.inv_dist_;

                    nn[1][i][j][k][l].dthp = nn[0][i][j][k][l].dthp;
                    nn[i][1][j][k][l].dwfr = nn[i][0][j][k][l].dwfr;
                    nn[i][j][1][k][l].dgfr = nn[i][j][0][k][l].dgfr;
                    nn[i][j][k][1][l].dalq = nn[i][j][k][0][l].dalq;
                    nn[i][j][k][l][1].dflo = nn[i][j][k][l][0].dflo;
                }
            }
        }
    }

    double t1, t2; //interpolation variables, so that t1 = (1-t) and t2 = t.

    // Remove dimensions one by one
    // Example: going from 3D to 2D to 1D, we start by interpolating along
    // the z axis first, leaving a 2D problem. Then interpolating along the y
    // axis, leaving a 1D, problem, etc.
    t2 = flo_i.factor_;
    t1 = (1.0-t2);
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                for (int a=0; a<=1; ++a) {
                    nn[t][w][g][a][0] = t1*nn[t][w][g][a][0] + t2*nn[t][w][g][a][1];
                }
            }
        }
    }

    t2 = alq_i.factor_;
    t1 = (1.0-t2);
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                nn[t][w][g][0][0] = t1*nn[t][w][g][0][0] + t2*nn[t][w][g][1][0];
            }
        }
    }

    t2 = gfr_i.factor_;
    t1 = (1.0-t2);
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            nn[t][w][0][0][0] = t1*nn[t][w][0][0][0] + t2*nn[t][w][1][0][0];
        }
    }

    t2 = wfr_i.factor_;
    t1 = (1.0-t2);
    for (int t=0; t<=1; ++t) {
        nn[t][0][0][0][0] = t1*nn[t][0][0][0][0] + t2*nn[t][1][0][0][0];
    }

    t2 = thp_i.factor_;
    t1 = (1.0-t2);
    nn[0][0][0][0][0] = t1*nn[0][0][0][0][0] + t2*nn[1][0][0][0][0];

    return nn[0][0][0][0][0];
}

#ifdef __GNUC__
#pragma GCC pop_options //unroll loops
#endif


double VFPProdProperties::find_x(const double& x0,
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



}
