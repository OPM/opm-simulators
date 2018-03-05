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


#ifndef OPM_AUTODIFF_VFPHELPERS_HPP_
#define OPM_AUTODIFF_VFPHELPERS_HPP_


#include <opm/parser/eclipse/EclipseState/Tables/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/VFPInjTable.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {
namespace detail {


typedef AutoDiffBlock<double> ADB;


/**
 * Returns zero if input value is NaN
 */
inline double zeroIfNan(const double& value) {
    return (std::isnan(value)) ? 0.0 : value;
}

/**
 * Returns zero if input value is NaN
 */
template <class EvalWell>
inline double zeroIfNan(const EvalWell& value) {
    return (std::isnan(value.value())) ? 0.0 : value.value();
}



/**
 * Returns zero for every entry in the ADB which is NaN
 */
inline ADB zeroIfNan(const ADB& values) {
    Selector<ADB::V::Scalar> not_nan_selector(values.value(), Selector<ADB::V::Scalar>::NotNaN);

    const ADB::V z = ADB::V::Zero(values.size());
    const ADB zero = ADB::constant(z, values.blockPattern());

    ADB retval = not_nan_selector.select(values, zero);
    return retval;
}





/**
 * Computes the flo parameter according to the flo_type_
 * for production tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
static T getFlo(const T& aqua, const T& liquid, const T& vapour,
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
    }
}



/**
 * Computes the flo parameter according to the flo_type_
 * for injection tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
static T getFlo(const T& aqua, const T& liquid, const T& vapour,
                  const VFPInjTable::FLO_TYPE& type) {
    switch (type) {
        case VFPInjTable::FLO_OIL:
            //Oil = liquid phase
            return liquid;
        case VFPInjTable::FLO_WAT:
            //Liquid = aqua phase
            return aqua;
        case VFPInjTable::FLO_GAS:
            //Gas = vapor phase
            return vapour;
        case VFPInjTable::FLO_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid FLO_TYPE: '" << type << "'");
    }
}







/**
 * Computes the wfr parameter according to the wfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
static T getWFR(const T& aqua, const T& liquid, const T& vapour,
                  const VFPProdTable::WFR_TYPE& type) {
    switch(type) {
        case VFPProdTable::WFR_WOR: {
            //Water-oil ratio = water / oil
            T wor = aqua / liquid;
            return zeroIfNan(wor);
        }
        case VFPProdTable::WFR_WCT:
            //Water cut = water / (water + oil)
            return zeroIfNan(aqua / (aqua + liquid));
        case VFPProdTable::WFR_WGR:
            //Water-gas ratio = water / gas
            return zeroIfNan(aqua / vapour);
        case VFPProdTable::WFR_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid WFR_TYPE: '" << type << "'");
    }
}






/**
 * Computes the gfr parameter according to the gfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
static T getGFR(const T& aqua, const T& liquid, const T& vapour,
                  const VFPProdTable::GFR_TYPE& type) {
    switch(type) {
        case VFPProdTable::GFR_GOR:
            // Gas-oil ratio = gas / oil
            return zeroIfNan(vapour / liquid);
        case VFPProdTable::GFR_GLR:
            // Gas-liquid ratio = gas / (oil + water)
            return zeroIfNan(vapour / (liquid + aqua));
        case VFPProdTable::GFR_OGR:
            // Oil-gas ratio = oil / gas
            return zeroIfNan(liquid / vapour);
        case VFPProdTable::GFR_INVALID: //Intentional fall-through
        default:
            OPM_THROW(std::logic_error, "Invalid GFR_TYPE: '" << type << "'");
    }
}






/**
 * Helper struct for linear interpolation
 */
struct InterpData {
    InterpData() : ind_{0, 0}, inv_dist_(0.0), factor_(0.0) {}
    int ind_[2]; //[First element greater than or equal to value, Last element smaller than or equal to value]
    double inv_dist_; // 1 / distance between the two end points of the segment. Used to calculate derivatives and uses 1.0 / 0.0 = 0.0 as a convention
    double factor_; // Interpolation factor
};






/**
 * Helper function to find indices etc. for linear interpolation and extrapolation
 *  @param value Value to find in values
 *  @param values Sorted list of values to search for value in.
 *  @return Data required to find the interpolated value
 */
inline InterpData findInterpData(const double& value, const std::vector<double>& values) {
    InterpData retval;

    const int nvalues = values.size();

    //If we only have one value in our vector, return that
    if (values.size() == 1) {
        retval.ind_[0] = 0;
        retval.ind_[1] = 0;
        retval.inv_dist_ = 0.0;
        retval.factor_ = 0.0;
    }
    // Else search in the vector
    else {
        //If value is less than all values, use first interval
        if (value < values.front()) {
            retval.ind_[0] = 0;
            retval.ind_[1] = 1;
        }
        //If value is greater than all values, use last interval
        else if (value >= values.back()) {
            retval.ind_[0] = nvalues-2;
            retval.ind_[1] = nvalues-1;
        }
        else {
            //Search internal intervals
            for (int i=1; i<nvalues; ++i) {
                if (values[i] >= value) {
                    retval.ind_[0] = i-1;
                    retval.ind_[1] = i;
                    break;
                }
            }
        }

        const double start = values[retval.ind_[0]];
        const double end   = values[retval.ind_[1]];

        //Find interpolation ratio
        if (end > start) {
            //FIXME: Possible source for floating point error here if value and floor are large,
            //but very close to each other
            retval.inv_dist_ = 1.0 / (end-start);
            retval.factor_ = (value-start) * retval.inv_dist_;
        }
        else {
            retval.inv_dist_ = 0.0;
            retval.factor_ = 0.0;
        }
    }

    return retval;
}






/**
 * An "ADB-like" structure with a single value and a set of derivatives
 */
struct VFPEvaluation {
    VFPEvaluation() : value(0.0), dthp(0.0), dwfr(0.0), dgfr(0.0), dalq(0.0), dflo(0.0) {};
    double value;
    double dthp;
    double dwfr;
    double dgfr;
    double dalq;
    double dflo;
};

inline VFPEvaluation operator+(
        VFPEvaluation lhs,
        const VFPEvaluation& rhs) {
    lhs.value += rhs.value;
    lhs.dthp += rhs.dthp;
    lhs.dwfr += rhs.dwfr;
    lhs.dgfr += rhs.dgfr;
    lhs.dalq += rhs.dalq;
    lhs.dflo += rhs.dflo;
    return lhs;
}

inline VFPEvaluation operator-(
        VFPEvaluation lhs,
        const VFPEvaluation& rhs) {
    lhs.value -= rhs.value;
    lhs.dthp -= rhs.dthp;
    lhs.dwfr -= rhs.dwfr;
    lhs.dgfr -= rhs.dgfr;
    lhs.dalq -= rhs.dalq;
    lhs.dflo -= rhs.dflo;
    return lhs;
}

inline VFPEvaluation operator*(
        double lhs,
        const VFPEvaluation& rhs) {
    VFPEvaluation retval;
    retval.value = rhs.value * lhs;
    retval.dthp = rhs.dthp * lhs;
    retval.dwfr = rhs.dwfr * lhs;
    retval.dgfr = rhs.dgfr * lhs;
    retval.dalq = rhs.dalq * lhs;
    retval.dflo = rhs.dflo * lhs;
    return retval;
}






/**
 * Helper function which interpolates data using the indices etc. given in the inputs.
 */
inline VFPEvaluation interpolate(
        const VFPProdTable::array_type& array,
        const InterpData& flo_i,
        const InterpData& thp_i,
        const InterpData& wfr_i,
        const InterpData& gfr_i,
        const InterpData& alq_i) {

    //Values and derivatives in a 5D hypercube
    VFPEvaluation nn[2][2][2][2][2];


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





/**
 * This basically models interpolate(VFPProdTable::array_type, ...)
 * which performs 5D interpolation, but here for the 2D case only
 */
inline VFPEvaluation interpolate(
        const VFPInjTable::array_type& array,
        const InterpData& flo_i,
        const InterpData& thp_i) {

    //Values and derivatives in a 2D plane
    VFPEvaluation nn[2][2];


    //Pick out nearest neighbors (nn) to our evaluation point
    //The following ladder of for loops will presumably be unrolled by a reasonable compiler.
    for (int t=0; t<=1; ++t) {
        for (int f=0; f<=1; ++f) {
            //Shorthands for indexing
            const int ti = thp_i.ind_[t];
            const int fi = flo_i.ind_[f];

            //Copy element
            nn[t][f].value = array[ti][fi];
        }
    }

    //Calculate derivatives
    //Note that the derivative of the two end points of a line aligned with the
    //"axis of the derivative" are equal
    for (int i=0; i<=1; ++i) {
        nn[0][i].dthp = (nn[1][i].value - nn[0][i].value) * thp_i.inv_dist_;
        nn[i][0].dwfr = -1e100;
        nn[i][0].dgfr = -1e100;
        nn[i][0].dalq = -1e100;
        nn[i][0].dflo = (nn[i][1].value - nn[i][0].value) * flo_i.inv_dist_;

        nn[1][i].dthp = nn[0][i].dthp;
        nn[i][1].dwfr = nn[i][0].dwfr;
        nn[i][1].dgfr = nn[i][0].dgfr;
        nn[i][1].dalq = nn[i][0].dalq;
        nn[i][1].dflo = nn[i][0].dflo;
    }

    double t1, t2; //interpolation variables, so that t1 = (1-t) and t2 = t.

    // Go from 2D to 1D
    t2 = flo_i.factor_;
    t1 = (1.0-t2);
    nn[0][0] = t1*nn[0][0] + t2*nn[0][1];
    nn[1][0] = t1*nn[1][0] + t2*nn[1][1];

    // Go from line to point on line
    t2 = thp_i.factor_;
    t1 = (1.0-t2);
    nn[0][0] = t1*nn[0][0] + t2*nn[1][0];

    return nn[0][0];
}

inline VFPEvaluation bhp(const VFPProdTable* table,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp,
        const double& alq) {
    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = detail::getWFR(aqua, liquid, vapour, table->getWFRType());
    double gfr = detail::getGFR(aqua, liquid, vapour, table->getGFRType());

    //First, find the values to interpolate between
    //Recall that flo is negative in Opm, so switch sign.
    auto flo_i = detail::findInterpData(-flo, table->getFloAxis());
    auto thp_i = detail::findInterpData( thp, table->getTHPAxis());
    auto wfr_i = detail::findInterpData( wfr, table->getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table->getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table->getALQAxis());

    detail::VFPEvaluation retval = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

    return retval;
}





inline VFPEvaluation bhp(const VFPInjTable* table,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp) {
    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());

    //First, find the values to interpolate between
    auto flo_i = detail::findInterpData(flo, table->getFloAxis());
    auto thp_i = detail::findInterpData(thp, table->getTHPAxis());

    //Then perform the interpolation itself
    detail::VFPEvaluation retval = detail::interpolate(table->getTable(), flo_i, thp_i);

    return retval;
}








/**
 * Returns the table from the map if found, or throws an exception
 */
template <typename T>
const T* getTable(const std::map<int, T*> tables, int table_id) {
    auto entry = tables.find(table_id);
    if (entry == tables.end()) {
        OPM_THROW(std::invalid_argument, "Nonexistent table " << table_id << " referenced.");
    }
    else {
        return entry->second;
    }
}










/**
 * Sets block_pattern to be the "union of x.blockPattern() and block_pattern".
 */
inline void extendBlockPattern(const ADB& x, std::vector<int>& block_pattern) {
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

/**
 * Finds the common block pattern for all inputs
 */
inline std::vector<int> commonBlockPattern(
        const ADB& x1,
        const ADB& x2,
        const ADB& x3,
        const ADB& x4) {
    std::vector<int> block_pattern;

    extendBlockPattern(x1, block_pattern);
    extendBlockPattern(x2, block_pattern);
    extendBlockPattern(x3, block_pattern);
    extendBlockPattern(x4, block_pattern);

    return block_pattern;
}

inline std::vector<int> commonBlockPattern(
        const ADB& x1,
        const ADB& x2,
        const ADB& x3,
        const ADB& x4,
        const ADB& x5) {
    std::vector<int> block_pattern = commonBlockPattern(x1, x2, x3, x4);
    extendBlockPattern(x5, block_pattern);

    return block_pattern;
}












/**
 * Returns the type variable for FLO/GFR/WFR for production tables
 */
template <typename TYPE, typename TABLE>
TYPE getType(const TABLE* table);

template <>
inline
VFPProdTable::FLO_TYPE getType(const VFPProdTable* table) {
    return table->getFloType();
}

template <>
inline
VFPProdTable::WFR_TYPE getType(const VFPProdTable* table) {
    return table->getWFRType();
}

template <>
inline
VFPProdTable::GFR_TYPE getType(const VFPProdTable* table) {
    return table->getGFRType();
}


/**
 * Returns the type variable for FLO for injection tables
 */
template <>
inline
VFPInjTable::FLO_TYPE getType(const VFPInjTable* table) {
    return table->getFloType();
}




/**
 * Returns the actual ADB for the type of FLO/GFR/WFR type
 */
template <typename TYPE>
ADB getValue(
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour, TYPE type);

template <>
inline
ADB getValue(
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        VFPProdTable::FLO_TYPE type) {
    return detail::getFlo(aqua, liquid, vapour, type);
}

template <>
inline
ADB getValue(
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        VFPProdTable::WFR_TYPE type) {
    return detail::getWFR(aqua, liquid, vapour, type);
}

template <>
inline
ADB getValue(
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        VFPProdTable::GFR_TYPE type) {
    return detail::getGFR(aqua, liquid, vapour, type);
}

template <>
inline
ADB getValue(
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        VFPInjTable::FLO_TYPE type) {
    return detail::getFlo(aqua, liquid, vapour, type);
}

/**
 * Given m wells and n types of VFP variables (e.g., FLO = {FLO_OIL, FLO_LIQ}
 * this function combines the n types of ADB objects, so that each of the
 * m wells gets the right ADB.
 * @param TYPE Type of variable to return, e.g., FLO_TYPE, WFR_TYPE, GFR_TYPE
 * @param TABLE Type of table to use, e.g., VFPInjTable, VFPProdTable.
 */
template <typename TYPE, typename TABLE>
ADB combineADBVars(const std::vector<const TABLE*>& well_tables,
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour) {

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
        const TABLE* table = well_tables[i];

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
        const std::vector<int>& current = value;

        //Add these elements to retval
        retval = retval + superset(subset(values, current), current, values.size());
    }

    return retval;
}

/**
 * Helper function that finds x for a given value of y for a line
 * *NOTE ORDER OF ARGUMENTS*
 */
inline double findX(const double& x0,
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










/**
 * This function finds the value of THP given a specific BHP.
 * Essentially:
 *   Given the function f(thp_array(x)) = bhp_array(x), which is piecewise linear,
 *   find thp so that f(thp) = bhp.
 */
inline double findTHP(
        const std::vector<double>& bhp_array,
        const std::vector<double>& thp_array,
        double bhp) {
    int nthp = thp_array.size();

    double thp = -1e100;

    //Check that our thp axis is sorted
    assert(std::is_sorted(thp_array.begin(), thp_array.end()));

    /**
     * Our *interpolated* bhp_array will be montonic increasing for increasing
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
            thp = detail::findX(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
            thp = detail::findX(x0, x1, y0, y1, bhp);
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
            static_cast<void>(found); //Silence compiler warning

            const double& x0 = thp_array[i  ];
            const double& x1 = thp_array[i+1];
            const double& y0 = bhp_array[i  ];
            const double& y1 = bhp_array[i+1];
            thp = detail::findX(x0, x1, y0, y1, bhp);
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
            thp = detail::findX(x0, x1, y0, y1, bhp);
        }
        else if (bhp <= bhp_array[0]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[0];
            const double& x1 = thp_array[1];
            const double& y0 = bhp_array[0];
            const double& y1 = bhp_array[1];
            thp = detail::findX(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
            thp = detail::findX(x0, x1, y0, y1, bhp);
        }
        else {
            OPM_THROW(std::logic_error, "Programmer error: Unable to find THP in THP array");
        }
    }

    return thp;
}







} // namespace detail


} // namespace




#endif /* OPM_AUTODIFF_VFPHELPERS_HPP_ */
