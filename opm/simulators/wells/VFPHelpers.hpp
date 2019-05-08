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

#include <opm/common/OpmLog/OpmLog.hpp>

#include <cmath>
#include <opm/common/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {
namespace detail {


/**
 * Returns zero if input value is NaN of INF
 */
inline double zeroIfNanInf(const double& value) {
    const bool nan_or_inf = std::isnan(value) || std::isinf(value);

    if (nan_or_inf) {
        OpmLog::warning("NAN_OR_INF_VFP", "NAN or INF value encountered during VFP calculation, the value is set to zero");
    }

    return nan_or_inf ? 0.0 : value;
}


/**
 * Returns zero if input value is NaN or INF
 */
template <class EvalWell>
inline EvalWell zeroIfNanInf(const EvalWell& value) {
    const bool nan_or_inf = std::isnan(value.value()) || std::isinf(value.value());

    if (nan_or_inf) {
        OpmLog::warning("NAN_OR_INF_VFP_EVAL", "NAN or INF Evalution encountered during VFP calculation, the Evalution is set to zero");
    }

    using Toolbox = MathToolbox<EvalWell>;

    return nan_or_inf ? Toolbox::createBlank(value) : value;
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
            return zeroIfNanInf(wor);
        }
        case VFPProdTable::WFR_WCT:
            //Water cut = water / (water + oil)
            return zeroIfNanInf(aqua / (aqua + liquid));
        case VFPProdTable::WFR_WGR:
            //Water-gas ratio = water / gas
            return zeroIfNanInf(aqua / vapour);
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
            return zeroIfNanInf(vapour / liquid);
        case VFPProdTable::GFR_GLR:
            // Gas-liquid ratio = gas / (oil + water)
            return zeroIfNanInf(vapour / (liquid + aqua));
        case VFPProdTable::GFR_OGR:
            // Oil-gas ratio = oil / gas
            return zeroIfNanInf(liquid / vapour);
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
        OPM_THROW(std::invalid_argument, "Nonexistent VFP table " << table_id << " referenced.");
    }
    else {
        return entry->second;
    }
}

/**
 * Check whether we have a table with the table number
 */
template <typename T>
bool hasTable(const std::map<int, T*> tables, int table_id) {
    const auto entry = tables.find(table_id);
    return (entry != tables.end() );
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



// a data type use to do the intersection calculation to get the intial bhp under THP control
struct RateBhpPair {
    double rate;
    double bhp;
};


// looking for a intersection point a line segment and a line, they are both defined with two points
// it is copied from #include <opm/polymer/Point2D.hpp>, which should be removed since it is only required by the lagacy polymer
inline bool findIntersection(const std::array<RateBhpPair, 2>& line_segment, const std::array<RateBhpPair, 2>& line, double& bhp) {
    const double x1 = line_segment[0].rate;
    const double y1 = line_segment[0].bhp;
    const double x2 = line_segment[1].rate;
    const double y2 = line_segment[1].bhp;

    const double x3 = line[0].rate;
    const double y3 = line[0].bhp;
    const double x4 = line[1].rate;
    const double y4 = line[1].bhp;

    const double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    if (d == 0.) {
        return false;
    }

    const double x = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
    const double y = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

    if (x >= std::min(x1,x2) && x <= std::max(x1,x2)) {
        bhp = y;
        return true;
    } else {
        return false;
    }
}

// calculating the BHP from thp through the intersection of VFP curves and inflow performance relationship
inline bool findIntersectionForBhp(const std::vector<RateBhpPair>& ratebhp_samples,
                                   const std::array<RateBhpPair, 2>& ratebhp_twopoints_ipr,
                                   double& obtained_bhp)
{
    // there possibly two intersection point, then we choose the one corresponding with the bigger rate

    const double bhp1 = ratebhp_twopoints_ipr[0].bhp;
    const double rate1 = ratebhp_twopoints_ipr[0].rate;

    const double bhp2 = ratebhp_twopoints_ipr[1].bhp;
    const double rate2 = ratebhp_twopoints_ipr[1].rate;

    assert(rate1 != rate2);

    const double line_slope = (bhp2 - bhp1) / (rate2 - rate1);

    // line equation will be
    // bhp - bhp1 - line_slope * (flo_rate - flo_rate1) = 0
    auto flambda = [&](const double flo_rate, const double bhp) {
        return bhp - bhp1 - line_slope * (flo_rate - rate1);
    };

    int number_intersection_found = 0;
    int index_segment = 0; // the intersection segment that intersection happens
    const size_t num_samples = ratebhp_samples.size();
    for (size_t i = 0; i < num_samples - 1; ++i) {
        const double temp1 = flambda(ratebhp_samples[i].rate, ratebhp_samples[i].bhp);
        const double temp2 = flambda(ratebhp_samples[i+1].rate, ratebhp_samples[i+1].bhp);
        if (temp1 * temp2 <= 0.) { // intersection happens
            // in theory there should be maximum two intersection points
            // while considering the situation == 0. here, we might find more
            // we always use the last one, which is the one corresponds to the biggest rate,
            // which we assume is the more stable one
            ++number_intersection_found;
            index_segment = i;
        }
    }

    if (number_intersection_found == 0) { // there is not intersection point
        return false;
    }

    // then we pick the segment from the VFP curve to do the line intersection calculation
    const std::array<RateBhpPair, 2> line_segment{ ratebhp_samples[index_segment], ratebhp_samples[index_segment + 1] };

    const bool intersection_found = findIntersection(line_segment, ratebhp_twopoints_ipr, obtained_bhp);

    return intersection_found;
}


} // namespace detail


} // namespace




#endif /* OPM_AUTODIFF_VFPHELPERS_HPP_ */
