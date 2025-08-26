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

#include <config.h>
#include <opm/simulators/wells/VFPHelpers.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

/**
 * Helper function that finds x for a given value of y for a line
 * *NOTE ORDER OF ARGUMENTS*
 */
template<class Scalar>
Scalar findX(const Scalar x0,
             const Scalar x1,
             const Scalar y0,
             const Scalar y1,
             const Scalar y)
{
    const Scalar dx = x1 - x0;
    const Scalar dy = y1 - y0;

    /**
     *       y = y0 + (dy / dx) * (x - x0)
     *   =>  x = x0 + (y - y0) * (dx / dy)
     *
     * If dy is zero, use x1 as the value.
     */

    Scalar x = 0.0;

    if (dy != 0.0) {
        x = x0 + (y-y0) * (dx/dy);
    }
    else {
        x = x1;
    }

    return x;
}

/**
 * Returns zero if input value is negative
 */
template <typename T>
static T chopNegativeValues(const T& value) {
    return Opm::max(0.0, value);
}

}

namespace Opm {

template<class Scalar>
detail::InterpData<Scalar> VFPHelpers<Scalar>::findInterpData(const Scalar value_in,
                                                              const std::vector<double>& values)
{
    detail::InterpData<Scalar> retval;

    const int nvalues = values.size();

    // chopping the value to be zero, which means we do not
    // extrapolate the table towards nagative ranges
    const Scalar value = value_in < 0.? 0. : value_in;

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

        const Scalar start = values[retval.ind_[0]];
        const Scalar end   = values[retval.ind_[1]];

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

template<class Scalar>
detail::VFPEvaluation<Scalar> VFPHelpers<Scalar>::
interpolate(const VFPProdTable& table,
            const detail::InterpData<Scalar>& flo_i,
            const detail::InterpData<Scalar>& thp_i,
            const detail::InterpData<Scalar>& wfr_i,
            const detail::InterpData<Scalar>& gfr_i,
            const detail::InterpData<Scalar>& alq_i)
{
    //Values and derivatives in a 5D hypercube
    detail::VFPEvaluation<Scalar> nn[2][2][2][2][2];

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
                        nn[t][w][g][a][f].value = table(ti,wi,gi,ai,fi);
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

    Scalar t1, t2; //interpolation variables, so that t1 = (1-t) and t2 = t.

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

template<class Scalar>
detail::VFPEvaluation<Scalar> VFPHelpers<Scalar>::
interpolate(const VFPInjTable& table,
            const detail::InterpData<Scalar>& flo_i,
            const detail::InterpData<Scalar>& thp_i)
{
    //Values and derivatives in a 2D plane
    detail::VFPEvaluation<Scalar> nn[2][2];

    //Pick out nearest neighbors (nn) to our evaluation point
    //The following ladder of for loops will presumably be unrolled by a reasonable compiler.
    for (int t=0; t<=1; ++t) {
        for (int f=0; f<=1; ++f) {
            //Shorthands for indexing
            const int ti = thp_i.ind_[t];
            const int fi = flo_i.ind_[f];

            //Copy element
            nn[t][f].value = table(ti,fi);
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

    Scalar t1, t2; //interpolation variables, so that t1 = (1-t) and t2 = t.

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

template<class Scalar>
detail::VFPEvaluation<Scalar> VFPHelpers<Scalar>::
bhp(const VFPProdTable& table,
    const Scalar aqua,
    const Scalar liquid,
    const Scalar vapour,
    const Scalar thp,
    const Scalar alq,
    const Scalar explicit_wfr,
    const Scalar explicit_gfr,
    const bool   use_vfpexplicit)
{
    //Find interpolation variables
    Scalar flo = detail::getFlo(table, aqua, liquid, vapour);
    Scalar wfr = detail::getWFR(table, aqua, liquid, vapour);
    Scalar gfr = detail::getGFR(table, aqua, liquid, vapour);
    if (use_vfpexplicit || -flo < table.getFloAxis().front()) {
        wfr = explicit_wfr;
        gfr = explicit_gfr;
    }

    //First, find the values to interpolate between
    //Recall that flo is negative in Opm, so switch sign.
    auto flo_i = findInterpData(-flo, table.getFloAxis());
    auto thp_i = findInterpData( thp, table.getTHPAxis());
    auto wfr_i = findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = findInterpData( gfr, table.getGFRAxis());
    auto alq_i = findInterpData( alq, table.getALQAxis());

    detail::VFPEvaluation retval = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);

    return retval;
}

template<class Scalar>
detail::VFPEvaluation<Scalar> VFPHelpers<Scalar>::
bhp(const VFPInjTable& table,
    const Scalar aqua,
    const Scalar liquid,
    const Scalar vapour,
    const Scalar thp)
{
    //Find interpolation variables
    Scalar flo = detail::getFlo(table, aqua, liquid, vapour);

    //First, find the values to interpolate between
    auto flo_i = findInterpData(flo, table.getFloAxis());
    auto thp_i = findInterpData(thp, table.getTHPAxis());

    //Then perform the interpolation itself
    detail::VFPEvaluation retval = interpolate(table, flo_i, thp_i);

    return retval;
}

template<class Scalar>
Scalar VFPHelpers<Scalar>::
findTHP(const std::vector<Scalar>& bhp_array,
        const std::vector<double>& thp_array,
        Scalar bhp, 
        const bool find_largest)
{
    int nthp = thp_array.size();

    if (!std::isfinite(bhp)) {
        throw NumericalProblem("findTHP: Error bhp is not finite");
    }
    Scalar thp = -1e100;

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
            const Scalar& x0 = thp_array[0];
            const Scalar& x1 = thp_array[1];
            const Scalar& y0 = bhp_array[0];
            const Scalar& y1 = bhp_array[1];
            thp = findX(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const Scalar& x0 = thp_array[nthp-2];
            const Scalar& x1 = thp_array[nthp-1];
            const Scalar& y0 = bhp_array[nthp-2];
            const Scalar& y1 = bhp_array[nthp-1];
            thp = findX(x0, x1, y0, y1, bhp);
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
                const Scalar& y0 = bhp_array[i  ];
                const Scalar& y1 = bhp_array[i+1];

                if (y0 < bhp && bhp <= y1) {
                    found = true;
                    break;
                }
            }
            //Canary in a coal mine: shouldn't really be required
            assert(found == true);
            static_cast<void>(found); //Silence compiler warning

            const Scalar& x0 = thp_array[i  ];
            const Scalar& x1 = thp_array[i+1];
            const Scalar& y0 = bhp_array[i  ];
            const Scalar& y1 = bhp_array[i+1];
            thp = findX(x0, x1, y0, y1, bhp);
        }
    }
    //bhp_array not sorted, raw search.
    else {
        //Here we're into damage prevention territory, and there may be any number of  
        //solutions (including zero). The well is currently not controlled by THP, and
        //since we're doing severe extrapolaton we would also like, if possible, to prevent 
        //it from switcing to THP. Accordingly, if there are multiple solutions, we return 
        //the value for the intersection corresponding to the largest (smallest) THP-value 
        //for producers (injectors). 

        //first check which extrapolations are valid
        const bool first_slope_positive = bhp_array[1] >= bhp_array[0];
        const bool valid_low = (bhp < bhp_array[0] && first_slope_positive) || (bhp >= bhp_array[0] && !first_slope_positive);
        const bool last_slope_positive = bhp_array[nthp-1] >= bhp_array[nthp-2];
        const bool valid_high = (bhp > bhp_array[nthp-1] && last_slope_positive) || (bhp <= bhp_array[nthp-1] && !last_slope_positive);

        bool found = false;
        int array_ix = 0;
        if (find_largest){//find intersection corresponding to the largest thp
            // high extrap -> table interp -> low extrap
            if (valid_high) {
                found = true;
                array_ix = nthp-2;
            } else {
                //search backward within table
                for (int i = nthp-2; i>=0; --i) {
                    const Scalar& y0 = bhp_array[i  ];
                    const Scalar& y1 = bhp_array[i+1];
                    if (std::min(y0, y1) < bhp && bhp <= std::max(y0, y1)) {
                        found = true;
                        array_ix = i;
                        break;
                    }
                }
                if (!found && valid_low) {
                    found = true;
                    array_ix = 0;
                }
            }
        } else {//find intersection corresponding to the smallest thp
            //low extrap -> table interp -> high extrap
            if (valid_low) {
                found = true;
                array_ix = 0;
            } else {
                //search forward within table
                for (int i = 0; i<nthp-1; ++i) {
                    const Scalar& y0 = bhp_array[i  ];
                    const Scalar& y1 = bhp_array[i+1];
                    if (std::min(y0, y1) < bhp && bhp <= std::max(y0, y1)) {
                        found = true;
                        array_ix = i;
                        break;
                    }
                }
                if (!found && valid_high) {
                    found = true;
                    array_ix = nthp-2;
                }
            }
        }
        if (found) {
            const Scalar& x0 = thp_array[array_ix  ];
            const Scalar& x1 = thp_array[array_ix+1];
            const Scalar& y0 = bhp_array[array_ix  ];
            const Scalar& y1 = bhp_array[array_ix+1];
            thp = findX(x0, x1, y0, y1, bhp);
        } else {
            // no intersection, just return largest/smallest value in table
            if (find_largest) {
                thp = thp_array[nthp-1];
            } else {
                thp = thp_array[0];
            }
        }
    }
    return thp;
}

template<class Scalar>
std::pair<Scalar, Scalar> VFPHelpers<Scalar>::
getMinimumBHPCoordinate(const VFPProdTable& table,
                        const Scalar thp,
                        const Scalar wfr,
                        const Scalar gfr,
                        const Scalar alq)
{
    // Given fixed thp, wfr, gfr and alq, this function finds the minimum bhp and returns
    // the corresponding pair (-flo_at_bhp_min, bhp_min). No assumption is taken on the
    // shape of the function bhp(flo), so all points in the flo-axis is checked.
    Scalar flo_at_bhp_min = 0.0; // start by checking flo=0
    auto flo_i = findInterpData(flo_at_bhp_min, table.getFloAxis());
    auto thp_i = findInterpData( thp, table.getTHPAxis());
    auto wfr_i = findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = findInterpData( gfr, table.getGFRAxis());
    auto alq_i = findInterpData( alq, table.getALQAxis());

    detail::VFPEvaluation bhp_i = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
    Scalar bhp_min = bhp_i.value;
    const std::vector<double>& flos = table.getFloAxis();
    for (size_t i = 0; i < flos.size(); ++i) {
        flo_i = findInterpData(flos[i], flos);
        bhp_i = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        if (bhp_i.value < bhp_min){
            bhp_min = bhp_i.value;
            flo_at_bhp_min = flos[i];
        }
    }
    // return negative flo
    return std::make_pair(-flo_at_bhp_min, bhp_min);
}

template<class Scalar>
std::optional<std::pair<Scalar, Scalar>> VFPHelpers<Scalar>::
intersectWithIPR(const VFPProdTable& table,
                 const Scalar thp,
                 const Scalar wfr,
                 const Scalar gfr,
                 const Scalar alq,
                 const Scalar ipr_a,
                 const Scalar ipr_b,
                 const std::function<Scalar(const Scalar)>& adjust_bhp)
{
    // Given fixed thp, wfr, gfr and alq, this function finds a stable (-flo, bhp)-intersection
    // between the ipr-line and bhp(flo) from table, if such an intersection exists. For multiple 
    // stable intersections, the one corresponding the largest flo is returned as long as this intersection
    // lies within the tabulated values. If the ipr-line lies above all (flo, bhp) points, the intersection
    // is determined by extrapolation based on the last two points. 
    // The adjust_bhp-function is used to adjust the vfp-table bhp-values to actual bhp-values due
    // to vfp/well ref-depth differences and/or WVFPDP-related pressure adjustments.

    // NOTE: ipr-line is q=b*bhp - a!
    // ipr is given for negative flo, so
    // flo = -b*bhp + a, i.e., bhp = -(flo-a)/b
    auto thp_i = findInterpData( thp, table.getTHPAxis());
    auto wfr_i = findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = findInterpData( gfr, table.getGFRAxis());
    auto alq_i = findInterpData( alq, table.getALQAxis());

    if (ipr_b == 0.0) {
        // this shouldn't happen, but deal with it to be safe
        auto flo_i = findInterpData(ipr_a, table.getFloAxis());
        detail::VFPEvaluation bhp_i = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        return std::make_pair(-ipr_a, adjust_bhp(bhp_i.value));
    }
    // find largest flo (flo_x) for which y = bhp(flo) + (flo-a)/b = 0 and dy/dflo > 0
    Scalar flo_x = -1.0;
    Scalar flo0;
    Scalar y0, y1;
    flo0 = 0.0; // start by checking flo=0
    auto flo_i = findInterpData(flo0, table.getFloAxis());
    detail::VFPEvaluation bhp_i = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
    y0 = adjust_bhp(bhp_i.value) - ipr_a/ipr_b; // +0.0/ipr_b

    const std::vector<double>& flos = table.getFloAxis();
    for (size_t i = 0; i < flos.size(); ++i) {
        const auto flo1 = flos[i];
        flo_i = findInterpData(flo1, flos);
        bhp_i = interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        y1 = adjust_bhp(bhp_i.value) + (flo1 - ipr_a)/ipr_b;
        if (y0 < 0 && y1 >= 0){
            // crossing with positive slope
            Scalar w = -y0/(y1-y0);
            w = std::clamp(w, Scalar{0.0}, Scalar{1.0}); // just to be safe (if y0~y1~0)
            flo_x = flo0 + w*(flo1 - flo0);
        }
        if (i < flos.size()-1) { // check next interval
            flo0 = flo1;
            y0 = y1;
        } else if (y1 < 0 && y0 < y1 && flo_x < 0) { // at last interval
            // If y0 < y1 < 0, there is a stable intersection above the largest flo-value by 
            // extrapolation. If no previous stable intersections were found, i.e., ipr-line lies 
            // above all (flo, bhp) points, then we return this intersection. Otherwise, we don't 
            // trust it (avoid vfp-extrapolation whenever possible)
            Scalar w = -y0/(y1-y0); // w > 1.0
            flo_x = flo0 + w*(flo1 - flo0);
        }
    }
    // return (last) intersection if found (negative flo)
    if (flo_x >= 0) {
        return std::make_pair(-flo_x, -(flo_x - ipr_a)/ipr_b);
    } else {
        return std::nullopt;
    }
}

namespace detail {

template<class Scalar>
VFPEvaluation<Scalar> operator+(VFPEvaluation<Scalar> lhs, const VFPEvaluation<Scalar>& rhs)
{
    lhs.value += rhs.value;
    lhs.dthp += rhs.dthp;
    lhs.dwfr += rhs.dwfr;
    lhs.dgfr += rhs.dgfr;
    lhs.dalq += rhs.dalq;
    lhs.dflo += rhs.dflo;
    return lhs;
}

template<class Scalar>
VFPEvaluation<Scalar> operator-(VFPEvaluation<Scalar> lhs, const VFPEvaluation<Scalar>& rhs)
{
    lhs.value -= rhs.value;
    lhs.dthp -= rhs.dthp;
    lhs.dwfr -= rhs.dwfr;
    lhs.dgfr -= rhs.dgfr;
    lhs.dalq -= rhs.dalq;
    lhs.dflo -= rhs.dflo;
    return lhs;
}

template<class Scalar>
VFPEvaluation<Scalar> operator*(Scalar lhs, const VFPEvaluation<Scalar>& rhs)
{
    VFPEvaluation<Scalar> retval;
    retval.value = rhs.value * lhs;
    retval.dthp = rhs.dthp * lhs;
    retval.dwfr = rhs.dwfr * lhs;
    retval.dgfr = rhs.dgfr * lhs;
    retval.dalq = rhs.dalq * lhs;
    retval.dflo = rhs.dflo * lhs;
    return retval;
}

template <typename T>
T getFlo(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour)
{
    auto type = table.getFloType();
    switch (type) {
    case VFPProdTable::FLO_TYPE::FLO_OIL:
        //Oil = liquid phase
        return liquid;
    case VFPProdTable::FLO_TYPE::FLO_LIQ:
        //Liquid = aqua + liquid phases
        return aqua + liquid;
    case VFPProdTable::FLO_TYPE::FLO_GAS:
        //Gas = vapor phase
        return vapour;
    default:
        throw std::logic_error("Invalid FLO_TYPE");
    }
}

template <typename T>
T getFlo(const VFPInjTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour)
{
    auto type = table.getFloType();
    switch (type) {
    case VFPInjTable::FLO_TYPE::FLO_OIL:
        //Oil = liquid phase
        return liquid;
    case VFPInjTable::FLO_TYPE::FLO_WAT:
        //Liquid = aqua phase
        return aqua;
    case VFPInjTable::FLO_TYPE::FLO_GAS:
        //Gas = vapor phase
        return vapour;
    default:
        throw std::logic_error("Invalid FLO_TYPE");
    }
}

static constexpr double threshold = 1e-12;

template <typename T>
T getWFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour)
{
    auto type = table.getWFRType();
    switch(type) {
    case VFPProdTable::WFR_TYPE::WFR_WOR: {
        //Water-oil ratio = water / oil
        return chopNegativeValues(-aqua) / max(threshold, chopNegativeValues(-liquid));
    }
    case VFPProdTable::WFR_TYPE::WFR_WCT:
        //Water cut = water / (water + oil)
        return chopNegativeValues(-aqua) / max(threshold, chopNegativeValues(-aqua - liquid));
    case VFPProdTable::WFR_TYPE::WFR_WGR:
        //Water-gas ratio = water / gas
        return chopNegativeValues(-aqua) / max(threshold, chopNegativeValues(-vapour));
    default:
        throw std::logic_error("Invalid WFR_TYPE");
    }
}

template <typename T>
T getGFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour)
{
    auto type = table.getGFRType();
    switch(type) {
    case VFPProdTable::GFR_TYPE::GFR_GOR:
        // Gas-oil ratio = gas / oil
        return chopNegativeValues(-vapour) / max(threshold, chopNegativeValues(-liquid));
    case VFPProdTable::GFR_TYPE::GFR_GLR:
        // Gas-liquid ratio = gas / (oil + water)
        return chopNegativeValues(-vapour) / max(threshold, chopNegativeValues(-liquid - aqua));
    case VFPProdTable::GFR_TYPE::GFR_OGR:
        // Oil-gas ratio = oil / gas
        return chopNegativeValues(-liquid) / max(threshold, chopNegativeValues(-vapour));
    default:
        throw std::logic_error("Invalid GFR_TYPE");
    }
}

template <typename T>
const T& getTable(const std::map<int, std::reference_wrapper<const T>>& tables, int table_id)
{
    auto entry = tables.find(table_id);
    if (entry == tables.end()) {
        OPM_THROW(std::invalid_argument,
                  "Nonexistent VFP table " +
                  std::to_string(table_id) + " referenced.");
    }
    else {
        return entry->second.get();
    }
}

template <>
VFPProdTable::FLO_TYPE getType(const VFPProdTable& table)
{
    return table.getFloType();
}

template <>
VFPProdTable::WFR_TYPE getType(const VFPProdTable& table)
{
    return table.getWFRType();
}

template <>
VFPProdTable::GFR_TYPE getType(const VFPProdTable& table)
{
    return table.getGFRType();
}

/**
 * Returns the type variable for FLO for injection tables
 */
template <>
VFPInjTable::FLO_TYPE getType(const VFPInjTable& table)
{
    return table.getFloType();
}

template const VFPInjTable&
getTable(const std::map<int, std::reference_wrapper<const VFPInjTable>>&, int);
template const VFPProdTable&
getTable(const std::map<int, std::reference_wrapper<const VFPProdTable>>&, int);

#define INSTANTIATE(...)                            \
    template __VA_ARGS__                            \
    getFlo(const VFPInjTable&, const __VA_ARGS__&,  \
           const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__                            \
    getFlo(const VFPProdTable&, const __VA_ARGS__&, \
           const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__                            \
    getGFR(const VFPProdTable&, const __VA_ARGS__&, \
           const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__                            \
    getWFR(const VFPProdTable&, const __VA_ARGS__&, \
           const __VA_ARGS__&, const __VA_ARGS__&);

#define INSTANTIATE_TYPE(T)                      \
    INSTANTIATE(T)                               \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 4u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 5u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 6u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 7u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 8u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 9u>)  \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 10u>) \
    INSTANTIATE(DenseAd::Evaluation<T, -1, 11u>) \
    INSTANTIATE(DenseAd::Evaluation<T, 3, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 4, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 5, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 6, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 7, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 8, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 9, 0u>)   \
    INSTANTIATE(DenseAd::Evaluation<T, 10, 0u>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace detail

template class VFPHelpers<double>;

#if FLOW_INSTANTIATE_FLOAT
template class VFPHelpers<float>;
#endif

} // namespace Opm
