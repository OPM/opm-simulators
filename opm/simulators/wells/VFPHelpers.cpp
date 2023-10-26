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
double findX(const double  x0,
             const double  x1,
             const double  y0,
             const double  y1,
             const double  y)
{
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
 * Returns zero if input value is negative
 */
template <typename T>
static T chopNegativeValues(const T& value) {
    return Opm::max(0.0, value);
}

}

namespace Opm {
namespace detail {

InterpData findInterpData(const double value_in, const std::vector<double>& values)
{
    InterpData retval;

    const int nvalues = values.size();

    // chopping the value to be zero, which means we do not
    // extrapolate the table towards nagative ranges
    const double value = value_in < 0.? 0. : value_in;

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

VFPEvaluation operator+(VFPEvaluation lhs, const VFPEvaluation& rhs)
{
    lhs.value += rhs.value;
    lhs.dthp += rhs.dthp;
    lhs.dwfr += rhs.dwfr;
    lhs.dgfr += rhs.dgfr;
    lhs.dalq += rhs.dalq;
    lhs.dflo += rhs.dflo;
    return lhs;
}

VFPEvaluation operator-(VFPEvaluation lhs, const VFPEvaluation& rhs)
{
    lhs.value -= rhs.value;
    lhs.dthp -= rhs.dthp;
    lhs.dwfr -= rhs.dwfr;
    lhs.dgfr -= rhs.dgfr;
    lhs.dalq -= rhs.dalq;
    lhs.dflo -= rhs.dflo;
    return lhs;
}

VFPEvaluation operator*(double lhs, const VFPEvaluation& rhs)
{
    VFPEvaluation retval;
    retval.value = rhs.value * lhs;
    retval.dthp = rhs.dthp * lhs;
    retval.dwfr = rhs.dwfr * lhs;
    retval.dgfr = rhs.dgfr * lhs;
    retval.dalq = rhs.dalq * lhs;
    retval.dflo = rhs.dflo * lhs;
    return retval;
}

VFPEvaluation interpolate(const VFPProdTable& table,
                          const InterpData& flo_i,
                          const InterpData& thp_i,
                          const InterpData& wfr_i,
                          const InterpData& gfr_i,
                          const InterpData& alq_i)
{
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

VFPEvaluation interpolate(const VFPInjTable& table,
                          const InterpData& flo_i,
                          const InterpData& thp_i)
{
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

VFPEvaluation bhp(const VFPProdTable& table,
                  const double  aqua,
                  const double  liquid,
                  const double  vapour,
                  const double  thp,
                  const double  alq,
                  const double  explicit_wfr,
                  const double  explicit_gfr,
                  const bool    use_vfpexplicit)
{
    //Find interpolation variables
    double flo = detail::getFlo(table, aqua, liquid, vapour);
    double wfr = detail::getWFR(table, aqua, liquid, vapour);
    double gfr = detail::getGFR(table, aqua, liquid, vapour);
    if (use_vfpexplicit || -flo < table.getFloAxis().front()) {
        wfr = explicit_wfr;
        gfr = explicit_gfr;
    }

    //First, find the values to interpolate between
    //Recall that flo is negative in Opm, so switch sign.
    auto flo_i = detail::findInterpData(-flo, table.getFloAxis());
    auto thp_i = detail::findInterpData( thp, table.getTHPAxis());
    auto wfr_i = detail::findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table.getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table.getALQAxis());

    detail::VFPEvaluation retval = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);

    return retval;
}

VFPEvaluation bhp(const VFPInjTable& table,
                  const double  aqua,
                  const double  liquid,
                  const double  vapour,
                  const double  thp)
{
    //Find interpolation variables
    double flo = detail::getFlo(table, aqua, liquid, vapour);

    //First, find the values to interpolate between
    auto flo_i = detail::findInterpData(flo, table.getFloAxis());
    auto thp_i = detail::findInterpData(thp, table.getTHPAxis());

    //Then perform the interpolation itself
    detail::VFPEvaluation retval = detail::interpolate(table, flo_i, thp_i);

    return retval;
}

double findTHP(const std::vector<double>& bhp_array,
               const std::vector<double>& thp_array,
               double bhp)
{
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
            thp = findX(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
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
            thp = findX(x0, x1, y0, y1, bhp);
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
            thp = findX(x0, x1, y0, y1, bhp);
        }
        else if (bhp <= bhp_array[0]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[0];
            const double& x1 = thp_array[1];
            const double& y0 = bhp_array[0];
            const double& y1 = bhp_array[1];
            thp = findX(x0, x1, y0, y1, bhp);
        }
        //Target bhp greater than all values in array, extrapolate
        else if (bhp > bhp_array[nthp-1]) {
            //TODO: LOG extrapolation
            const double& x0 = thp_array[nthp-2];
            const double& x1 = thp_array[nthp-1];
            const double& y0 = bhp_array[nthp-2];
            const double& y1 = bhp_array[nthp-1];
            thp = findX(x0, x1, y0, y1, bhp);
        }
        else {
            OPM_THROW(std::logic_error, "Programmer error: Unable to find THP in THP array");
        }
    }

    return thp;
}

std::pair<double, double> 
getMinimumBHPCoordinate(const VFPProdTable& table,
                        const double thp,
                        const double wfr,
                        const double gfr,
                        const double alq)
{   
    double flo_at_bhp_min = 0.0; // start by checking flo=0
    auto flo_i = detail::findInterpData(flo_at_bhp_min, table.getFloAxis());
    auto thp_i = detail::findInterpData( thp, table.getTHPAxis());
    auto wfr_i = detail::findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table.getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table.getALQAxis());

    detail::VFPEvaluation bhp_i = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
    double bhp_min = bhp_i.value;
    std::vector<double> flos = table.getFloAxis();
    for (size_t i = 0; i < flos.size(); ++i) {
        flo_i = detail::findInterpData(flos[i], table.getFloAxis());
        bhp_i = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        if (bhp_i.value < bhp_min){
            bhp_min = bhp_i.value;
            flo_at_bhp_min = flos[i];
        }
    }
    // return negative flo
    return std::make_pair(-flo_at_bhp_min, bhp_min); 
}      

std::optional<std::pair<double, double>> 
intersectWithIPR(const VFPProdTable& table,
                 const double thp,
                 const double wfr,
                 const double gfr,
                 const double alq, 
                 const double ipr_a,
                 const double ipr_b)
{   
    // NOTE: ipr-line is q=b*bhp - a!
    // ipr is given for negative flo, so
    // flo = -b*bhp + a, i.e., bhp = -(flo-a)/b  
    auto thp_i = detail::findInterpData( thp, table.getTHPAxis());
    auto wfr_i = detail::findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table.getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table.getALQAxis());

    if (ipr_b == 0.0) {
        // this shouldn't happen, but deal with it to be safe
        auto flo_i = detail::findInterpData(-ipr_a, table.getFloAxis());
        detail::VFPEvaluation bhp_i = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        return std::make_pair(ipr_b, bhp_i.value);
    }
    // find largest flo (flo_x) for which y = bhp(flo) + (flo-a)/b = 0 and dy/dflo > 0
    double flo_x = -1.0;
    double flo0, flo1;
    double y0, y1;
    flo0 = 0.0; // start by checking flo=0
    auto flo_i = detail::findInterpData(flo0, table.getFloAxis());
    detail::VFPEvaluation bhp_i = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
    y0 = bhp_i.value - ipr_a/ipr_b; // +0.0/ipr_b

    std::vector<double> flos = table.getFloAxis();
    for (size_t i = 0; i < flos.size(); ++i) {
        flo1 = flos[i];
        flo_i = detail::findInterpData(flo1, table.getFloAxis());
        bhp_i = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);
        y1 = bhp_i.value + (flo1 - ipr_a)/ipr_b;
        if (y0 < 0 && y1 >= 0){
            // crossing with positive slope
            double w = -y0/(y1-y0);
            w = std::clamp(w, 0.0, 1.0); // just to be safe (if y0~y1)
            flo_x = flo0 + w*(flo1 - flo0);
        }
        flo0 = flo1;
        y0 = y1;
    }
    // return (last) intersection if found (negative flo)
    if (flo_x >= 0) {
        return std::make_pair(-flo_x, -(flo_x - ipr_a)/ipr_b);
    } else {
        return std::nullopt;
    }
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

template const VFPInjTable& getTable<VFPInjTable>(const std::map<int, std::reference_wrapper<const VFPInjTable>>&, int);
template const VFPProdTable& getTable<VFPProdTable>(const std::map<int, std::reference_wrapper<const VFPProdTable>>&, int);

#define INSTANCE(...) \
    template __VA_ARGS__ getFlo(const VFPInjTable&, const __VA_ARGS__&, const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__ getFlo(const VFPProdTable&, const __VA_ARGS__&, const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__ getGFR(const VFPProdTable&, const __VA_ARGS__&, const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__ getWFR(const VFPProdTable&, const __VA_ARGS__&, const __VA_ARGS__&, const __VA_ARGS__&);

INSTANCE(double)
INSTANCE(DenseAd::Evaluation<double, -1, 4u>)
INSTANCE(DenseAd::Evaluation<double, -1, 5u>)
INSTANCE(DenseAd::Evaluation<double, -1, 6u>)
INSTANCE(DenseAd::Evaluation<double, -1, 7u>)
INSTANCE(DenseAd::Evaluation<double, -1, 8u>)
INSTANCE(DenseAd::Evaluation<double, -1, 9u>)
INSTANCE(DenseAd::Evaluation<double, -1, 10u>)
INSTANCE(DenseAd::Evaluation<double, -1, 11u>)
INSTANCE(DenseAd::Evaluation<double, 3, 0u>)
INSTANCE(DenseAd::Evaluation<double, 4, 0u>)
INSTANCE(DenseAd::Evaluation<double, 5, 0u>)
INSTANCE(DenseAd::Evaluation<double, 6, 0u>)
INSTANCE(DenseAd::Evaluation<double, 7, 0u>)
INSTANCE(DenseAd::Evaluation<double, 8, 0u>)
INSTANCE(DenseAd::Evaluation<double, 9, 0u>)
INSTANCE(DenseAd::Evaluation<double, 10, 0u>)

} // namespace detail
} // namespace Opm
