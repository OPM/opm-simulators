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


/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {
namespace detail {


typedef VFPProdProperties::ADB ADB;










/**
 * Returns zero if input value is NaN
 */
inline double zeroIfNan(const double& value) {
    return (std::isnan(value)) ? 0.0 : value;
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
 * Helper function to find indices etc. for linear interpolation
 */
inline InterpData findInterpData(const double& value, const std::vector<double>& values) {
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







/**
 * Helper function which interpolates data using the indices etc. given in the inputs.
 */
#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif
inline VFPProdProperties::adb_like interpolate(
        const VFPProdTable::array_type& array,
        const InterpData& flo_i,
        const InterpData& thp_i,
        const InterpData& wfr_i,
        const InterpData& gfr_i,
        const InterpData& alq_i) {

    //Values and derivatives in a 5D hypercube
    VFPProdProperties::adb_like nn[2][2][2][2][2];


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





} // namespace detail


} // namespace




#endif /* OPM_AUTODIFF_VFPHELPERS_HPP_ */
