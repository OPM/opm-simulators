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

#include <array>
#include <functional>
#include <map>
#include <vector>
#include <optional>

/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {

class VFPInjTable;
class VFPProdTable;

namespace detail {

/**
 * Computes the flo parameter according to the flo_type_
 * for production tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getFlo(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the flo parameter according to the flo_type_
 * for injection tables
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getFlo(const VFPInjTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the wfr parameter according to the wfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getWFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

/**
 * Computes the gfr parameter according to the gfr_type_
 * @return Production rate of oil, gas or liquid.
 */
template <typename T>
T getGFR(const VFPProdTable& table,
         const T& aqua,
         const T& liquid,
         const T& vapour);

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
 *  @param value_in Value to find in values
 *  @param values Sorted list of values to search for value in.
 *  @return Data required to find the interpolated value
 */
InterpData findInterpData(const double value_in, const std::vector<double>& values);

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

VFPEvaluation operator+(VFPEvaluation lhs, const VFPEvaluation& rhs);
VFPEvaluation operator-(VFPEvaluation lhs, const VFPEvaluation& rhs);
VFPEvaluation operator*(double lhs, const VFPEvaluation& rhs);

/**
 * Helper function which interpolates data using the indices etc. given in the inputs.
 */
VFPEvaluation interpolate(const VFPProdTable& table,
                          const InterpData& flo_i,
                          const InterpData& thp_i,
                          const InterpData& wfr_i,
                          const InterpData& gfr_i,
                          const InterpData& alq_i);

/**
 * This basically models interpolate(VFPProdTable::array_type, ...)
 * which performs 5D interpolation, but here for the 2D case only
 */
VFPEvaluation interpolate(const VFPInjTable& table,
                          const InterpData& flo_i,
                          const InterpData& thp_i);

VFPEvaluation bhp(const VFPProdTable& table,
                  const double aqua,
                  const double liquid,
                  const double vapour,
                  const double thp,
                  const double alq,
                  const double explicit_wfr,
                  const double explicit_gfr,
                  const bool   use_vfpexplicit);

VFPEvaluation bhp(const VFPInjTable& table,
                  const double aqua,
                  const double liquid,
                  const double vapour,
                  const double thp);


/**
 * Returns the table from the map if found, or throws an exception
 */
template <typename T>
const T& getTable(const std::map<int, std::reference_wrapper<const T>>& tables, int table_id);

/**
 * Check whether we have a table with the table number
 */
template <typename T>
bool hasTable(const std::map<int, std::reference_wrapper<const T>>& tables, int table_id) {
    const auto entry = tables.find(table_id);
    return (entry != tables.end() );
}


/**
 * Returns the type variable for FLO/GFR/WFR for production tables
 */
template <typename TYPE, typename TABLE>
TYPE getType(const TABLE& table);


/**
 * This function finds the value of THP given a specific BHP.
 * Essentially:
 *   Given the function f(thp_array(x)) = bhp_array(x), which is piecewise linear,
 *   find thp so that f(thp) = bhp.
 */
double findTHP(const std::vector<double>& bhp_array,
               const std::vector<double>& thp_array,
               double bhp);

/**
* Get (flo, bhp) at minimum bhp for given thp,wfr,gfr,alq
*/
std::pair<double, double> 
getMinimumBHPCoordinate(const VFPProdTable& table,
                        const double thp,
                        const double wfr,
                        const double gfr,
                        const double alq);

/**
* Get (flo, bhp) at largest occuring stable vfp/ipr-intersection
* if it exists
*/  
std::optional<std::pair<double, double>> 
intersectWithIPR(const VFPProdTable& table,
                 const double thp,
                 const double wfr,
                 const double gfr,
                 const double alq, 
                 const double ipr_a,
                 const double ipr_b);                        

} // namespace detail


} // namespace




#endif /* OPM_AUTODIFF_VFPHELPERS_HPP_ */
