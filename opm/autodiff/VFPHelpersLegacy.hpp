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


#ifndef OPM_AUTODIFF_VFPHELPERSLEGACY_HPP_
#define OPM_AUTODIFF_VFPHELPERSLEGACY_HPP_


#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

/**
 * This file contains a set of helper functions used by VFPProd / VFPInj.
 */
namespace Opm {
namespace detail {


typedef AutoDiffBlock<double> ADB;


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

}
}

#include "VFPHelpers.hpp"


namespace Opm {
namespace detail {

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

} // namespace detail


} // namespace


#endif /* OPM_AUTODIFF_VFPHELPERSLEGACY_HPP_ */
