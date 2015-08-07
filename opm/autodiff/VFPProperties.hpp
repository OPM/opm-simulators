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

#ifndef OPM_AUTODIFF_VFPPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPPROPERTIES_HPP_

#include <opm/core/wells.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/VFPInjTable.hpp>
#include <boost/multi_array.hpp>

#include <map>

namespace Opm {

class VFPProdProperties;
class VFPInjProperties;

/**
 * A thin wrapper class that holds one VFPProdProperties and one
 * VFPInjProperties object.
 */
class VFPProperties {
public:
    VFPProperties();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_table  A *single* VFPINJ table or NULL (no table)
     * @param prod_table A *single* VFPPROD table or NULL (no table)
     */
    VFPProperties(const VFPInjTable* inj_table, const VFPProdTable* prod_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     * @param prod_tables A map of different VFPPROD tables.
     */
    VFPProperties(const std::map<int, VFPInjTable>& inj_tables,
                  const std::map<int, VFPProdTable>& prod_tables);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     */
    VFPProperties(const std::map<int, VFPInjTable>& inj_tables);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_tables A map of different VFPPROD tables.
     */
    VFPProperties(const std::map<int, VFPProdTable>& prod_tables);

    const VFPInjProperties* getInj() const {
        return m_inj.get();
    }

    const VFPProdProperties* getProd() const {
        return m_prod.get();
    }

private:
    std::shared_ptr<VFPInjProperties> m_inj;
    std::shared_ptr<VFPProdProperties> m_prod;
};

/**
 * Class which linearly interpolates BHP as a function of rate, tubing head pressure,
 * water fraction, gas fraction, and artificial lift for production VFP tables, and similarly
 * the BHP as a function of the rate and tubing head pressure.
 */
class VFPProdProperties {
public:
    typedef AutoDiffBlock<double> ADB;

    /**
     * An "ADB-like" structure with a single value and a set of derivatives
     */
    struct adb_like {
        adb_like() : value(0.0), dthp(0.0), dwfr(0.0), dgfr(0.0), dalq(0.0), dflo(0.0) {};
        double value;
        double dthp;
        double dwfr;
        double dgfr;
        double dalq;
        double dflo;
    };

    /**
     * Empty constructor
     */
    VFPProdProperties();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_table A *single* VFPPROD table
     */
    VFPProdProperties(const VFPProdTable* prod_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_tables A map of different VFPPROD tables.
     */
    VFPProdProperties(const std::map<int, VFPProdTable>& prod_tables);

    /**
     * Linear interpolation of bhp as function of the input parameters.
     * @param table_id Table number to use
     * @param wells Wells structure with information about wells in qs
     * @param qs Flow quantities
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    ADB bhp(const std::vector<int>& table_id,
            const Wells& wells,
            const ADB& qs,
            const ADB& thp,
            const ADB& alq) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters given as ADBs
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    ADB bhp(const std::vector<int>& table_id,
            const ADB& aqua,
            const ADB& liquid,
            const ADB& vapour,
            const ADB& thp,
            const ADB& alq) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    adb_like bhp(int table_id,
            const double& aqua,
            const double& liquid,
            const double& vapour,
            const double& thp,
            const double& alq) const;

    /**
     * Linear interpolation of thp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param bhp Bottom hole pressure
     * @param alq Artificial lift or other parameter
     *
     * @return The tubing hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double thp(int table_id,
            const double& aqua,
            const double& liquid,
            const double& vapour,
            const double& bhp,
            const double& alq) const;

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


private:
    // Map which connects the table number with the table itself
    std::map<int, const VFPProdTable*> m_tables;

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
    static InterpData find_interp_data(const double& value, const std::vector<double>& values);

    /**
     * Helper function which interpolates data using the indices etc. given in the inputs.
     */
    static adb_like interpolate(const VFPProdTable::array_type& array,
            const InterpData& flo_i,
            const InterpData& thp_i,
            const InterpData& wfr_i,
            const InterpData& gfr_i,
            const InterpData& alq_i);

    /**
     * Helper function that finds x for a given value of y for a line
     * *NOTE ORDER OF ARGUMENTS*
     */
    static double find_x(const double& x0,
            const double& x1,
            const double& y0,
            const double& y1,
            const double& y);


    /**
     * Initialization routine
     */
    void init(const std::map<int, VFPProdTable>& prod_tables);

    /**
     * Misc helper functions
     */
    const VFPProdTable* getProdTable(int table_id) const;

    static inline double zeroIfNan(const double& value) {
        return (std::isnan(value)) ? 0.0 : value;
    }

    static inline ADB zeroIfNan(const ADB& values) {
        Selector<ADB::V::Scalar> not_nan_selector(values.value(), Selector<ADB::V::Scalar>::NotNaN);

        const ADB::V z = ADB::V::Zero(values.size());
        const ADB zero = ADB::constant(z, values.blockPattern());

        ADB retval = not_nan_selector.select(values, zero);
        return retval;
    }
};



inline VFPProdProperties::adb_like operator+(
        VFPProdProperties::adb_like lhs,
        const VFPProdProperties::adb_like& rhs) {
    lhs.value += rhs.value;
    lhs.dthp += rhs.dthp;
    lhs.dwfr += rhs.dwfr;
    lhs.dgfr += rhs.dgfr;
    lhs.dalq += rhs.dalq;
    lhs.dflo += rhs.dflo;
    return lhs;
}

inline VFPProdProperties::adb_like operator-(
        VFPProdProperties::adb_like lhs,
        const VFPProdProperties::adb_like& rhs) {
    lhs.value -= rhs.value;
    lhs.dthp -= rhs.dthp;
    lhs.dwfr -= rhs.dwfr;
    lhs.dgfr -= rhs.dgfr;
    lhs.dalq -= rhs.dalq;
    lhs.dflo -= rhs.dflo;
    return lhs;
}

inline VFPProdProperties::adb_like operator*(
        double lhs,
        const VFPProdProperties::adb_like& rhs) {
    VFPProdProperties::adb_like retval;
    retval.value = rhs.value * lhs;
    retval.dthp = rhs.dthp * lhs;
    retval.dwfr = rhs.dwfr * lhs;
    retval.dgfr = rhs.dgfr * lhs;
    retval.dalq = rhs.dalq * lhs;
    retval.dflo = rhs.dflo * lhs;
    return retval;
}

} //Namespace

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
