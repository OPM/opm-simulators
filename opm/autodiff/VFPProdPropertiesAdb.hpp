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

#ifndef OPM_AUTODIFF_VFPPRODPROPERTIES_ADB_HPP_
#define OPM_AUTODIFF_VFPPRODPROPERTIES_ADB_HPP_

#include <opm/parser/eclipse/EclipseState/Tables/VFPProdTable.hpp>
#include <opm/core/wells.h>
#include <opm/autodiff/VFPHelpersAdb.hpp>

#include <vector>
#include <map>


namespace Opm {

template <class Scalar>
class AutoDiffBlock;

/**
 * Class which linearly interpolates BHP as a function of rate, tubing head pressure,
 * water fraction, gas fraction, and artificial lift for production VFP tables, and similarly
 * the BHP as a function of the rate and tubing head pressure.
 */
class VFPProdPropertiesAdb {
public:
    typedef AutoDiffBlock<double> ADB;

    /**
     * Empty constructor
     */
    VFPProdPropertiesAdb();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_table A *single* VFPPROD table
     */
    explicit VFPProdPropertiesAdb(const VFPProdTable* prod_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_tables A map of different VFPPROD tables.
     */
    explicit VFPProdPropertiesAdb(const std::map<int, VFPProdTable>& prod_tables);

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
    double bhp(int table_id,
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
     * Returns the table associated with the ID, or throws an exception if
     * the table does not exist
     */
    const VFPProdTable* getTable(const int table_id) const;

    /**
     * Returns true if no vfp tables are in the current map
     */
    bool empty() const {
        return m_tables.empty();
    }

private:
    // Map which connects the table number with the table itself
    std::map<int, const VFPProdTable*> m_tables;
};




} //namespace


#endif /* OPM_AUTODIFF_VFPPRODPROPERTIES_ADB_HPP_ */
