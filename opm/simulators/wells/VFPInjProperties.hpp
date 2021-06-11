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

#ifndef OPM_AUTODIFF_VFPINJPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPINJPROPERTIES_HPP_


#include <functional>
#include <map>
#include <vector>


namespace Opm {

class VFPInjTable;

class VFPInjProperties {
public:
    VFPInjProperties() = default;
    /**
     * Takes *no* ownership of data.
     */
    void addTable(const VFPInjTable& new_table);

    /**
     * Linear interpolation of bhp as a function of the input parameters given as
     * Evaluation
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    template <class EvalWell>
    EvalWell bhp(const int table_id,
                 const EvalWell& aqua,
                 const EvalWell& liquid,
                 const EvalWell& vapour,
                 const double& thp) const;

    /**
     * Returns the table associated with the ID, or throws an exception if
     * the table does not exist
     */
    const VFPInjTable& getTable(const int table_id) const;

    /**
     * Check whether there is table associated with ID
     */
    bool hasTable(const int table_id) const;

    /**
     * Returns true if no vfp tables are in the current map
     */
    bool empty() const {
        return m_tables.empty();
    }

    /**
     * Linear interpolation of bhp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double bhp(int table_id,
               const double& aqua,
               const double& liquid,
               const double& vapour,
               const double& thp) const;

    /**
     * Linear interpolation of thp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param bhp Bottom hole pressure
     *
     * @return The tubing hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double thp(int table_id,
               const double& aqua,
               const double& liquid,
               const double& vapour,
               const double& bhp) const;

protected:
    // Map which connects the table number with the table itself
    std::map<int, std::reference_wrapper<const VFPInjTable>> m_tables;
};



} //namespace



#endif /* OPM_AUTODIFF_VFPINJPROPERTIES_HPP_ */
