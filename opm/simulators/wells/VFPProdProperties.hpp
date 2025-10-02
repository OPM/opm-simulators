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

#ifndef OPM_AUTODIFF_VFPPRODPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPPRODPROPERTIES_HPP_

#include <functional>
#include <map>
#include <vector>


namespace Opm {

class VFPProdTable;

/**
 * Class which linearly interpolates BHP as a function of rate, tubing head pressure,
 * water fraction, gas fraction, and artificial lift for production VFP tables, and similarly
 * the BHP as a function of the rate and tubing head pressure.
 */
template<class Scalar>
class VFPProdProperties {
public:
    /**
     * Takes *no* ownership of data.
     */
    void addTable(const VFPProdTable& new_table);

    /**
     * Linear interpolation of bhp as a function of the input parameters given as
     * Evalutions
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     * @param explicit_wfr Explicit wfr
     * @param explicit_gfr Explicit gfr
     * @param use_expvfp True to use explicit VFP calculations
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    template <class EvalWell>
    EvalWell bhp(const int       table_id,
                 const EvalWell& aqua,
                 const EvalWell& liquid,
                 const EvalWell& vapour,
                 const Scalar    thp,
                 const Scalar    alq,
                 const Scalar    explicit_wfr,
                 const Scalar    explicit_gfr,
                 const bool      use_expvfp) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     * @param explicit_wfr Explicit wfr
     * @param explicit_gfr Explicit gfr
     * @param use_expvfp True to use explicit VFP calculations
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    Scalar bhp(const int    table_id,
               const Scalar aqua,
               const Scalar liquid,
               const Scalar vapour,
               const Scalar thp,
               const Scalar alq,
               const Scalar explicit_wfr,
               const Scalar explicit_gfr,
               const bool   use_expvfp) const;

    /**
     * Linear interpolation of thp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param bhp Bottom hole pressure
     * @param alq Artificial lift or other parameter
     * @param explicit_wfr Explicit wfr
     * @param explicit_gfr Explicit gfr
     * @param use_expvfp True to use explicit VFP calculations
     *
     * @return The tubing hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    Scalar thp(const int table_id,
               const Scalar aqua,
               const Scalar liquid,
               const Scalar vapour,
               const Scalar bhp,
               const Scalar alq,
               const Scalar explicit_wfr,
               const Scalar explicit_gfr,
               const bool use_expvfp) const;

    /**
     * Returns the table associated with the ID, or throws an exception if
     * the table does not exist
     */
    const VFPProdTable& getTable(const int table_id) const;

    /**
     * Check whether there is table associated with ID
     */
    bool hasTable(const int table_id) const;

    /**
     * Returns true if no vfp tables are in the current map
     */
    bool empty() const
    {
        return m_tables.empty();
    }

    /**
     * Returns minimum bhp for given thp, wfr, gfr and alq
     */
    Scalar minimumBHP(const int table_id, const Scalar thp,
                      const Scalar wfr, const Scalar gfr, const Scalar alq) const;

protected:
    // calculate a group bhp values with a group of flo rate values
    std::vector<Scalar> bhpwithflo(const std::vector<Scalar>& flos,
                                   const int table_id,
                                   const Scalar wfr,
                                   const Scalar gfr,
                                   const Scalar thp,
                                   const Scalar alq,
                                   const Scalar dp) const;

    // Map which connects the table number with the table itself
    std::map<int, std::reference_wrapper<const VFPProdTable>> m_tables;
};

} // namespace Opm

#endif /* OPM_AUTODIFF_VFPPRODPROPERTIES_HPP_ */
