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

#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/core/wells.h>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>

#include <vector>
#include <map>


namespace Opm {

/**
 * Class which linearly interpolates BHP as a function of rate, tubing head pressure,
 * water fraction, gas fraction, and artificial lift for production VFP tables, and similarly
 * the BHP as a function of the rate and tubing head pressure.
 */
class VFPProdProperties {
public:
    /**
     * Empty constructor
     */
    VFPProdProperties();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_table A *single* VFPPROD table
     */
    explicit VFPProdProperties(const VFPProdTable* prod_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param prod_tables A map of different VFPPROD tables.
     */
    using ProdTable = std::map<int, std::shared_ptr<const VFPProdTable> >;
    explicit VFPProdProperties(const ProdTable& prod_tables);

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
                 const double& thp,
                 const double& alq) const {

        //Get the table
        const VFPProdTable* table = detail::getTable(m_tables, table_id);
        EvalWell bhp = 0.0 * aqua;

        //Find interpolation variables
        EvalWell flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());
        EvalWell wfr = detail::getWFR(aqua, liquid, vapour, table->getWFRType());
        EvalWell gfr = detail::getGFR(aqua, liquid, vapour, table->getGFRType());

        //Compute the BHP for each well independently
        if (table != nullptr) {
            //First, find the values to interpolate between
            //Value of FLO is negative in OPM for producers, but positive in VFP table
            auto flo_i = detail::findInterpData(-flo.value(), table->getFloAxis());
            auto thp_i = detail::findInterpData( thp, table->getTHPAxis()); // assume constant
            auto wfr_i = detail::findInterpData( wfr.value(), table->getWFRAxis());
            auto gfr_i = detail::findInterpData( gfr.value(), table->getGFRAxis());
            auto alq_i = detail::findInterpData( alq, table->getALQAxis()); //assume constant

            detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

            bhp = (bhp_val.dwfr * wfr) + (bhp_val.dgfr * gfr) - (bhp_val.dflo * flo);
            bhp.setValue(bhp_val.value);
        }
        else {
            bhp.setValue(-1e100); //Signal that this value has not been calculated properly, due to "missing" table
        }
        return bhp;
    }

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
     * Calculate the Bhp value from the THP target/constraint value
     * based on inflow performance relationship and VFP curves
     */
     double
     calculateBhpWithTHPTarget(const std::vector<double>& ipr_a,
                               const std::vector<double>& ipr_b,
                               const double bhp_limit,
                               const double thp_table_id,
                               const double thp_limit,
                               const double alq,
                               const double dp) const;

protected:
    // calculate a group bhp values with a group of flo rate values
    std::vector<double> bhpwithflo(const std::vector<double>& flos,
                                   const int table_id,
                                   const double wfr,
                                   const double gfr,
                                   const double thp,
                                   const double alq,
                                   const double dp) const;

    // Map which connects the table number with the table itself
    std::map<int, const VFPProdTable*> m_tables;
};




} //namespace


#endif /* OPM_AUTODIFF_VFPPRODPROPERTIES_HPP_ */
