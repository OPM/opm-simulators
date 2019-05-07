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

#include "config.h"

#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>



namespace Opm {




VFPProdProperties::VFPProdProperties() {

}


VFPProdProperties::VFPProdProperties(const VFPProdTable* table){
    m_tables[table->getTableNum()] = table;
}


VFPProdProperties::VFPProdProperties(const VFPProdProperties::ProdTable& tables) {
    for (const auto& table : tables) {
        m_tables[table.first] = table.second.get();
    }
}


double VFPProdProperties::thp(int table_id,
                              const double& aqua,
                              const double& liquid,
                              const double& vapour,
                              const double& bhp_arg,
                              const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);
    const VFPProdTable::array_type& data = table->getTable();

    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = detail::getWFR(aqua, liquid, vapour, table->getWFRType());
    double gfr = detail::getGFR(aqua, liquid, vapour, table->getGFRType());

    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     * Recall that flo is negative in Opm, so switch the sign
     */
    auto flo_i = detail::findInterpData(-flo, table->getFloAxis());
    auto wfr_i = detail::findInterpData( wfr, table->getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table->getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table->getALQAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(data, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}


double VFPProdProperties::bhp(int table_id,
                              const double& aqua,
                              const double& liquid,
                              const double& vapour,
                              const double& thp_arg,
                              const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg, alq);
    return retval.value;
}


const VFPProdTable* VFPProdProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}

bool VFPProdProperties::hasTable(const int table_id) const {
    return detail::hasTable(m_tables, table_id);
}


std::vector<double>
VFPProdProperties::
bhpwithflo(const std::vector<double>& flos,
           const int table_id,
           const double wfr,
           const double gfr,
           const double thp,
           const double alq,
           const double dp) const
{
    // Get the table
    const VFPProdTable* table = detail::getTable(m_tables, table_id);
    const auto thp_i = detail::findInterpData( thp, table->getTHPAxis()); // assume constant
    const auto wfr_i = detail::findInterpData( wfr, table->getWFRAxis());
    const auto gfr_i = detail::findInterpData( gfr, table->getGFRAxis());
    const auto alq_i = detail::findInterpData( alq, table->getALQAxis()); //assume constant

    std::vector<double> bhps(flos.size(), 0.);
    for (size_t i = 0; i < flos.size(); ++i) {
        // Value of FLO is negative in OPM for producers, but positive in VFP table
        const auto flo_i = detail::findInterpData(-flos[i], table->getFloAxis());
        const detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

        // TODO: this kind of breaks the conventions for the functions here by putting dp within the function
        bhps[i] = bhp_val.value - dp;
    }

    return bhps;
}





double
VFPProdProperties::
calculateBhpWithTHPTarget(const std::vector<double>& ipr_a,
                          const std::vector<double>& ipr_b,
                          const double bhp_limit,
                          const double thp_table_id,
                          const double thp_limit,
                          const double alq,
                          const double dp) const
{
    // For producers, bhp_safe_limit is the highest BHP value that can still produce based on IPR
    double bhp_safe_limit = 1.e100;
    for (size_t i = 0; i < ipr_a.size(); ++i) {
        if (ipr_b[i] == 0.) continue;

        const double bhp = ipr_a[i] / ipr_b[i];
        if (bhp < bhp_safe_limit) {
            bhp_safe_limit = bhp;
        }
    }

    // Here, we use the middle point between the bhp_limit and bhp_safe_limit to calculate the ratio of the flow
    // and the middle point serves one of the two points to describe inflow performance relationship line
    const double bhp_middle = (bhp_limit + bhp_safe_limit) / 2.0;

    // FLO is the rate based on the type specified with the VFP table
    // The two points correspond to the bhp values of bhp_limit, and the middle of bhp_limit and bhp_safe_limit
    // for producers, the rates are negative
    std::vector<double> rates_bhp_limit(ipr_a.size());
    std::vector<double> rates_bhp_middle(ipr_a.size());
    for (size_t i = 0; i < rates_bhp_limit.size(); ++i) {
        rates_bhp_limit[i] = bhp_limit * ipr_b[i] - ipr_a[i];
        rates_bhp_middle[i] = bhp_middle * ipr_b[i] - ipr_a[i];
    }

    // TODO: we need to be careful that there is nothings wrong related to the indices here
    const int Water = BlackoilPhases::Aqua;
    const int Oil = BlackoilPhases::Liquid;
    const int Gas = BlackoilPhases::Vapour;

    const VFPProdTable* table = detail::getTable(m_tables, thp_table_id);
    const double aqua_bhp_limit = rates_bhp_limit[Water];
    const double liquid_bhp_limit = rates_bhp_limit[Oil];
    const double vapour_bhp_limit = rates_bhp_limit[Gas];
    const double flo_bhp_limit = detail::getFlo(aqua_bhp_limit, liquid_bhp_limit, vapour_bhp_limit, table->getFloType() );

    const double aqua_bhp_middle = rates_bhp_middle[Water];
    const double liquid_bhp_middle = rates_bhp_middle[Oil];
    const double vapour_bhp_middle = rates_bhp_middle[Gas];
    const double flo_bhp_middle = detail::getFlo(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table->getFloType() );

    // we use the ratios based on the middle value of bhp_limit and bhp_safe_limit
    const double wfr = detail::getWFR(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table->getWFRType());
    const double gfr = detail::getGFR(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table->getGFRType());

    // we get the flo sampling points from the table,
    // then extend it with zero and rate under bhp_limit for extrapolation
    std::vector<double> flo_samples = table->getFloAxis();

    if (flo_samples[0] > 0.) {
        flo_samples.insert(flo_samples.begin(), 0.);
    }

    if (flo_samples.back() < std::abs(flo_bhp_limit)) {
        flo_samples.push_back(std::abs(flo_bhp_limit));
    }

    // kind of unncessarily following the tradation that producers should have negative rates
    // the key is here that it should be consistent with the function bhpwithflo
    for (double& value : flo_samples) {
        value = -value;
    }

    // get the bhp sampling values based on the flo sample values
    const std::vector<double> bhp_flo_samples = bhpwithflo(flo_samples, thp_table_id, wfr, gfr, thp_limit, alq, dp);

    std::vector<detail::RateBhpPair> ratebhp_samples;
    for (size_t i = 0; i < flo_samples.size(); ++i) {
        ratebhp_samples.push_back( detail::RateBhpPair{flo_samples[i], bhp_flo_samples[i]} );
    }

    const std::array<detail::RateBhpPair, 2> ratebhp_twopoints_ipr {detail::RateBhpPair{flo_bhp_middle, bhp_middle},
                                                                    detail::RateBhpPair{flo_bhp_limit, bhp_limit} };

    double obtain_bhp = 0.;
    const bool can_obtain_bhp_with_thp_limit = detail::findIntersectionForBhp(ratebhp_samples, ratebhp_twopoints_ipr, obtain_bhp);

    // \Note: assuming that negative BHP does not make sense
    if (can_obtain_bhp_with_thp_limit && obtain_bhp > 0.) {
        // getting too high bhp that might cause negative rates (rates in the undesired direction)
        if (obtain_bhp >= bhp_safe_limit) {
            const std::string msg (" We are getting a too high BHP value from the THP constraint, which may "
                                   " cause problems later ");
            OpmLog::info("TOO_HIGH_BHP_FOUND_THP_TARGET", msg);

            const std::string debug_msg = " obtain_bhp " + std::to_string(obtain_bhp)
                                        + " bhp_safe_limit " + std::to_string(bhp_safe_limit)
                                        + " thp limit " + std::to_string(thp_limit);
            OpmLog::debug(debug_msg);
        }
        return obtain_bhp;
    } else {
        OpmLog::warning("NO_BHP_FOUND_THP_TARGET", " we could not find a bhp value with thp target.");
        return -100.;
    }
}


}
