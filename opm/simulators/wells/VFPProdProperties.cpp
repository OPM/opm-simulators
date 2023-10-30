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

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>

#include <cstddef>

namespace Opm {




double VFPProdProperties::thp(int table_id,
                              const double& aqua,
                              const double& liquid,
                              const double& vapour,
                              const double& bhp_arg,
                              const double& alq) const {
    const VFPProdTable& table = detail::getTable(m_tables, table_id);

    // Find interpolation variables.
    double flo = 0.0;
    double wfr = 0.0;
    double gfr = 0.0;
    if (aqua == 0.0 && liquid == 0.0 && vapour == 0.0) {
        // All zero, likely at initial state.
        // Set FLO variable to minimum to avoid extrapolation.
        // The water and gas fractions are kept at 0.0.
        flo = table.getFloAxis().front();
    } else {
        // The usual case.
        // Recall that production rate is negative in Opm, so switch the sign.
        flo = -detail::getFlo(table, aqua, liquid, vapour);
        wfr = detail::getWFR(table, aqua, liquid, vapour);
        gfr = detail::getGFR(table, aqua, liquid, vapour);
    }

    const std::vector<double> thp_array = table.getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small.
     */
    auto flo_i = detail::findInterpData( flo, table.getFloAxis());
    auto wfr_i = detail::findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table.getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table.getALQAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}


double VFPProdProperties::bhp(int table_id,
                              const double& aqua,
                              const double& liquid,
                              const double& vapour,
                              const double& thp_arg,
                              const double& alq,
                              const double& explicit_wfr,
                              const double& explicit_gfr,
                              const bool    use_expvfp, 
                              const double ipr_slope) const {
    const VFPProdTable& table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg, alq, explicit_wfr,explicit_gfr, use_expvfp);
    return retval.value;
}


const VFPProdTable& VFPProdProperties::getTable(const int table_id) const {
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
    const VFPProdTable& table = detail::getTable(m_tables, table_id);
    const auto thp_i = detail::findInterpData( thp, table.getTHPAxis()); // assume constant
    const auto wfr_i = detail::findInterpData( wfr, table.getWFRAxis());
    const auto gfr_i = detail::findInterpData( gfr, table.getGFRAxis());
    const auto alq_i = detail::findInterpData( alq, table.getALQAxis()); //assume constant

    std::vector<double> bhps(flos.size(), 0.);
    for (std::size_t i = 0; i < flos.size(); ++i) {
        // Value of FLO is negative in OPM for producers, but positive in VFP table
        const auto flo_i = detail::findInterpData(-flos[i], table.getFloAxis());
        const detail::VFPEvaluation bhp_val = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);

        // TODO: this kind of breaks the conventions for the functions here by putting dp within the function
        bhps[i] = bhp_val.value - dp;
    }

    return bhps;
}

double
VFPProdProperties::
minimumBHP(const int table_id,
           const double thp,
           const double wfr,
           const double gfr,
           const double alq) const
{
    // Get the table
    const VFPProdTable& table = detail::getTable(m_tables, table_id);
    const auto retval = detail::getMinimumBHPCoordinate(table, thp, wfr, gfr, alq);
    // returned pair is (flo, bhp)
    return retval.second;
}



void VFPProdProperties::addTable(const VFPProdTable& new_table) {
    this->m_tables.emplace( new_table.getTableNum(), new_table );
}

template <class EvalWell>
EvalWell VFPProdProperties::bhp(const int table_id,
                                const EvalWell& aqua,
                                const EvalWell& liquid,
                                const EvalWell& vapour,
                                const double& thp,
                                const double& alq,
                                const double& explicit_wfr,
                                const double& explicit_gfr,
                                const bool use_expvfp,
                                const double ipr_slope /*=0*/) const
{
    //Get the table
    const VFPProdTable& table = detail::getTable(m_tables, table_id);
    EvalWell bhp = 0.0 * aqua;

    //Find interpolation variables
    EvalWell flo = detail::getFlo(table, aqua, liquid, vapour);
    EvalWell wfr = detail::getWFR(table, aqua, liquid, vapour);
    EvalWell gfr = detail::getGFR(table, aqua, liquid, vapour);
    if (use_expvfp || -flo.value() < table.getFloAxis().front()) {
        wfr = explicit_wfr;
        gfr = explicit_gfr;
    }

    //First, find the values to interpolate between
    //Value of FLO is negative in OPM for producers, but positive in VFP table
    auto flo_i = detail::findInterpData(-flo.value(), table.getFloAxis());
    auto thp_i = detail::findInterpData( thp, table.getTHPAxis()); // assume constant
    auto wfr_i = detail::findInterpData( wfr.value(), table.getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr.value(), table.getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table.getALQAxis()); //assume constant

    detail::VFPEvaluation bhp_val = detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);

   
   //if (bhp_val.dflo < ipr_slope) {
   //     bhp = (bhp_val.dwfr * wfr) + (bhp_val.dgfr * gfr) - (std::max(0.0, bhp_val.dflo) * flo);
   //} else {
    bhp = (bhp_val.dwfr * wfr) + (bhp_val.dgfr * gfr) - (bhp_val.dflo* flo);
   //}

    bhp.setValue(bhp_val.value);
    return bhp;
}

#define INSTANCE(...) \
    template __VA_ARGS__ VFPProdProperties::bhp<__VA_ARGS__>(const int, \
                                                             const __VA_ARGS__&, const __VA_ARGS__&, const __VA_ARGS__&, \
                                                             const double&, const double&, const double&, const double&, const bool, const double) const;

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

}
