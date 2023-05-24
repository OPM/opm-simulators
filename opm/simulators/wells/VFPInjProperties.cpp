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
#include <opm/simulators/wells/VFPInjProperties.hpp>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>

#include <limits>

namespace Opm {

double VFPInjProperties::bhp(int table_id,
                                 const double& aqua,
                                 const double& liquid,
                                 const double& vapour,
                                 const double& thp_arg) const {
    const VFPInjTable& table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg);
    return retval.value;
}

double VFPInjProperties::thp(int table_id,
                             const double& aqua,
                             const double& liquid,
                             const double& vapour,
                             const double& bhp_arg) const {
    const VFPInjTable& table = detail::getTable(m_tables, table_id);

    //Find interpolation variables
    const double flo = detail::getFlo(table, aqua, liquid, vapour);
    if (std::abs(flo) < std::numeric_limits<double>::epsilon()) {
        return 0.;
    }

    const std::vector<double> thp_array = table.getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    auto flo_i = detail::findInterpData(flo, table.getFloAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(table, flo_i, thp_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}

const VFPInjTable& VFPInjProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}

bool VFPInjProperties::hasTable(const int table_id) const {
    return detail::hasTable(m_tables, table_id);
}

void VFPInjProperties::addTable(const VFPInjTable& new_table) {
    this->m_tables.emplace( new_table.getTableNum(), new_table );
}

template <class EvalWell>
EvalWell VFPInjProperties::bhp(const int table_id,
                               const EvalWell& aqua,
                               const EvalWell& liquid,
                               const EvalWell& vapour,
                               const double& thp) const
{
    //Get the table
    const VFPInjTable& table = detail::getTable(m_tables, table_id);
    EvalWell bhp = 0.0 * aqua;

    //Find interpolation variables
    EvalWell flo = detail::getFlo(table, aqua, liquid, vapour);

    //First, find the values to interpolate between
    //Value of FLO is negative in OPM for producers, but positive in VFP table
    auto flo_i = detail::findInterpData(flo.value(), table.getFloAxis());
    auto thp_i = detail::findInterpData( thp, table.getTHPAxis()); // assume constant

    detail::VFPEvaluation bhp_val = detail::interpolate(table, flo_i, thp_i);

    bhp = bhp_val.dflo * flo;
    bhp.setValue(bhp_val.value); // thp is assumed constant i.e.
    return bhp;
}

#define INSTANCE(...) \
    template __VA_ARGS__ VFPInjProperties::bhp<__VA_ARGS__>(const int, \
                                                            const __VA_ARGS__&, \
                                                            const __VA_ARGS__&, \
                                                            const __VA_ARGS__&, \
                                                            const double&) const;

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
} //Namespace Opm
