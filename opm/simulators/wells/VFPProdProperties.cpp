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

template<class Scalar>
Scalar VFPProdProperties<Scalar>::
thp(const int    table_id,
    const Scalar aqua,
    const Scalar liquid,
    const Scalar vapour,
    const Scalar bhp_arg,
    const Scalar alq,
    const Scalar explicit_wfr,
    const Scalar explicit_gfr,
    const bool   use_expvfp) const
{
    const VFPProdTable& table = detail::getTable(m_tables, table_id);

    Scalar flo = detail::getFlo(table, aqua, liquid, vapour);
    Scalar wfr = detail::getWFR(table, aqua, liquid, vapour);
    Scalar gfr = detail::getGFR(table, aqua, liquid, vapour);
    if (use_expvfp || -flo < table.getFloAxis().front()) {
        wfr = explicit_wfr;
        gfr = explicit_gfr;
    }

    const std::vector<double> thp_array = table.getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small.
     */
    auto flo_i = VFPHelpers<Scalar>::findInterpData(-flo, table.getFloAxis());
    auto wfr_i = VFPHelpers<Scalar>::findInterpData( wfr, table.getWFRAxis());
    auto gfr_i = VFPHelpers<Scalar>::findInterpData( gfr, table.getGFRAxis());
    auto alq_i = VFPHelpers<Scalar>::findInterpData( alq, table.getALQAxis());
    std::vector<Scalar> bhp_array(nthp);
    for (int i = 0; i < nthp; ++i) {
        auto thp_i = VFPHelpers<Scalar>::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = VFPHelpers<Scalar>::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
    }

    return VFPHelpers<Scalar>::findTHP(bhp_array, thp_array, bhp_arg, /*find_largest*/ true);
}

template<class Scalar>
Scalar VFPProdProperties<Scalar>::
bhp(const int     table_id,
     const Scalar aqua,
     const Scalar liquid,
     const Scalar vapour,
     const Scalar thp_arg,
     const Scalar alq,
     const Scalar explicit_wfr,
     const Scalar explicit_gfr,
     const bool   use_expvfp) const
{
    const VFPProdTable& table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = VFPHelpers<Scalar>::bhp(table, aqua, liquid, vapour,
                                                           thp_arg, alq, explicit_wfr,
                                                           explicit_gfr, use_expvfp);
    return retval.value;
}

template<class Scalar>
const VFPProdTable&
VFPProdProperties<Scalar>::getTable(const int table_id) const
{
    return detail::getTable(m_tables, table_id);
}

template<class Scalar>
bool VFPProdProperties<Scalar>::hasTable(const int table_id) const
{
    return detail::hasTable(m_tables, table_id);
}

template<class Scalar>
std::vector<Scalar>
VFPProdProperties<Scalar>::
bhpwithflo(const std::vector<Scalar>& flos,
           const int table_id,
           const Scalar wfr,
           const Scalar gfr,
           const Scalar thp,
           const Scalar alq,
           const Scalar dp) const
{
    // Get the table
    const VFPProdTable& table = detail::getTable(m_tables, table_id);
    const auto thp_i = VFPHelpers<Scalar>::findInterpData( thp, table.getTHPAxis()); // assume constant
    const auto wfr_i = VFPHelpers<Scalar>::findInterpData( wfr, table.getWFRAxis());
    const auto gfr_i = VFPHelpers<Scalar>::findInterpData( gfr, table.getGFRAxis());
    const auto alq_i = VFPHelpers<Scalar>::findInterpData( alq, table.getALQAxis()); //assume constant

    std::vector<Scalar> bhps(flos.size(), 0.);
    for (std::size_t i = 0; i < flos.size(); ++i) {
        // Value of FLO is negative in OPM for producers, but positive in VFP table
        const auto flo_i = VFPHelpers<Scalar>::findInterpData(-flos[i], table.getFloAxis());
        const detail::VFPEvaluation bhp_val = VFPHelpers<Scalar>::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i);

        // TODO: this kind of breaks the conventions for the functions here by putting dp within the function
        bhps[i] = bhp_val.value - dp;
    }

    return bhps;
}

template<class Scalar>
Scalar VFPProdProperties<Scalar>::
minimumBHP(const int table_id,
           const Scalar thp,
           const Scalar wfr,
           const Scalar gfr,
           const Scalar alq) const
{
    // Get the table
    const VFPProdTable& table = detail::getTable(m_tables, table_id);
    const auto retval = VFPHelpers<Scalar>::getMinimumBHPCoordinate(table, thp, wfr, gfr, alq);
    // returned pair is (flo, bhp)
    return retval.second;
}

template<class Scalar>
void VFPProdProperties<Scalar>::addTable(const VFPProdTable& new_table)
{
    this->m_tables.emplace( new_table.getTableNum(), new_table );
}

template<class Scalar>
template <class EvalWell>
EvalWell VFPProdProperties<Scalar>::
bhp(const int       table_id,
    const EvalWell& aqua,
    const EvalWell& liquid,
    const EvalWell& vapour,
    const Scalar    thp,
    const Scalar    alq,
    const Scalar    explicit_wfr,
    const Scalar    explicit_gfr,
    const bool      use_expvfp) const
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
    auto flo_i = VFPHelpers<Scalar>::findInterpData(-flo.value(), table.getFloAxis());
    auto thp_i = VFPHelpers<Scalar>::findInterpData( thp, table.getTHPAxis()); // assume constant
    auto wfr_i = VFPHelpers<Scalar>::findInterpData( wfr.value(), table.getWFRAxis());
    auto gfr_i = VFPHelpers<Scalar>::findInterpData( gfr.value(), table.getGFRAxis());
    auto alq_i = VFPHelpers<Scalar>::findInterpData( alq, table.getALQAxis()); //assume constant

    detail::VFPEvaluation bhp_val = VFPHelpers<Scalar>::interpolate(table, flo_i, thp_i, wfr_i,
                                                                    gfr_i, alq_i);

    bhp = (bhp_val.dwfr * wfr) + (bhp_val.dgfr * gfr) - (bhp_val.dflo * flo);
    bhp.setValue(bhp_val.value);
    return bhp;
}

#define INSTANTIATE(T,...)                        \
    template __VA_ARGS__                          \
    VFPProdProperties<T>::bhp(const int,          \
                              const __VA_ARGS__&, \
                              const __VA_ARGS__&, \
                              const __VA_ARGS__&, \
                              const T ,           \
                              const T ,           \
                              const T ,           \
                              const T ,           \
                              const bool) const;

#define INSTANTIATE_TYPE(T)                        \
    template class VFPProdProperties<T>;           \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 4u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 5u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 6u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 7u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 8u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 9u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 10u>) \
    INSTANTIATE(T,DenseAd::Evaluation<T, -1, 11u>) \
    INSTANTIATE(T,DenseAd::Evaluation<T, 3, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 4, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 5, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 6, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 7, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 8, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 9, 0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T, 10, 0u>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
