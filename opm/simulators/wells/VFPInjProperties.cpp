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

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>

#include <limits>

namespace Opm {

template<class Scalar>
Scalar VFPInjProperties<Scalar>::
bhp(const int    table_id,
    const Scalar aqua,
    const Scalar liquid,
    const Scalar vapour,
    const Scalar thp_arg) const
{
    const VFPInjTable& table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = VFPHelpers<Scalar>::bhp(table, aqua, liquid, vapour, thp_arg);
    return retval.value;
}

template<class Scalar>
Scalar VFPInjProperties<Scalar>::
thp(const int    table_id,
    const Scalar aqua,
    const Scalar liquid,
    const Scalar vapour,
    const Scalar bhp_arg) const
{
    const VFPInjTable& table = detail::getTable(m_tables, table_id);

    //Find interpolation variables
    const Scalar flo = detail::getFlo(table, aqua, liquid, vapour);
    if (std::abs(flo) < std::numeric_limits<Scalar>::epsilon()) {
        return 0.;
    }

    const auto thp_array = table.getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    const auto flo_i = VFPHelpers<Scalar>::findInterpData(flo, table.getFloAxis());
    std::vector<Scalar> bhp_array(nthp);
    for (int i = 0; i < nthp; ++i) {
        auto thp_i = VFPHelpers<Scalar>::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = VFPHelpers<Scalar>::interpolate(table, flo_i, thp_i).value;
    }

    return VFPHelpers<Scalar>::findTHP(bhp_array, thp_array, bhp_arg, /*find_largest*/ false);
}

template<class Scalar>
const VFPInjTable&
VFPInjProperties<Scalar>::getTable(const int table_id) const
{
    return detail::getTable(m_tables, table_id);
}

template<class Scalar>
bool VFPInjProperties<Scalar>::hasTable(const int table_id) const
{
    return detail::hasTable(m_tables, table_id);
}

template<class Scalar>
void VFPInjProperties<Scalar>::addTable(const VFPInjTable& new_table)
{
    this->m_tables.emplace( new_table.getTableNum(), new_table );
}

template<class Scalar>
template <class EvalWell>
EvalWell VFPInjProperties<Scalar>::bhp(const int       table_id,
                                       const EvalWell& aqua,
                                       const EvalWell& liquid,
                                       const EvalWell& vapour,
                                       const Scalar    thp) const
{
    //Get the table
    const VFPInjTable& table = detail::getTable(m_tables, table_id);
    EvalWell bhp = 0.0 * aqua;

    //Find interpolation variables
    EvalWell flo = detail::getFlo(table, aqua, liquid, vapour);

    //First, find the values to interpolate between
    //Value of FLO is negative in OPM for producers, but positive in VFP table
    const auto flo_i = VFPHelpers<Scalar>::findInterpData(flo.value(), table.getFloAxis());
    const auto thp_i = VFPHelpers<Scalar>::findInterpData(thp, table.getTHPAxis()); // assume constant

    detail::VFPEvaluation bhp_val = VFPHelpers<Scalar>::interpolate(table, flo_i, thp_i);

    bhp = bhp_val.dflo * flo;
    bhp.setValue(bhp_val.value); // thp is assumed constant i.e.
    return bhp;
}

#define INSTANTIATE(T,...)                           \
    template __VA_ARGS__                             \
        VFPInjProperties<T>::bhp(const int,          \
                                 const __VA_ARGS__&, \
                                 const __VA_ARGS__&, \
                                 const __VA_ARGS__&, \
                                 const T) const;

#define INSTANTIATE_TYPE(T)                        \
    template class VFPInjProperties<T>;            \
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

} //Namespace Opm
