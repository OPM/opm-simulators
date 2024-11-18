/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <config.h>
#include <opm/simulators/wells/RatioCalculator.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/EvaluationFormat.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <fmt/format.h>

namespace {
    template<class dValue, class Value>
    auto dValueError(const dValue& d,
                     const std::string& name,
                     const std::string& methodName,
                     const Value& Rs,
                     const Value& Rv,
                     const Value& pressure)
    {
        return fmt::format("Problematic d value {} obtained for well {}"
                           " during {} calculations with rs {}"
                           ", rv {} and pressure {}."
                           " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                           " for this connection.", d, name, methodName, Rs, Rv, pressure);
    }
}

namespace Opm {

template<class Value>
RatioCalculator<Value>::
RatioCalculator(unsigned gasCompIdx,
                unsigned oilCompIdx,
                unsigned waterCompIdx,
                std::string_view name)
    : gasComp_{gasCompIdx}
    , oilComp_(oilCompIdx)
    , waterComp_{waterCompIdx}
    , name_(name)
{
}

template<class Value>
void
RatioCalculator<Value>::
disOilVapWatVolumeRatio(Value& volumeRatio,
                        const Value& rvw,
                        const Value& rsw,
                        const Value& pressure,
                        const std::vector<Value>& cmix_s,
                        const std::vector<Value>& b_perfcells_dense,
                        DeferredLogger& deferred_logger) const
{
    // Incorporate RSW/RVW factors if both water and gas active
    const Value d = 1.0 - rvw * rsw;

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, name_,
                                          "disOilVapWatVolumeRatio",
                                          rsw, rvw, pressure));
    }
    const Value tmp_wat = d > 0.0 ? (cmix_s[waterComp_] - rvw * cmix_s[gasComp_]) / d
                                  : cmix_s[waterComp_];
    volumeRatio += tmp_wat / b_perfcells_dense[waterComp_];

    const Value tmp_gas =  d > 0.0 ? (cmix_s[gasComp_] - rsw * cmix_s[waterComp_]) / d
                                   : cmix_s[gasComp_];
    volumeRatio += tmp_gas / b_perfcells_dense[gasComp_];
}

template<class Value>
void
RatioCalculator<Value>::
gasOilVolumeRatio(Value& volumeRatio,
                  const Value& rv,
                  const Value& rs,
                  const Value& pressure,
                  const std::vector<Value>& cmix_s,
                  const std::vector<Value>& b_perfcells_dense,
                  DeferredLogger& deferred_logger) const
{
    // Incorporate RS/RV factors if both oil and gas active
    const Value d = 1.0 - rv * rs;

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, name_,
                                          "gasOilVolumeRatio",
                                          rs, rv, pressure));
    }
    const Value tmp_oil = d > 0.0 ? (cmix_s[oilComp_] - rv * cmix_s[gasComp_]) / d
                                  : cmix_s[oilComp_];
    volumeRatio += tmp_oil / b_perfcells_dense[oilComp_];

    const Value tmp_gas =  d > 0.0 ? (cmix_s[gasComp_] - rs * cmix_s[oilComp_]) / d
                                   : cmix_s[gasComp_];
    volumeRatio += tmp_gas / b_perfcells_dense[gasComp_];
}

#define INSTANTIATE_TYPE(T)                                          \
    template class RatioCalculator<T>;                               \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 4u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 5u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 6u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 7u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 8u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 9u>>;  \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 10u>>; \
    template class RatioCalculator<DenseAd::Evaluation<T, -1, 11u>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
