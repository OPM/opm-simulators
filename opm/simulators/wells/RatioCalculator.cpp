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
#include <opm/simulators/wells/PerforationData.hpp>

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
RatioCalculator(int gasCompIdx,
                int oilCompIdx,
                int waterCompIdx,
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
gasOilPerfRateInj(const std::vector<Value>& cq_s,
                  PerforationRates<Scalar>& perf_rates,
                  const Value& rv,
                  const Value& rs,
                  const Value& pressure,
                  const Value& rvw,
                  const bool waterActive,
                  DeferredLogger& deferred_logger) const
{
    // TODO: the formulations here remain to be tested with cases with strong crossflow through production wells
    // s means standard condition, r means reservoir condition
    // q_os = q_or * b_o + rv * q_gr * b_g
    // q_gs = q_gr * b_g + rs * q_or * b_o
    // d = 1.0 - rs * rv
    // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
    // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)

    const Scalar d = 1.0 - getValue(rv) * getValue(rs);

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, name_,
                                          "gasOilPerfRateInj",
                                          rs, rv, pressure));
    } else {
        // vaporized oil into gas
        // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
        perf_rates.vap_oil = getValue(rv) * (getValue(cq_s[gasComp_]) -
                                             getValue(rs) * getValue(cq_s[oilComp_])) / d;
        // dissolved of gas in oil
        // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
        perf_rates.dis_gas = getValue(rs) * (getValue(cq_s[oilComp_]) -
                                             getValue(rv) * getValue(cq_s[gasComp_])) / d;
    }

    if (waterActive) {
        // q_ws = q_wr * b_w + rvw * q_gr * b_g
        // q_wr = 1 / b_w * (q_ws - rvw * q_gr * b_g) = 1 / b_w * (q_ws - rvw * 1 / d  (q_gs - rs * q_os))
        // vaporized water in gas
        // rvw * q_gr * b_g = q_ws -q_wr *b_w = rvw * (q_gs -rs *q_os) / d
        perf_rates.vap_wat = getValue(rvw) * (getValue(cq_s[gasComp_]) -
                                              getValue(rs) * getValue(cq_s[oilComp_])) / d;
    }
}

template<class Value>
void
RatioCalculator<Value>::
gasOilPerfRateProd(std::vector<Value>& cq_s,
                   PerforationRates<Scalar>& perf_rates,
                   const Value& rv,
                   const Value& rs,
                   const Value& rvw,
                   const bool waterActive,
                   const bool isProducer) const
{
    const Value cq_sOil = cq_s[oilComp_];
    const Value cq_sGas = cq_s[gasComp_];
    const Value dis_gas = rs * cq_sOil;
    const Value vap_oil = rv * cq_sGas;

    cq_s[gasComp_] += dis_gas;
    cq_s[oilComp_] += vap_oil;

    // recording the perforation solution gas rate and solution oil rates
    if (isProducer) {
        perf_rates.dis_gas = getValue(dis_gas);
        perf_rates.vap_oil = getValue(vap_oil);
    }

    if (waterActive) {
        const Value vap_wat = rvw * cq_sGas;
        cq_s[waterComp_] += vap_wat;
        if (isProducer) {
            perf_rates.vap_wat = getValue(vap_wat);
        }
    }
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

template<class Value>
void
RatioCalculator<Value>::
gasWaterPerfRateInj(const std::vector<Value>& cq_s,
                     PerforationRates<Scalar>& perf_rates,
                     const Value& rvw,
                     const Value& rsw,
                     const Value& pressure,
                     DeferredLogger& deferred_logger) const
{
    const Scalar dw = 1.0 - getValue(rvw) * getValue(rsw);

    if (dw <= 0.0) {
        deferred_logger.debug(dValueError(dw, name_,
                                          "gasWaterPerfRateInj",
                                          rsw, rvw, pressure));
    } else {
        // vaporized water into gas
        // rvw * q_gr * b_g = rvw * (q_gs - rsw * q_ws) / dw
        perf_rates.vap_wat = getValue(rvw) * (getValue(cq_s[gasComp_]) -
                                              getValue(rsw) * getValue(cq_s[waterComp_])) / dw;
        // dissolved gas in water
        // rsw * q_wr * b_w = rsw * (q_ws - rvw * q_gs) / dw
        perf_rates.dis_gas_in_water = getValue(rsw) * (getValue(cq_s[waterComp_]) -
                                                       getValue(rvw) * getValue(cq_s[gasComp_])) / dw;
    }
}

template<class Value>
void
RatioCalculator<Value>::
gasWaterPerfRateProd(std::vector<Value>& cq_s,
                     PerforationRates<Scalar>& perf_rates,
                     const Value& rvw,
                     const Value& rsw,
                     const bool isProducer) const
{
    const Value cq_sWat = cq_s[waterComp_];
    const Value cq_sGas = cq_s[gasComp_];
    const Value vap_wat = rvw * cq_sGas;
    const Value dis_gas_wat = rsw * cq_sWat;
    cq_s[waterComp_] += vap_wat;
    cq_s[gasComp_]   += dis_gas_wat;
    if (isProducer) {
        perf_rates.vap_wat = getValue(vap_wat);
        perf_rates.dis_gas_in_water = getValue(dis_gas_wat);
    }
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
