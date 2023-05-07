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
#include <opm/simulators/wells/WellPerforations.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/EvaluationFormat.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

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

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasOilRateInj(const std::vector<Value>& cq_s,
              PerforationRates& perf_rates,
              const Value& rv,
              const Value& rs,
              const Value& pressure,
              const Value& rvw,
              DeferredLogger& deferred_logger) const
{
    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    // TODO: the formulations here remain to be tested with cases with strong crossflow through production wells
    // s means standard condition, r means reservoir condition
    // q_os = q_or * b_o + rv * q_gr * b_g
    // q_gs = q_gr * b_g + rs * q_or * b_o
    // d = 1.0 - rs * rv
    // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
    // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)

    const double d = 1.0 - getValue(rv) * getValue(rs);

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, well_.name(),
                                          "gasOilRateInj",
                                          rs, rv, pressure));
    } else {
        // vaporized oil into gas
        // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
        perf_rates.vap_oil = getValue(rv) * (getValue(cq_s[gasCompIdx]) - getValue(rs) * getValue(cq_s[oilCompIdx])) / d;
        // dissolved of gas in oil
        // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
        perf_rates.dis_gas = getValue(rs) * (getValue(cq_s[oilCompIdx]) - getValue(rv) * getValue(cq_s[gasCompIdx])) / d;
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        // q_ws = q_wr * b_w + rvw * q_gr * b_g
        // q_wr = 1 / b_w * (q_ws - rvw * q_gr * b_g) = 1 / b_w * (q_ws - rvw * 1 / d  (q_gs - rs * q_os))
        // vaporized water in gas
        // rvw * q_gr * b_g = q_ws -q_wr *b_w = rvw * (q_gs -rs *q_os) / d
        perf_rates.vap_wat = getValue(rvw) * (getValue(cq_s[gasCompIdx]) - getValue(rs) * getValue(cq_s[oilCompIdx])) / d;
    }
}


template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasOilRateProd(std::vector<Value>& cq_s,
               PerforationRates& perf_rates,
               const Value& rv,
               const Value& rs,
               const Value& rvw) const
{
    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    const Value cq_sOil = cq_s[oilCompIdx];
    const Value cq_sGas = cq_s[gasCompIdx];
    const Value dis_gas = rs * cq_sOil;
    const Value vap_oil = rv * cq_sGas;

    cq_s[gasCompIdx] += dis_gas;
    cq_s[oilCompIdx] += vap_oil;

    // recording the perforation solution gas rate and solution oil rates
    if (well_.isProducer()) {
        perf_rates.dis_gas = getValue(dis_gas);
        perf_rates.vap_oil = getValue(vap_oil);
    }

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        const Value vap_wat = rvw * cq_sGas;
        cq_s[waterCompIdx] += vap_wat;
        if (well_.isProducer()) {
            perf_rates.vap_wat = getValue(vap_wat);
        }
    }
}

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasWaterRateInj(const std::vector<Value>& cq_s,
                PerforationRates& perf_rates,
                const Value& rvw,
                const Value& rsw) const
{
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    perf_rates.vap_wat = getValue(rvw) * getValue(cq_s[gasCompIdx]);
    perf_rates.dis_gas_in_water = getValue(rsw) * getValue(cq_s[waterCompIdx]);
}

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasWaterRateProd(std::vector<Value>& cq_s,
                 PerforationRates& perf_rates,
                 const Value& rvw,
                 const Value& rsw) const
{
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    const Value cq_sWat = cq_s[waterCompIdx];
    const Value cq_sGas = cq_s[gasCompIdx];
    const Value vap_wat = rvw * cq_sGas;
    const Value dis_gas_wat = rsw * cq_sWat;
    cq_s[waterCompIdx] += vap_wat;
    cq_s[gasCompIdx]   += dis_gas_wat;
    if (well_.isProducer()) {
        perf_rates.vap_wat = getValue(vap_wat);
        perf_rates.dis_gas_in_water = getValue(dis_gas_wat);
    }
}

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
disOilVapWatVolumeRatio(Value& volumeRatio,
                        const Value& rvw,
                        const Value& rsw,
                        const Value& pressure,
                        const std::vector<Value>& cmix_s,
                        const std::vector<Value>& b_perfcells_dense,
                        DeferredLogger& deferred_logger) const
{
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    // Incorporate RSW/RVW factors if both water and gas active
    const Value d = 1.0 - rvw * rsw;

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, well_.name(),
                                          "disOilVapWatVolumeRatio",
                                          rsw, rvw, pressure));
    }
    const Value tmp_wat = d > 0.0 ? (cmix_s[waterCompIdx] - rvw * cmix_s[gasCompIdx]) / d : cmix_s[waterCompIdx];
    volumeRatio += tmp_wat / b_perfcells_dense[waterCompIdx];

    const Value tmp_gas = d > 0.0 ? (cmix_s[gasCompIdx] - rsw * cmix_s[waterCompIdx]) / d : cmix_s[waterCompIdx];
    volumeRatio += tmp_gas / b_perfcells_dense[gasCompIdx];
}

template<class FluidSystem, class Indices, class Scalar, class Value>
void WellPerforations<FluidSystem,Indices,Scalar,Value>::
gasOilVolumeRatio(Value& volumeRatio,
                  const Value& rv,
                  const Value& rs,
                  const Value& pressure,
                  const std::vector<Value>& cmix_s,
                  const std::vector<Value>& b_perfcells_dense,
                  DeferredLogger& deferred_logger) const
{
    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    // Incorporate RS/RV factors if both oil and gas active
    const Value d = 1.0 - rv * rs;

    if (d <= 0.0) {
        deferred_logger.debug(dValueError(d, well_.name(),
                                          "gasOilVolumeRatio",
                                          rs, rv, pressure));
    }
    const Value tmp_oil = d > 0.0? (cmix_s[oilCompIdx] - rv * cmix_s[gasCompIdx]) / d : cmix_s[oilCompIdx];
    volumeRatio += tmp_oil / b_perfcells_dense[oilCompIdx];

    const Value tmp_gas =  d > 0.0? (cmix_s[gasCompIdx] - rs * cmix_s[oilCompIdx]) / d : cmix_s[gasCompIdx];
    volumeRatio += tmp_gas / b_perfcells_dense[gasCompIdx];
}

#define INSTANCE(Dim, ...) \
template class WellPerforations<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                                __VA_ARGS__, double, double>; \
template class WellPerforations<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, \
                                __VA_ARGS__, double, \
                                DenseAd::Evaluation<double,-1,Dim>>;

INSTANCE(4u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(5u, BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(9u, BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(6u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(7u, BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(8u, BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

INSTANCE(8u, BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(9u, BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(10u, BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(10u, BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)

} // namespace Opm
