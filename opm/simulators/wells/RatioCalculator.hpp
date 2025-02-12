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

#ifndef RATIO_CALCULATOR_HPP
#define RATIO_CALCULATOR_HPP

#include <opm/material/densead/Math.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace Opm {

class DeferredLogger;
template<class Scalar> struct PerforationRates;

template<class Value>
class RatioCalculator
{
public:
    using Scalar = decltype(getValue(Value{}));

    RatioCalculator(int gasCompIdx,
                    int oilCompIdx,
                    int waterCompIdx,
                    std::string_view name);

    void disOilVapWatVolumeRatio(Value& volumeRatio,
                                 const Value& rvw,
                                 const Value& rsw,
                                 const Value& pressure,
                                 const std::vector<Value>& cmix_s,
                                 const std::vector<Value>& b_perfcells_dense,
                                 DeferredLogger& deferred_logger) const;

    void gasOilPerfRateInj(const std::vector<Value>& cq_s,
                           PerforationRates<Scalar>& perf_rates,
                           const Value& rv,
                           const Value& rs,
                           const Value& pressure,
                           const Value& rvw,
                           const bool waterActive,
                           DeferredLogger& deferred_logger) const;

    void gasOilPerfRateProd(std::vector<Value>& cq_s,
                            PerforationRates<Scalar>& perf_rates,
                            const Value& rv,
                            const Value& rs,
                            const Value& rvw,
                            const bool waterActive,
                            const bool isProducer) const;

    void gasOilVolumeRatio(Value& volumeRatio,
                           const Value& rv,
                           const Value& rs,
                           const Value& pressure,
                           const std::vector<Value>& cmix_s,
                           const std::vector<Value>& b_perfcells_dense,
                           DeferredLogger& deferred_logger) const;

    void gasWaterPerfRateInj(const std::vector<Value>& cq_s,
                             PerforationRates<Scalar>& perf_rates,
                             const Value& rvw,
                             const Value& rsw,
                             const Value& pressure,
                             DeferredLogger& deferred_logger) const;

    void gasWaterPerfRateProd(std::vector<Value>& cq_s,
                              PerforationRates<Scalar>& perf_rates,
                              const Value& rvw,
                              const Value& rsw,
                              const bool isProducer) const;

private:
    int gasComp_;
    int oilComp_;
    int waterComp_;
    std::string name_;
};

} // namespace Opm

#endif // RATIO_CALCULATOR_HPP
