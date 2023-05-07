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

#ifndef OPM_WELL_PERFORATIONS_HPP
#define OPM_WELL_PERFORATIONS_HPP

#include <opm/simulators/wells/WellInterfaceIndices.hpp>

#include <vector>

namespace Opm {

class DeferredLogger;
class PerforationRates;

template<class FluidSystem, class Indices, class Scalar, class Value>
class WellPerforations {
public:
    using WellIfIndices = WellInterfaceIndices<FluidSystem,Indices,Scalar>;

    //! \brief Constructor sets reference to well.
    WellPerforations(const WellIfIndices& well) : well_(well) {}

    //! \brief Compute rate for an injecting perforation in system with gas and oil.
    void gasOilRateInj(const std::vector<Value>& cq_s,
                       PerforationRates& perf_rates,
                       const Value& rv,
                       const Value& rs,
                       const Value& pressure,
                       const Value& rvw,
                       DeferredLogger& deferred_logger) const;


    //! \brief Compute rate for a producing perforation in system with gas and oil.
    void gasOilRateProd(std::vector<Value>& cq_s,
                        PerforationRates& perf_rates,
                        const Value& rv,
                        const Value& rs,
                        const Value& rvw) const;

    //! \brief Compute rate for an injecting perforation in a system with gas and water.
    void gasWaterRateInj(const std::vector<Value>& cq_s,
                               PerforationRates& perf_rates,
                               const Value& rvw,
                               const Value& rsw) const;

    //! \brief Compute rate for a producing perforation in a system with gas and water.
    void gasWaterRateProd(std::vector<Value>& cq_s,
                          PerforationRates& perf_rates,
                          const Value& rvw,
                          const Value& rsw) const;

    //! \brief Compute volume ratio for a perforation in a system with dissolved oil and vaporized water.
    void disOilVapWatVolumeRatio(Value& volumeRatio,
                                 const Value& rvw,
                                 const Value& rsw,
                                 const Value& pressure,
                                 const std::vector<Value>& cmix_s,
                                 const std::vector<Value>& b_perfcells_dense,
                                 DeferredLogger& deferred_logger) const;

    //! \brief Compute volume ratio for perforation in a system with gas and oil.
    void gasOilVolumeRatio(Value& volumeRatio,
                           const Value& rv,
                           const Value& rs,
                           const Value& pressure,
                           const std::vector<Value>& cmix_s,
                           const std::vector<Value>& b_perfcells_dense,
                           DeferredLogger& deferred_logger) const;

private:
    const WellIfIndices& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_PERFORATIONS_HPP
