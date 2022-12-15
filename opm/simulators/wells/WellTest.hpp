/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELL_TEST_HEADER_INCLUDED
#define OPM_WELL_TEST_HEADER_INCLUDED

#include <limits>
#include <vector>

namespace Opm
{

class DeferredLogger;
struct PhaseUsage;
class SingleWellState;
class WellEconProductionLimits;
class WellInterfaceGeneric;
class WellTestState;

//! \brief Class for performing well tests.
class WellTest {
public:
    //! \brief Constructor sets reference to well.
    WellTest(const WellInterfaceGeneric& well) : well_(well) {}

    void updateWellTestStateEconomic(const SingleWellState& ws,
                                     const double simulation_time,
                                     const bool write_message_to_opmlog,
                                     WellTestState& well_test_state,
                                     DeferredLogger& deferred_logger) const;

    void updateWellTestStatePhysical(const double simulation_time,
                                     const bool write_message_to_opmlog,
                                     WellTestState& well_test_state,
                                     DeferredLogger& deferred_logger) const;

private:
    struct RatioLimitCheckReport {
        static constexpr int INVALIDCOMPLETION = std::numeric_limits<int>::max();
        bool ratio_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = 0.0;
    };

    void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          RatioLimitCheckReport& report) const;

    void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          RatioLimitCheckReport& report) const;

    void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                               const SingleWellState& ws,
                               RatioLimitCheckReport& report) const;

    template<class RatioFunc>
    bool checkMaxRatioLimitWell(const SingleWellState& ws,
                                const double max_ratio_limit,
                                const RatioFunc& ratioFunc) const;

    template<class RatioFunc>
    void checkMaxRatioLimitCompletions(const SingleWellState& ws,
                                       const double max_ratio_limit,
                                       const RatioFunc& ratioFunc,
                                       RatioLimitCheckReport& report) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                             const std::vector<double>& rates_or_potentials,
                             DeferredLogger& deferred_logger) const;

    RatioLimitCheckReport
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const SingleWellState& ws,
                         DeferredLogger& deferred_logger) const;


    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_TEST_HEADER_INCLUDED
