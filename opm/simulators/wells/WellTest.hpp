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
template<typename Scalar, typename IndexTraits> class SingleWellState;
class WellEconProductionLimits;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
class WellTestState;

//! \brief Class for performing well tests.
template<typename Scalar, typename IndexTraits>
class WellTest {
public:
    //! \brief Constructor sets reference to well.
    explicit WellTest(const WellInterfaceGeneric<Scalar, IndexTraits>& well) : well_(well) {}

    void updateWellTestStateEconomic(const SingleWellState<Scalar, IndexTraits>& ws,
                                     const double simulation_time,
                                     const bool write_message_to_opmlog,
                                     WellTestState& well_test_state,
                                     bool zero_group_target,
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
        Scalar violation_extent = 0.0;
    };

    void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState<Scalar, IndexTraits>& ws,
                          RatioLimitCheckReport& report) const;

    void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState<Scalar, IndexTraits>& ws,
                          RatioLimitCheckReport& report) const;

    void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                               const SingleWellState<Scalar, IndexTraits>& ws,
                               RatioLimitCheckReport& report) const;

    template<class RatioFunc>
    bool checkMaxRatioLimitWell(const SingleWellState<Scalar, IndexTraits>& ws,
                                const Scalar max_ratio_limit,
                                const RatioFunc& ratioFunc) const;

    template<class RatioFunc>
    void checkMaxRatioLimitCompletions(const SingleWellState<Scalar, IndexTraits>& ws,
                                       const Scalar max_ratio_limit,
                                       const RatioFunc& ratioFunc,
                                       RatioLimitCheckReport& report) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                             const std::vector<Scalar>& rates_or_potentials,
                             DeferredLogger& deferred_logger) const;

    RatioLimitCheckReport
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const SingleWellState<Scalar, IndexTraits>& ws,
                         DeferredLogger& deferred_logger) const;


    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_TEST_HEADER_INCLUDED
