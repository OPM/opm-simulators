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

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <cstddef>
#include <ctime>
#include <limits>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
template<typename Scalar, typename IndexTraits> class SingleWellState;
class UnitSystem;
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
                                     const UnitSystem& unit_system,
                                     const std::time_t start_time,
                                     DeferredLogger& deferred_logger) const;

    void updateWellTestStateCECON(const SingleWellState<Scalar, IndexTraits>& ws,
                                  const double simulation_time,
                                  const bool write_message_to_opmlog,
                                  WellTestState& well_test_state,
                                  const UnitSystem& unit_system,
                                  const std::time_t start_time,
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
        //! \brief Metadata describing the most-violating ratio limit, used to
        //!        build the human-readable workover message.
        std::string ratio_name{};
        UnitSystem::measure ratio_measure = UnitSystem::measure::identity;
        Scalar ratio_value = 0.0;
        Scalar ratio_limit = 0.0;
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

    //! \brief Check whether the well-level ratio exceeds \p max_ratio_limit;
    //!        \p well_ratio_value returns the computed ratio (reported by WECON).
    template<class RatioFunc>
    bool checkMaxRatioLimitWell(const SingleWellState<Scalar, IndexTraits>& ws,
                                const Scalar max_ratio_limit,
                                const RatioFunc& ratioFunc,
                                Scalar& well_ratio_value) const;

    template<class RatioFunc>
    void checkMaxRatioLimitCompletions(const SingleWellState<Scalar, IndexTraits>& ws,
                                       const Scalar max_ratio_limit,
                                       const Scalar well_ratio_value,
                                       const RatioFunc& ratioFunc,
                                       const std::string& ratio_name,
                                       const UnitSystem::measure ratio_measure,
                                       RatioLimitCheckReport& report) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                             const std::vector<Scalar>& rates_or_potentials,
                             DeferredLogger& deferred_logger) const;

    RatioLimitCheckReport
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const SingleWellState<Scalar, IndexTraits>& ws,
                         DeferredLogger& deferred_logger) const;

    //! \brief Describe a completion as "Completion N - block (i, j, k)" when it
    //!        owns a single connection, or just "Completion N" when it spans
    //!        several connections/blocks.
    std::string completionDescriptor(int complnum) const;

    //! \brief Apply the CON / +CON (close-connection) workover procedure.
    //!
    //! \param when    Human-readable "at time ... (date = ...)" clause shared
    //!                with the triggering economic-limit message.
    //! \param reason  Human-readable ratio-violation clause (e.g.
    //!                "water-gas ratio 1.0353e-06 SM3/SM3 exceeds the limit ...").
    //! \param ratio_subject  Owner of the violated ratio in the closing message,
    //!                inserted as "Because \p ratio_subject \p reason". WECON
    //!                reports a well-level ratio ("the well"); CECON reports the
    //!                completion's own ratio ("its").
    void closeOffendingCompletion(int offending_completion,
                                  bool close_connections_below,
                                  double simulation_time,
                                  bool write_message_to_opmlog,
                                  WellTestState& well_test_state,
                                  const std::string& when,
                                  const std::string& reason,
                                  const std::string& ratio_subject,
                                  DeferredLogger& deferred_logger) const;

    //! \brief A line of \p sep_length repetitions of \p sep_char used to frame
    //!        economic-limit workover messages (mirrors GroupEconomicLimitsChecker).
    std::string message_separator(const char sep_char = '*',
                                  const std::size_t sep_length = 110) const
    { return std::string(sep_length, sep_char); }

    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_TEST_HEADER_INCLUDED
