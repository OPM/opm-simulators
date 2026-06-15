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
#include <unordered_set>
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

    //! \param during_well_test  true when called from the WTEST re-open loop in
    //!        WellInterface::wellTesting(), which re-solves the well after every
    //!        completion closure; a CON/+CON workover then applies one workover
    //!        round per call (for +CON this can close multiple completions) and
    //!        defers further rounds to the caller's re-converged rates. false for
    //!        the regular timestep update, which is not re-solved between closures
    //!        and applies the whole workover here.
    void updateWellTestStateEconomic(const SingleWellState<Scalar, IndexTraits>& ws,
                                     const double simulation_time,
                                     const bool write_message_to_opmlog,
                                     const bool during_well_test,
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
        //! \brief Ratio value used when the denominator phase rate is
        //!        non-positive while the numerator is positive, i.e. the ratio
        //!        is effectively infinite. Always violates the limit; messages
        //!        must not present it as a physical value.
        static constexpr Scalar INFINITE_RATIO = 1.0e30;
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

    //! \brief Check one active max-ratio limit (water cut, GOR or WGR): if the
    //!        well-level ratio computed by \p ratioFunc exceeds \p max_ratio_limit,
    //!        mark the report violated and record the worst-offending completion
    //!        together with the ratio metadata (\p ratio_name, \p ratio_measure).
    template<class RatioFunc>
    void checkMaxRatioLimit(const SingleWellState<Scalar, IndexTraits>& ws,
                            const Scalar max_ratio_limit,
                            const RatioFunc& ratioFunc,
                            const std::unordered_set<int>& excluded_completions,
                            const std::string& ratio_name,
                            const UnitSystem::measure ratio_measure,
                            RatioLimitCheckReport& report) const;

    //! \brief Check whether the well-level ratio exceeds \p max_ratio_limit;
    //!        \p well_ratio_value returns the computed ratio (reported by WECON).
    //!
    //! With an empty \p excluded_completions the well rates are taken directly
    //! from the well state. During a CON / +CON workover event the completions
    //! closed so far still carry their converged rates, so the well rates are
    //! instead re-accumulated from the remaining completions.
    //!
    //! Note that the two bases are not identical: ws.surface_rates is the
    //! converged well-level rate, while the re-accumulated rates are sums of
    //! the perforation rates of the remaining completions, and the two can
    //! differ slightly (e.g. under cross-flow). Consequently, the ratio used to
    //! trigger the workover event and the ratios steering the subsequent
    //! cascade are on slightly different bases. This is considered acceptable
    //! because the cascade rates are frozen (not re-converged) approximations
    //! anyway; using the well-state rate for the trigger keeps the behavior of
    //! a plain (non-cascading) WECON check unchanged.
    template<class RatioFunc>
    bool checkMaxRatioLimitWell(const SingleWellState<Scalar, IndexTraits>& ws,
                                const Scalar max_ratio_limit,
                                const RatioFunc& ratioFunc,
                                const std::unordered_set<int>& excluded_completions,
                                Scalar& well_ratio_value) const;

    template<class RatioFunc>
    void checkMaxRatioLimitCompletions(const SingleWellState<Scalar, IndexTraits>& ws,
                                       const Scalar max_ratio_limit,
                                       const Scalar well_ratio_value,
                                       const RatioFunc& ratioFunc,
                                       const std::unordered_set<int>& excluded_completions,
                                       const std::string& ratio_name,
                                       const UnitSystem::measure ratio_measure,
                                       RatioLimitCheckReport& report) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                             const std::vector<Scalar>& rates_or_potentials,
                             DeferredLogger& deferred_logger) const;

    //! \brief Check all active ratio limits, ignoring \p excluded_completions
    //!        (completions already closed by the ongoing workover event).
    RatioLimitCheckReport
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const SingleWellState<Scalar, IndexTraits>& ws,
                         const std::unordered_set<int>& excluded_completions,
                         DeferredLogger& deferred_logger) const;

    //! \brief Describe a completion as "Completion N - block (i, j, k)" when it
    //!        owns a single connection, or just "Completion N" when it spans
    //!        several connections/blocks.
    std::string completionDescriptor(int complnum) const;

    //! \brief Apply one round of the CON / +CON (close-connection) workover:
    //!        close \p offending_completion (and, for +CON when
    //!        \p close_connections_below is set, every completion below it in the
    //!        wellbore), adding them to \p closed_this_event. \p when / \p reason
    //!        are the pre-formatted clauses of the logged closure message.
    //!        Returns true if this left the well with no open completions, i.e.
    //!        the well has been shut.
    bool closeOffendingCompletion(int offending_completion,
                                  bool close_connections_below,
                                  double simulation_time,
                                  bool write_message_to_opmlog,
                                  WellTestState& well_test_state,
                                  const std::string& when,
                                  const std::string& reason,
                                  std::unordered_set<int>& closed_this_event,
                                  DeferredLogger& deferred_logger) const;

    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_TEST_HEADER_INCLUDED
