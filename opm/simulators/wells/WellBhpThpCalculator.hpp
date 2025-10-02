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


#ifndef OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED
#define OPM_WELL_BPH_THP_CALCULATOR_HEADER_INCLUDED

#include <functional>
#include <optional>
#include <vector>

namespace Opm {

class DeferredLogger;
class SummaryState;
class Well;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
template<typename Scalar, typename IndexTraits> class WellState;

//! \brief Class for computing BHP limits.
template<typename Scalar, typename IndexTraits>
class WellBhpThpCalculator {
public:
    //! \brief Constructor sets reference to well.
    explicit WellBhpThpCalculator(const WellInterfaceGeneric<Scalar, IndexTraits>& well) : well_(well) {}

    //! \brief Checks if well has THP constraints.
    bool wellHasTHPConstraints(const SummaryState& summaryState) const;

    //! \brief Get THP constraint for well.
    Scalar getTHPConstraint(const SummaryState& summaryState) const;

    //! \brief Obtain the most strict BHP from BHP limits.
    Scalar mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const;

    //! \brief Calculates THP from BHP.
    Scalar calculateThpFromBhp(const std::vector<Scalar>& rates,
                               const Scalar bhp,
                               const Scalar rho,
                               const std::optional<Scalar>& alq,
                               const Scalar thp_limit,
                               DeferredLogger& deferred_logger) const;

    //! \brief Compute BHP from THP limit for a producer.
    std::optional<Scalar>
    computeBhpAtThpLimitProd(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                             const SummaryState& summary_state,
                             const Scalar maxPerfPress,
                             const Scalar rho,
                             const Scalar alq_value,
                             const Scalar thp_limit,
                             DeferredLogger& deferred_logger) const;

    //! \brief Compute BHP from THP limit for an injector.
    std::optional<Scalar>
    computeBhpAtThpLimitInj(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                            const SummaryState& summary_state,
                            const Scalar rho,
                            const Scalar flo_rel_tol,
                            const int max_iteration,
                            const bool throwOnError,
                            DeferredLogger& deferred_logger) const;

    //! \brief Update THP.
    void updateThp(const Scalar rho,
                   const std::function<Scalar()>& alq_value,
                   WellState<Scalar, IndexTraits>& well_state,
                   const SummaryState& summary_state,
                   DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    EvalWell calculateBhpFromThp(const WellState<Scalar, IndexTraits>& well_state,
                                 const std::vector<EvalWell>& rates,
                                 const Well& well,
                                 const SummaryState& summaryState,
                                 const Scalar rho,
                                 DeferredLogger& deferred_logger) const;

    Scalar calculateMinimumBhpFromThp(const WellState<Scalar, IndexTraits>& well_state,
                                      const Well& well,
                                      const SummaryState& summaryState,
                                      const Scalar rho) const;

  bool isStableSolution(const WellState<Scalar, IndexTraits>& well_state,
                        const Well& well,
                        const std::vector<Scalar>& rates,
                        const SummaryState& summaryState) const;

  std::optional<Scalar>
  estimateStableBhp (const WellState<Scalar, IndexTraits>& well_state,
                    const Well& well,
                    const std::vector<Scalar>& rates,
                    const Scalar rho,
                    const SummaryState& summaryState) const;

  std::pair<Scalar, Scalar>
  getFloIPR(const WellState<Scalar, IndexTraits>& well_state,
            const Well& well,
            const SummaryState& summary_state) const;

  //! \brief Find limits using brute-force solver.
  static bool bruteForceBracketCommonTHP(const std::function<Scalar(const Scalar)>& eq,
                                const std::array<Scalar, 2>& range,
                                Scalar& low, Scalar& high,
                                std::optional<Scalar>& approximate_solution,
                                const Scalar& limit,
                                DeferredLogger& deferred_logger);

  //! \brief Find limits using brute-force solver.
  static bool bruteForceBracketCommonTHP(const std::function<Scalar(const Scalar)>& eq,
                                Scalar& min_thp, Scalar& max_thp);

private:
    //! \brief Compute BHP from THP limit for an injector - implementation.
    template<class ErrorPolicy>
    std::optional<Scalar>
    computeBhpAtThpLimitInjImpl(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                                const SummaryState& summary_state,
                                const Scalar rho,
                                const Scalar flo_rel_tol,
                                const int max_iteration,
                                DeferredLogger& deferred_logger) const;

    //! \brief Calculate max BHP.
    std::optional<Scalar>
    bhpMax(const std::function<Scalar(const Scalar)>& fflo,
           const Scalar bhp_limit,
           const Scalar maxPerfPress,
           const Scalar vfp_flo_front,
           DeferredLogger& deferred_logger) const;

    //! \brief Common code for finding BHP from THP limit for producers/injectors.
    std::optional<Scalar>
    computeBhpAtThpLimit(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                         const std::function<Scalar(const std::vector<Scalar>)>& fbhp,
                         const std::array<Scalar, 2>& range,
                         DeferredLogger& deferred_logger) const;

    //! \brief Get pressure adjustment to the bhp calculated from VFP table
    Scalar getVfpBhpAdjustment(const Scalar bph_tab, const Scalar thp_limit) const;

    //! \brief Find limits using bisection.
    bool bisectBracket(const std::function<Scalar(const Scalar)>& eq,
                       const std::array<Scalar, 2>& range,
                       Scalar& low, Scalar& high,
                       std::optional<Scalar>& approximate_solution,
                       DeferredLogger& deferred_logger) const;

 //! \brief Find limits using brute-force solver.
  static bool bruteForceBracket(const std::function<Scalar(const Scalar)>& eq,
                                const std::array<Scalar, 2>& range,
                                Scalar& low, Scalar& high,
                                DeferredLogger& deferred_logger);


    Scalar findThpFromBhpIteratively(const std::function<Scalar(const Scalar, const Scalar)>& thp_func,
                                     const Scalar bhp,
                                     const Scalar thp_limit,
                                     const Scalar dp,
                                     DeferredLogger& deferred_logger) const;

    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_BHP_THP_CALCULATOR_HEADER_INCLUDED
